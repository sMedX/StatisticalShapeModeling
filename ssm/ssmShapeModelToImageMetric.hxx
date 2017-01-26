#pragma once

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

#include "ssmShapeModelToImageMetric.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TShapeModel, typename TImage >
ShapeModelToImageMetric<TShapeModel, TImage>::ShapeModelToImageMetric()
{
  m_ShapeModel = nullptr;
  m_PointsContainer = nullptr;
  m_Image = nullptr;
  m_SpatialTransform = nullptr;
  m_Interpolator = nullptr;
  m_ComputeGradient = true;
  m_GradientImage = nullptr;
  m_NumberOfSamplesCounted = 0;
  m_NumberOfThreads = 0;
  m_MaximalNumberOfThreads = 0;
}

/**
 * Set the parameters that define a unique transform
 */
template< typename TShapeModel, typename TImage >
void ShapeModelToImageMetric<TShapeModel, TImage>::SetTransformParameters(const ParametersType & parameters) const
{
  if( !m_Transform ) {
    itkExceptionMacro(<< "Transform has not been assigned");
  }
  m_Transform->SetParameters(parameters);
}

/**
 * PrintSelf
 */
template<typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::Initialize(void)
throw ( itk::ExceptionObject )
{
  if( !m_SpatialTransform ) {
    itkExceptionMacro(<< "Spatial transform is not present");
  }
  m_NumberOfSpatialParameters = m_SpatialTransform->GetNumberOfParameters();

  if( !m_Image ) {
    itkExceptionMacro(<< "Image is not present");
  }

  if (m_Image->GetSource()) {
    m_Image->GetSource()->Update();
  }

  if (!m_Interpolator) {
    typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;
    m_Interpolator = InterpolatorType::New();
  }
  m_Interpolator->SetInputImage(m_Image);

  if (!m_ShapeModel) {
    itkExceptionMacro(<< "ShapeModel is not present");
  }

  m_NumberOfComponents = m_ShapeModel->GetNumberOfPrincipalComponents();
  m_ShapeTransform.set_size(m_NumberOfComponents);

  m_NumberOfParameters = m_NumberOfComponents + m_NumberOfSpatialParameters;
  m_SpatialParameters.set_size(m_NumberOfParameters - m_NumberOfComponents);

  if ( m_ComputeGradient ) {
    double sigma = m_Image->GetSpacing().Get_vnl_vector().max_value();

    GradientImageFilterPointer gradient = GradientImageFilterType::New();
    gradient->SetInput(m_Image);
    gradient->SetSigma(sigma);
    gradient->SetNormalizeAcrossScale(true);
    gradient->Update();
    m_GradientImage = gradient->GetOutput();
  }

  m_PointsContainer = m_ShapeModel->GetRepresenter()->GetReference()->GetPoints();

  if (m_PointsContainer->Size() < 1) {
    itkExceptionMacro(<< "number of points is zero");
  }

  this->MultiThreadingInitialize();
}

template<typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::MultiThreadingInitialize()
{
  m_MaximalNumberOfThreads = omp_get_max_threads();

  if (m_NumberOfThreads == 0 || m_NumberOfThreads > m_MaximalNumberOfThreads) {
    m_NumberOfThreads = m_MaximalNumberOfThreads;
  }

  m_Threads.clear();

  for (size_t t = 0; t < m_NumberOfThreads; ++t) {
    PerThreadData thread;

    thread.m_Derivative = DerivativeType(m_NumberOfParameters);
    thread.m_Jacobian = TransformJacobianType(PointDimension, m_NumberOfParameters);
    thread.m_JacobianCache = TransformJacobianType(PointDimension, PointDimension);
    thread.m_ModelJacobian = TransformJacobianType(PointDimension, m_NumberOfComponents);
    thread.m_SpatialJacobian = TransformJacobianType(PointDimension, m_NumberOfSpatialParameters);

    m_Threads.push_back(thread);
  }

  size_t numberOfSamplesPerThread = itk::Math::Ceil<size_t, double>(m_PointsContainer->Size() / (double)m_NumberOfThreads);
  PointIteratorType iter = m_PointsContainer->Begin();

  for (size_t t = 0; t < m_NumberOfThreads; ++t) {
    m_Threads[t].m_Begin = iter;
    m_Threads[t].m_End = (iter += numberOfSamplesPerThread);
  }

  m_Threads[m_NumberOfThreads - 1].m_End = m_PointsContainer->End();
}

/**
* Get the match Measure
*/
template <typename TShapeModel, typename TImage>
typename ShapeModelToImageMetric<TShapeModel, TImage>::MeasureType 
ShapeModelToImageMetric<TShapeModel, TImage>::GetValue(const TransformParametersType & parameters ) const
{
  MeasureType value = itk::NumericTraits<MeasureType>::ZeroValue();
  DerivativeType derivative = DerivativeType(m_NumberOfParameters);
  derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());
  this->GetValueAndDerivative(parameters, value, derivative);
  return value;
}

/**
* Get the Derivative Measure
*/
template <typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  itkExceptionMacro(<< "not implemented");
}

/*
* Get both the match Measure and theDerivative Measure
*/
template <typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  for (size_t i = 0; i < m_NumberOfComponents; ++i) {
    m_ShapeTransform[i] = parameters[i];
  }

  for (size_t i = 0; i < m_NumberOfSpatialParameters; ++i) {
    m_SpatialParameters[i] = parameters[m_NumberOfComponents + i];
  }

  m_SpatialTransform->SetParameters(m_SpatialParameters);

  #pragma omp parallel num_threads(m_NumberOfThreads)
  {
    #pragma omp for
    for (int t = 0; t < m_NumberOfThreads; ++t) {
      GetValueAndDerivativeThreadProcessSample(m_Threads[t]);
    }
  }

  m_NumberOfSamplesCounted = 0;
  value = itk::NumericTraits<MeasureType>::ZeroValue();
  derivative = DerivativeType(m_NumberOfParameters);
  derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (size_t t = 0; t < m_NumberOfThreads; ++t ) {
    value += m_Threads[t].m_Value;
    m_NumberOfSamplesCounted += m_Threads[t].m_NumberOfSamplesCounted;

    for (size_t i = 0; i < m_NumberOfParameters; ++i) {
      derivative[i] += m_Threads[t].m_Derivative[i];
    }
  }

  if (!m_NumberOfSamplesCounted) {
    itkExceptionMacro(<< "All the points mapped to outside of the image");
  }
  else {
    value /= m_NumberOfSamplesCounted;
    for (size_t i = 0; i < m_NumberOfParameters; i++) {
      derivative[i] /= m_NumberOfSamplesCounted;
    }
  }

  this->CalculateValueAndDerivativePenalty(parameters, value, derivative);
}

template <typename TShapeModel, typename TImage>
inline void ShapeModelToImageMetric<TShapeModel, TImage>::GetValueAndDerivativeThreadProcessSample(PerThreadData & thread) const
{
  thread.m_NumberOfSamplesCounted = 0;
  thread.m_Value = itk::NumericTraits<MeasureType>::ZeroValue();
  thread.m_Derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (PointIteratorType iter = thread.m_Begin; iter != thread.m_End; ++iter) {
    InputPointType point = iter.Value();

    const OutputPointType modelTransformedPoint = m_ShapeModel->DrawSampleAtPoint(m_ShapeTransform, iter.Index());
    const OutputPointType transformedPoint = m_SpatialTransform->TransformPoint(modelTransformedPoint);

    if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
      thread.m_NumberOfSamplesCounted++;

      // compute the derivatives
      const typename TShapeModel::MatrixType & jacobian = m_ShapeModel->GetJacobian(iter.Index());
      for (size_t i = 0; i < PointDimension; ++i) {
        for (size_t k = 0; k < m_NumberOfComponents; ++k) {
          thread.m_ModelJacobian[i][k] = jacobian[i][k];
        }
      }

      m_SpatialTransform->ComputeJacobianWithRespectToParameters(modelTransformedPoint, thread.m_SpatialJacobian);
      thread.m_Jacobian.update(thread.m_SpatialJacobian, 0, m_NumberOfComponents);

      m_SpatialTransform->ComputeJacobianWithRespectToPosition(modelTransformedPoint, thread.m_JacobianCache);
      thread.m_Jacobian.update(thread.m_JacobianCache * thread.m_ModelJacobian, 0, 0);

      // get the gradient by NearestNeighboorInterpolation, which is equivalent to round up the point components.
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef typename itk::ContinuousIndex<CoordRepType, ImageDimension> MovingImageContinuousIndexType;
      MovingImageContinuousIndexType index;

      m_Image->TransformPhysicalPointToContinuousIndex(transformedPoint, index);
      typename ImageType::IndexType mappedIndex;
      mappedIndex.CopyWithRound(index);
      const GradientPixelType gradient = m_GradientImage->GetPixel(mappedIndex);

      // compute image value
      const RealType value = m_Interpolator->Evaluate(transformedPoint);
      thread.m_Value += std::pow(value, m_Degree);

      for (size_t i = 0; i < m_NumberOfParameters; ++i) {
        RealType sum = itk::NumericTraits<RealType>::ZeroValue();

        for (size_t d = 0; d < Self::PointDimension; ++d) {
          sum += m_Degree * std::pow(value, m_Degree - 1) * thread.m_Jacobian[d][i] * gradient[d];
        }

        thread.m_Derivative[i] += sum;
      }
    }
  }
}

/**
* Compute penalty
*/
template <typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::CalculateValuePenalty(const TransformParametersType & parameters, MeasureType & value) const
{
  MeasureType penaltyValue = 0;
  for (size_t n = 0; n < m_NumberOfComponents; ++n) {
    penaltyValue += parameters[n] * parameters[n];
  }
  value += penaltyValue * m_RegularizationParameter;
}

template<typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::CalculateDerivativePenalty(const TransformParametersType & parameters, DerivativeType  & derivative) const
{
  for (size_t n = 0; n < m_NumberOfComponents; ++n) {
    derivative[n] += 2 * parameters[n] * m_RegularizationParameter;
  }
}

template<typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::CalculateValueAndDerivativePenalty(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  MeasureType penaltyValue = 0;

  #pragma omp parallel reduction (+: penaltyValue) num_threads(m_NumberOfThreads)
  {
    # pragma omp for
    for (int n = 0; n < m_NumberOfComponents; ++n) {
      penaltyValue += parameters[n] * parameters[n];
      derivative[n] += 2 * parameters[n] * m_RegularizationParameter;
    }
  }

  value += penaltyValue * m_RegularizationParameter;
}

/**
 * PrintSelf 
 */
template< typename TShapeModel, typename TImage >
void ShapeModelToImageMetric<TShapeModel, TImage>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving image      " << m_Image.GetPointer()  << std::endl;
  os << indent << "Gradient image    " << m_GradientImage.GetPointer()   << std::endl;
  os << indent << "Spatial transform " << m_SpatialTransform.GetPointer()    << std::endl;
  os << indent << "Interpolator      " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Number of samples " << m_NumberOfSamplesCounted << std::endl;
  os << indent << "Compute gradient  " << m_ComputeGradient << std::endl;
  os << indent << "Number of threads " << m_NumberOfThreads << " / " << m_MaximalNumberOfThreads << std::endl;
}
}
