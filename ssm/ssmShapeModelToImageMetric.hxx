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
  m_Transform = nullptr;
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
  if( !m_Transform ) {
    itkExceptionMacro(<< "Transform is not present");
  }

  m_NumberOfParameters = this->GetNumberOfParameters();

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

  if ( m_ComputeGradient ) {
    GradientImageFilterPointer gradient = GradientImageFilterType::New();
    gradient->SetInput(m_Image);

    const typename ImageType::SpacingType & spacing = m_Image->GetSpacing();
    double maximumSpacing = 0.0;
    for ( unsigned int i = 0; i < ImageDimension; i++ ) {
      if ( spacing[i] > maximumSpacing ) {
        maximumSpacing = spacing[i];
      }
    }
    gradient->SetSigma(maximumSpacing);
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

  for (unsigned int t = 0; t < m_NumberOfThreads; ++t) {
    PerThreadData thread;

    thread.m_Derivative = DerivativeType(m_NumberOfParameters);
    thread.m_Jacobian = TransformJacobianType(TImage::ImageDimension, m_Transform->GetNumberOfParameters());
    thread.m_JacobianCache = TransformJacobianType(TImage::ImageDimension, TImage::ImageDimension);

    m_Threads.push_back(thread);
  }

  NumberOfSamplesPerThread = itk::Math::Ceil<unsigned int, double>(m_PointsContainer->Size() / (double)m_NumberOfThreads);
  PointIteratorType iter = m_PointsContainer->Begin();

  for (unsigned int t = 0; t < m_NumberOfThreads; ++t) {
    m_Threads[t].m_Begin = iter;
    m_Threads[t].m_End = (iter += NumberOfSamplesPerThread);
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
  m_Transform->SetParameters(parameters);

  #pragma omp parallel num_threads(m_NumberOfThreads)
  {
    #pragma omp for
    for (int t = 0; t < m_NumberOfThreads; ++t) {
      GetValueAndDerivativeThreadProcessSample(m_Threads[t], parameters, value, derivative);
    }
  }

  m_NumberOfSamplesCounted = 0;
  value = itk::NumericTraits<MeasureType>::ZeroValue();
  derivative = DerivativeType(m_NumberOfParameters);
  derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (unsigned int t = 0; t < m_NumberOfThreads; ++t ) {
    value += m_Threads[t].m_Value;
    m_NumberOfSamplesCounted += m_Threads[t].m_NumberOfSamplesCounted;

    for (unsigned int i = 0; i < m_NumberOfParameters; ++i) {
      derivative[i] += m_Threads[t].m_Derivative[i];
    }
  }

  if (!m_NumberOfSamplesCounted) {
    itkExceptionMacro(<< "All the points mapped to outside of the image");
  }
  else {
    value /= m_NumberOfSamplesCounted;
    for (unsigned int i = 0; i < m_NumberOfParameters; i++) {
      derivative[i] /= m_NumberOfSamplesCounted;
    }
  }

  this->CalculatePenalty(parameters, value, derivative);
}

template <typename TShapeModel, typename TImage>
inline void ShapeModelToImageMetric<TShapeModel, TImage>::GetValueAndDerivativeThreadProcessSample(PerThreadData & thread, const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  thread.m_NumberOfSamplesCounted = 0;
  thread.m_Value = itk::NumericTraits<MeasureType>::ZeroValue();
  thread.m_Derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  for (PointIteratorType iter = thread.m_Begin; iter != thread.m_End; ++iter) {
    InputPointType point = iter.Value();
    OutputPointType transformedPoint = m_Transform->TransformPoint(point);

    if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
      thread.m_NumberOfSamplesCounted++;

      // compute the derivatives
      m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(point, thread.m_Jacobian, thread.m_JacobianCache);

      // get the gradient by NearestNeighboorInterpolation, which is equivalent to round up the point components.
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef itk::ContinuousIndex<CoordRepType, ImageType::ImageDimension> MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      m_Image->TransformPhysicalPointToContinuousIndex(transformedPoint, tempIndex);
      typename ImageType::IndexType mappedIndex;
      mappedIndex.CopyWithRound(tempIndex);
      const GradientPixelType gradient = m_GradientImage->GetPixel(mappedIndex);

      // compute image value
      const RealType imageValue = m_Interpolator->Evaluate(transformedPoint);
      thread.m_Value += std::pow(imageValue, m_Degree);

      for (unsigned int i = 0; i < m_NumberOfParameters; i++) {
        RealType sum = itk::NumericTraits<RealType>::ZeroValue();

        for (unsigned int d = 0; d < Self::PointDimension; d++) {
          sum += m_Degree * std::pow(imageValue, m_Degree - 1) * thread.m_Jacobian(d, i) * gradient[d];
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
  for (unsigned int n = 0; n < m_NumberOfComponents; ++n) {
    penaltyValue += parameters[n] * parameters[n];
  }
  value += penaltyValue * m_RegularizationParameter;
}

template<typename TShapeModel, typename TImage>
void ShapeModelToImageMetric<TShapeModel, TImage>::CalculateDerivativePenalty(const TransformParametersType & parameters, DerivativeType  & derivative) const
{
  for (unsigned int n = 0; n < m_NumberOfComponents; ++n) {
    derivative[n] += 2 * parameters[n] * m_RegularizationParameter;
  }
}

template<typename TShapeModel, typename TImage>
inline void ShapeModelToImageMetric<TShapeModel, TImage>::CalculatePenalty(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
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
  os << indent << "Transform         " << m_Transform.GetPointer()    << std::endl;
  os << indent << "Interpolator      " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Number of samples " << m_NumberOfSamplesCounted << std::endl;
  os << indent << "Compute gradient  " << m_ComputeGradient << std::endl;
  os << indent << "Number of threads " << m_NumberOfThreads << " / " << m_MaximalNumberOfThreads << std::endl;
}
}
