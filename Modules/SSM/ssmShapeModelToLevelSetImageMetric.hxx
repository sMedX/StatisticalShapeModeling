#pragma once

#include <itkLinearInterpolateImageFunction.h>

#include "ssmShapeModelToLevelSetImageMetric.h"

namespace ssm
{
/**
 * Constructor
 */
template <typename TShapeModel, typename TTransform, typename TImage>
ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::ShapeModelToLevelSetImageMetric()
{
  m_ShapeModel = nullptr;
  m_PointsContainer = nullptr;
  m_LevelSetImage = nullptr;
  m_Transform = nullptr;
  m_Interpolator = nullptr;
  m_ComputeGradient = true;
  m_GradientImage = nullptr;
  m_NumberOfSamplesCounted = 0;
  m_NumberOfThreads = 0;
  m_MaximalNumberOfThreads = 0;
  m_NumberOfEvaluations = 0;
  m_IsInitialized = false;
  m_NumberOfUsedPoints = 0;

  m_RegularizationParameter = 0.1;
  m_Degree = 2;

  m_Logger = itk::Logger::New();
}

/**
* Sparse points
*/
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::SparsePoints()
{
  if (m_NumberOfUsedPoints < 1) {
    m_PointsContainer = m_ShapeModel->GetRepresenter()->GetReference()->GetPoints();
  }
  else {
    auto points = m_ShapeModel->GetRepresenter()->GetReference()->GetPoints();
    m_PointsContainer = PointsContainerType::New();

    if (m_NumberOfUsedPoints > points->Size()) {
      itkExceptionMacro(<< "number of used points is too large")
    }

    size_t sparce = itk::Math::Ceil<size_t, float> ( (float)points->Size() / m_NumberOfUsedPoints);
    size_t count = 0;

    for (typename PointsContainerType::Iterator iter = points->Begin(); iter != points->End(); ++iter) {
      if (iter.Index() % sparce == 0) {
        m_PointsContainer->InsertElement(count++, iter.Value());
      }
    }
  }

  m_NumberOfUsedPoints = m_PointsContainer->Size();

  if (m_NumberOfUsedPoints < 1) {
    itkExceptionMacro(<< "number of points is too low.");
  }
}

/**
 * Initialize
 */
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::Initialize(void)
throw ( itk::ExceptionObject )
{
  // define numbers of used components
  if (!m_ShapeModel) {
    itkExceptionMacro(<< "ShapeModel is not initialized");
  }

  // check spatial transform
  if (!m_Transform) {
    itkExceptionMacro(<< "Spatial transform is not initialized");
  }

  // check image
  if( !m_LevelSetImage ) {
    itkExceptionMacro(<< "Level set image is not initialized");
  }
  if (m_LevelSetImage->GetSource()) {
    m_LevelSetImage->GetSource()->Update();
  }

  m_NumberOfParameters = m_Transform->GetNumberOfParameters();

  //if (!m_Interpolator) {
  //  typedef itk::LinearInterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;
  //  m_Interpolator = InterpolatorType::New();
  //}
  //m_Interpolator->SetInputImage(m_LevelSetImage);

  // compute gradient
  if ( m_ComputeGradient ) {
    auto sigma = m_LevelSetImage->GetSpacing().GetVnlVector().max_value();

    auto gradient = GradientImageFilterType::New();
    gradient->SetInput(m_LevelSetImage);
    gradient->SetSigma(sigma);
    gradient->SetNormalizeAcrossScale(true);
    gradient->Update();
    m_GradientImage = gradient->GetOutput();
  }

  // sparse points
  this->SparsePoints();

  // initialize multi threading data
  this->MultiThreadingInitialize();

  m_NumberOfEvaluations = 0;
  m_IsInitialized = true;
}

template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::MultiThreadingInitialize()
{
  m_MaximalNumberOfThreads = omp_get_max_threads();

  if (m_NumberOfThreads == 0 || m_NumberOfThreads > m_MaximalNumberOfThreads) {
    m_NumberOfThreads = m_MaximalNumberOfThreads;
  }

  m_Threads.clear();

  for (size_t t = 0; t < m_NumberOfThreads; ++t) {
    PerThreadData thread;

    thread.m_Derivative.SetSize(m_NumberOfParameters);
    thread.m_Jacobian.SetSize(PointDimension, m_NumberOfParameters);
    thread.m_ModelJacobian.SetSize(PointDimension, m_ShapeModel->GetNumberOfPrincipalComponents());
    thread.m_SpatialJacobian.SetSize(PointDimension, m_Transform->GetSpatialTransform()->GetNumberOfParameters());
    thread.m_JacobianCache.SetSize(PointDimension, PointDimension);

    m_Threads.push_back(thread);
  }

  auto numberOfSamplesPerThread = itk::Math::Ceil<size_t, double>(m_PointsContainer->Size() / (double) m_NumberOfThreads);
  auto iter = m_PointsContainer->Begin();

  for (size_t t = 0; t < m_NumberOfThreads; ++t) {
    m_Threads[t].m_Begin = iter;
    m_Threads[t].m_End = (iter += numberOfSamplesPerThread);
  }

  m_Threads[m_NumberOfThreads - 1].m_End = m_PointsContainer->End();
}

/**
* Get the match Measure
*/
template <typename TShapeModel, typename TTransform, typename TImage>
typename ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::MeasureType
ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::GetValue(const TransformParametersType & parameters ) const
{
  MeasureType value = itk::NumericTraits<MeasureType>::ZeroValue();
  DerivativeType derivative = DerivativeType(m_NumberOfParameters);
  this->GetValueAndDerivative(parameters, value, derivative);
  return value;
}

/**
* Get the Derivative Measure
*/
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::GetDerivative(const TransformParametersType & parameters, DerivativeType & derivative) const
{
  MeasureType value = itk::NumericTraits<MeasureType>::ZeroValue();
  this->GetValueAndDerivative(parameters, value, derivative);
}

/*
* Get both the match Measure and theDerivative Measure
*/
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  m_Transform->SetParameters(parameters);

  #pragma omp parallel num_threads(m_NumberOfThreads)
  {
    #pragma omp for
    for (int t = 0; t < m_NumberOfThreads; ++t) {
      this->GetValueAndDerivativeThreadProcessSample(m_Threads[t]);
    }
  }

  m_NumberOfSamplesCounted = 0;
  value = 0;

  for (size_t t = 0; t < m_NumberOfThreads; ++t ) {
    value += m_Threads[t].m_Value;
    m_NumberOfSamplesCounted += m_Threads[t].m_NumberOfSamplesCounted;
  }

  if (!m_NumberOfSamplesCounted) {
    itkExceptionMacro(<< "All the points mapped to outside of the image");
  }
  else {
    value /= m_NumberOfSamplesCounted;

    if (derivative.size() != m_NumberOfParameters) {
      derivative.set_size(m_NumberOfParameters);
    }

    #pragma omp parallel num_threads(m_NumberOfThreads)
    {
      #pragma omp for
      for (int i = 0; i < m_NumberOfParameters; ++i) {
        double sum = 0;

        for (size_t t = 0; t < m_NumberOfThreads; ++t) {
          sum += m_Threads[t].m_Derivative[i];
        }

        derivative[i] = sum / m_NumberOfSamplesCounted;
      }
    }
  }

  this->CalculateValueAndDerivativePenalty(parameters, value, derivative);
  m_NumberOfEvaluations++;
}

template <typename TShapeModel, typename TTransform, typename TImage>
inline void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::GetValueAndDerivativeThreadProcessSample(PerThreadData & thread) const
{
  thread.m_NumberOfSamplesCounted = 0;
  thread.m_Value = 0;
  thread.m_Derivative.Fill(0);

  typename ImageType::IndexType indexOfPoint;

  for (PointIteratorType iter = thread.m_Begin; iter != thread.m_End; ++iter) {
    const auto transformedPoint = m_Transform->TransformPoint(iter.Index());

    if (m_LevelSetImage->TransformPhysicalPointToIndex(transformedPoint, indexOfPoint)) {
      thread.m_NumberOfSamplesCounted++;

      // compute the derivatives
      m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(iter.Index(), thread.m_Jacobian, thread.m_ModelJacobian, thread.m_SpatialJacobian, thread.m_JacobianCache);

      // get the gradient by NearestNeighboorInterpolation, which is equivalent to round up the point components.
      const auto gradient = m_GradientImage->GetPixel(indexOfPoint);

      // compute image value
      // const double value = m_Interpolator->Evaluate(transformedPoint);
      const double value = m_LevelSetImage->GetPixel(indexOfPoint);
      thread.m_Value += std::pow(value, m_Degree);

      for (size_t i = 0; i < m_NumberOfParameters; ++i) {
        double sum = 0;

        for (size_t d = 0; d < PointDimension; ++d) {
          sum += thread.m_Jacobian[d][i] * gradient[d];
        }

        thread.m_Derivative[i] += m_Degree * std::pow(value, m_Degree - 1) * sum;
      }
    }
  }
}

/**
* Compute penalty
*/
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::CalculateValueAndDerivativePenalty(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const
{
  if (m_RegularizationParameter < itk::NumericTraits<double>::epsilon()) {
    return;
  }

  MeasureType penaltyValue = 0;

  //#pragma omp parallel reduction (+: penaltyValue) num_threads(m_NumberOfThreads)
  {
    //# pragma omp for
    for (int i = 0; i < m_Transform->GetNumberOfUsedComponents(); ++i) {
      penaltyValue += parameters[i] * parameters[i];
      derivative[i] += 2 * parameters[i] * m_RegularizationParameter;
    }
  }

  value += penaltyValue * m_RegularizationParameter;
}

/**
 * Print report 
 */
template <typename TShapeModel, typename TTransform, typename TImage>
void ShapeModelToLevelSetImageMetric<TShapeModel, TTransform, TImage>::PrintReport() const
{
  m_Message.str("");
  m_Message << this->GetNameOfClass() << std::endl;
  m_Message << std::endl;

  m_Message << "Number of used points " << m_NumberOfUsedPoints << ", number of counted samples " << m_NumberOfSamplesCounted << std::endl;
  m_Message << "Number of threads     " << m_NumberOfThreads << " / " << m_MaximalNumberOfThreads << std::endl;
  m_Message << "Compute gradient      " << m_ComputeGradient << std::endl;
  m_Message << "Regularization        " << m_RegularizationParameter << std::endl;
  m_Message << std::endl;

  m_Logger->Info(m_Message.str());
}
}
