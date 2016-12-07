#pragma once

#include "itkShapeModelToImageMetric.h"

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
}

/**
* Get the match Measure
*/
template <typename TShapeModel, typename TImage>
typename ShapeModelToImageMetric<TShapeModel, TImage>::MeasureType 
ShapeModelToImageMetric<TShapeModel, TImage>::GetValue(const TransformParametersType & parameters ) const
{
  m_Transform->SetParameters(parameters);

  m_NumberOfPixelsCounted = 0;
  PointIteratorType pointIter = m_PointsContainer->Begin();

  MeasureType value = itk::NumericTraits<MeasureType>::ZeroValue();

  for (; pointIter != m_PointsContainer->End(); ++pointIter) {
    InputPointType inputPoint;
    inputPoint.CastFrom(pointIter.Value());
    OutputPointType transformedPoint = m_Transform->TransformPoint(inputPoint);

    if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
      // compute image value
      const RealType imageValue = m_Interpolator->Evaluate(transformedPoint);
      value += imageValue * imageValue;
    
      m_NumberOfPixelsCounted++;
    }
  }

  if (!m_NumberOfPixelsCounted) {
    itkExceptionMacro(<< "All the points mapped to outside of the image");
  }
  else {
    value /= m_NumberOfPixelsCounted;
  }

  this->CalculateValuePenalty(parameters, value);

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

  m_NumberOfPixelsCounted = 0;
  value = itk::NumericTraits<MeasureType>::ZeroValue();

  derivative = DerivativeType(m_NumberOfParameters);
  derivative.Fill(itk::NumericTraits<typename DerivativeType::ValueType>::ZeroValue());

  PointIteratorType pointIter = m_PointsContainer->Begin();

  TransformJacobianType jacobian(TImage::ImageDimension, this->m_Transform->GetNumberOfParameters());
  TransformJacobianType jacobianCache(TImage::ImageDimension, TImage::ImageDimension);

  for (; pointIter != m_PointsContainer->End(); ++pointIter) {
    InputPointType inputPoint;
    inputPoint.CastFrom(pointIter.Value());
    OutputPointType transformedPoint = m_Transform->TransformPoint(inputPoint);

    if (this->m_Interpolator->IsInsideBuffer(transformedPoint)) {
      m_NumberOfPixelsCounted++;

      // compute the derivatives
      m_Transform->ComputeJacobianWithRespectToParametersCachedTemporaries(inputPoint, jacobian, jacobianCache);

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
      value += imageValue * imageValue;

      for (unsigned int par = 0; par < m_NumberOfParameters; par++) {
        RealType sum = itk::NumericTraits<RealType>::ZeroValue();

        for (unsigned int dim = 0; dim < Self::PointDimension; dim++) {
          sum += 2.0 *imageValue *jacobian(dim, par) * gradient[dim];
        }

        derivative[par] += sum;
      }
    }
  }

  if (!m_NumberOfPixelsCounted) {
    itkExceptionMacro(<< "All the points mapped to outside of the image");
  }
  else {
    value /= m_NumberOfPixelsCounted;

    for (unsigned int i = 0; i < m_NumberOfParameters; i++) {
      derivative[i] /= m_NumberOfPixelsCounted;
    }
  }

  this->CalculateValuePenalty(parameters, value);
  this->CalculateDerivativePenalty(parameters, derivative);
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

/**
 * PrintSelf 
 */
template< typename TShapeModel, typename TImage >
void ShapeModelToImageMetric<TShapeModel, TImage>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving Image: " << m_Image.GetPointer()  << std::endl;
  os << indent << "Gradient Image: " << m_GradientImage.GetPointer()   << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted << std::endl;
  os << indent << "Compute Gradient: " << m_ComputeGradient << std::endl;
}
}
