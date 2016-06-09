#ifndef __PointSetToImageMetrics_hxx
#define __PointSetToImageMetrics_hxx

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>

#include "PointSetToImageMetrics.h"

/**
 * Constructor
 */
template <typename TFixedPointSet, typename TMovingImage>
PointSetToImageMetrics<TFixedPointSet, TMovingImage>::PointSetToImageMetrics()
{
  typedef itk::LinearInterpolateImageFunction<TMovingImage, double> InterpolatorType;
  m_Interpolator = InterpolatorType::New();
}

/**
 * Get the match Measure
 */

template <typename TFixedPointSet, typename TMovingImage>
void PointSetToImageMetrics<TFixedPointSet, TMovingImage>::Compute()
{
  FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

  if( !fixedPointSet ) {
    itkExceptionMacro(<< "Fixed point set has not been assigned");
  }

  MovingImageConstPointer movingImage = this->GetMovingImage();

  if (!movingImage) {
    itkExceptionMacro(<< "Moving image has not been assigned");
  }

  m_Interpolator->SetInputImage(movingImage);

  PointIterator pointItr = fixedPointSet->GetPoints()->Begin();
  PointIterator pointEnd = fixedPointSet->GetPoints()->End();

  m_FixedValue = itk::NumericTraits<RealType>::Zero;

  m_MeanValue = itk::NumericTraits<MeasureType>::Zero;
  m_RMSEValue = itk::NumericTraits<MeasureType>::Zero;
  m_MaximalValue = itk::NumericTraits<MeasureType>::Zero;

  this->m_NumberOfPixelsCounted = 0;

  while( pointItr != pointEnd ) {
    InputPointType point = pointItr.Value();

    if (this->m_Interpolator->IsInsideBuffer(point)) {
      const RealType movingValue = this->m_Interpolator->Evaluate(point);
      const RealType error = std::abs(movingValue - m_FixedValue);

      m_MeanValue += error;
      m_RMSEValue += error * error;
      m_MaximalValue = std::max(error, m_MaximalValue);

      this->m_NumberOfPixelsCounted++;
    }

    ++pointItr;
  }

  if( !this->m_NumberOfPixelsCounted ) {
    itkExceptionMacro(<< "All the points mapped to outside of the moving image");
  }
  else {
    m_MeanValue = m_MeanValue / this->m_NumberOfPixelsCounted;
    m_RMSEValue = std::sqrt(m_RMSEValue/this->m_NumberOfPixelsCounted);
  }
}

template <typename TFixedPointSet, typename TMovingImage>
void PointSetToImageMetrics<TFixedPointSet, TMovingImage>::PrintReport(std::ostream& os) const
{
  std::string indent = "    ";

  os << "Metric values:" << std::endl;
  os << indent << "   Mean = " << m_MeanValue << std::endl;
  os << indent << "   RMSE = " << m_RMSEValue << std::endl;
  os << indent << "Maximal = " << m_MaximalValue << std::endl;
  os << std::endl;
}
#endif
