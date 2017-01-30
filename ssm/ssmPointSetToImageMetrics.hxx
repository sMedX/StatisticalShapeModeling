#pragma once

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSampleToHistogramFilter.h>
#include <itkListSample.h>
#include <itkHistogram.h>
#include <string>
#include <sstream>
#include <iostream>

#include "ssmPointSetToImageMetrics.h"

namespace ssm
{
  /**
   * Constructor
   */
  template <typename TFixedPointSet, typename TMovingImage>
  PointSetToImageMetrics<TFixedPointSet, TMovingImage>::PointSetToImageMetrics()
  {
    m_Interpolator = InterpolatorType::New();
    m_LevelOfQuantile = 0.95;
    m_HistogramSize = 1000;
  }

  /**
   * Get the match Measure
   */
  template <typename TFixedPointSet, typename TMovingImage>
  void PointSetToImageMetrics<TFixedPointSet, TMovingImage>::Compute()
  {
    FixedPointSetConstPointer fixedPointSet = this->GetFixedPointSet();

    if (!fixedPointSet) {
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

    typedef itk::Vector<double, 1> VectorType;
    typedef itk::Statistics::ListSample<VectorType> ListSampleType;
    ListSampleType::Pointer sample = ListSampleType::New();

    this->m_NumberOfPixelsCounted = 0;

    while (pointItr != pointEnd) {
      InputPointType point = pointItr.Value();

      if (this->m_Interpolator->IsInsideBuffer(point)) {
        const double movingValue = this->m_Interpolator->Evaluate(point);
        const double error = std::abs(movingValue - m_FixedValue);

        sample->PushBack(error);

        m_MeanValue += error;
        m_RMSEValue += error * error;
        m_MaximalValue = std::max(error, m_MaximalValue);

        this->m_NumberOfPixelsCounted++;
      }

      ++pointItr;
    }

    typedef itk::Statistics::Histogram<float, itk::Statistics::DenseFrequencyContainer2> HistogramType;
    typedef itk::Statistics::SampleToHistogramFilter<ListSampleType, HistogramType> SampleToHistogramFilterType;
    SampleToHistogramFilterType::Pointer sampleToHistogram = SampleToHistogramFilterType::New();
    sampleToHistogram->SetInput(sample);
    sampleToHistogram->SetAutoMinimumMaximum(true);

    SampleToHistogramFilterType::HistogramSizeType histogramSize(1);
    histogramSize.Fill(m_HistogramSize);

    sampleToHistogram->SetHistogramSize(histogramSize);

    try {
      sampleToHistogram->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      itkExceptionMacro(<< excep);
    }

    const HistogramType* histogram = sampleToHistogram->GetOutput();

    if (!this->m_NumberOfPixelsCounted) {
      itkExceptionMacro(<< "All the points mapped to outside of the moving image");
    }
    else {
      m_MeanValue = m_MeanValue / this->m_NumberOfPixelsCounted;
      m_RMSEValue = std::sqrt(m_RMSEValue / this->m_NumberOfPixelsCounted);
      m_QuantileValue = histogram->Quantile(0, m_LevelOfQuantile);
    }
  }

  template <typename TFixedPointSet, typename TMovingImage>
  void PointSetToImageMetrics<TFixedPointSet, TMovingImage>::PrintReport(std::ostream& os) const
  {
    std::string indent = "    ";

    os << "Metric values:" << std::endl;
    os << indent << "    Mean = " << m_MeanValue << std::endl;
    os << indent << "    RMSE = " << m_RMSEValue << std::endl;
    os << indent << "Quantile = " << m_QuantileValue << ", level = " << m_LevelOfQuantile << std::endl;
    os << indent << " Maximal = " << m_MaximalValue << std::endl;
    os << std::endl;
  }

  template <typename TFixedPointSet, typename TMovingImage>
  void PointSetToImageMetrics<TFixedPointSet, TMovingImage>::PrintReportToFile(const std::string & fileName, const std::string & datasetURI) const
  {
    std::string dlm = ";";

    std::string header = dlm;
    std::string scores = datasetURI + dlm;

    header += "Mean" + dlm;
    scores += std::to_string(m_MeanValue) + dlm;

    header += "RMSE" + dlm;
    scores += std::to_string(m_RMSEValue) + dlm;

    header += "Quantile " + std::to_string(m_LevelOfQuantile) + dlm;
    scores += std::to_string(m_QuantileValue) + dlm;

    header += "Maximal" + dlm;
    scores += std::to_string(m_MaximalValue) + dlm;

    header += dlm;
    scores += dlm;

    for (auto it = m_Info.begin(); it != m_Info.end(); ++it) {
      header += (*it).first + dlm;
      scores += (*it).second + dlm;
    }

    bool exist = boost::filesystem::exists(fileName);
    std::ofstream ofile;

    ofile.open(fileName, std::ofstream::out | std::ofstream::app);

    if (!exist) {
      ofile << header << std::endl;
    }

    ofile << scores << std::endl;
    ofile.close();
  }
}
