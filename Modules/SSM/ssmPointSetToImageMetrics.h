#pragma once
#include <boost/filesystem.hpp>

#include <itkObject.h>
#include <itkPointSet.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkSampleToHistogramFilter.h>
#include <itkListSample.h>
#include <itkHistogram.h>

namespace ssm
{
  template< typename TFixedPointSet, typename TMovingImage >
  class PointSetToImageMetrics : public itk::Object
  {
  public:

    /** Standard class typedefs. */
    typedef PointSetToImageMetrics                      Self;
    typedef itk::Object                                 Superclass;
    typedef itk::SmartPointer< Self >                   Pointer;
    typedef itk::SmartPointer< const Self >             ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(PointSetToImageMetrics, Object);

    /** Constants for the image dimensions */
    itkStaticConstMacro(Dimension, unsigned int, TFixedPointSet::PointDimension);
    static_assert(TFixedPointSet::PointDimension == TMovingImage::ImageDimension, "Invalid dimensions of the input data.");

    /** Types transferred from the base class */
    typedef double                                       MeasureType;
    typedef TFixedPointSet                               FixedPointSetType;
    typedef typename FixedPointSetType::PointType        FixedPointType;
    typedef typename FixedPointSetType::ConstPointer     FixedPointSetConstPointer;
    typedef typename FixedPointSetType::PointsContainer  PointsContainerType;
    typedef TMovingImage                                 MovingImageType;
    typedef typename MovingImageType::ConstPointer       MovingImageConstPointer;

    typedef typename FixedPointSetType::PointsContainerConstIterator     PointIterator;
    typedef itk::LinearInterpolateImageFunction<TMovingImage, double>    InterpolatorType;

    /**  Type of the additional information. */
    typedef std::pair<std::string, std::string> PairType;
    typedef std::vector<PairType> InfoType;
    
    void SetInfo(const InfoType& info)  { m_Info = info; }

    /** Set the point set as vtkPolyData.  */
    template< typename PolyData>
    void SetPointSetAsPolyData(PolyData * surface)
    {
      if (Dimension != 3) {
        itkExceptionMacro(<<"For the method SetSurfaceAsPolyData dimension 3 is supproted.")
      }

      auto points = FixedPointSetType::New();
      FixedPointType point;

      for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
        for (size_t i = 0; i < Dimension; ++i) {
          point[i] = surface->GetPoints()->GetPoint(n)[i];
        }
        points->SetPoint(n, point);
      }

      this->SetFixedPointSet(points);
    }

    template< typename MeshType>
    void SetPointSetAsMesh(MeshType * mesh)
    {
      auto points = FixedPointSetType::New();
      points->SetPoints(mesh->GetPoints());
      this->SetFixedPointSet(points);
    }

    /** Get/Set the fixed point set.  */
    itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
    itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

    /** Get/Set the moving image.  */
    itkSetConstObjectMacro(MovingImage, MovingImageType);
    itkGetConstObjectMacro(MovingImage, MovingImageType);

    /*Get/Set values to compute quantile. */
    itkSetMacro(LevelOfQuantile, double);
    itkGetMacro(LevelOfQuantile, double);

    itkSetMacro(HistogramSize, size_t);
    itkGetMacro(HistogramSize, size_t);

    /*Get metrics values. */
    itkGetMacro(MeanValue, MeasureType);
    itkGetMacro(RMSEValue, MeasureType);
    itkGetMacro(QuantileValue, MeasureType);
    itkGetMacro(MaximalValue, MeasureType);

    void PrintReport(std::ostream& os) const
    {
      std::string indent = "    ";

      std::cout << "Information" << std::endl;
      for (const auto & pair : m_Info) {
        std::cout << pair.first << " " << pair.second << std::endl;
      }
      os << std::endl;

      os << "Metric values" << std::endl;
      os << indent << "    Mean = " << m_MeanValue << std::endl;
      os << indent << "    RMSE = " << m_RMSEValue << std::endl;
      os << indent << "Quantile = " << m_QuantileValue << ", level = " << m_LevelOfQuantile << std::endl;
      os << indent << " Maximal = " << m_MaximalValue << std::endl;
      os << std::endl;
    }

    void PrintReportToFile(const std::string & fileName, const std::string & datasetURI) const
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
      std::ofstream file(fileName, std::ofstream::out | std::ofstream::app);

      if (!exist) {
        file << header << std::endl;
      }

      file << scores << std::endl;
      file.close();
    }

    /**  Compute values. */
    void Compute()
    {
      if (!m_FixedPointSet) {
        itkExceptionMacro(<< "Fixed point set has not been assigned");
      }

      if (!m_MovingImage) {
        itkExceptionMacro(<< "Moving image has not been assigned");
      }
      if (m_MovingImage->GetSource()) {
        m_MovingImage->GetSource()->Update();
      }

      m_Interpolator = InterpolatorType::New();
      m_Interpolator->SetInputImage(m_MovingImage);

      m_MeanValue = itk::NumericTraits<MeasureType>::Zero;
      m_RMSEValue = itk::NumericTraits<MeasureType>::Zero;
      m_MaximalValue = itk::NumericTraits<MeasureType>::Zero;

      typedef itk::Vector<MeasureType, 1> VectorType;
      typedef itk::Statistics::ListSample<VectorType> ListSampleType;
      ListSampleType::Pointer sample = ListSampleType::New();

      m_NumberOfPixelsCounted = 0;

      for (PointIterator it = m_FixedPointSet->GetPoints()->Begin(); it != m_FixedPointSet->GetPoints()->End(); ++it) {
        FixedPointType point = it.Value();

        if (m_Interpolator->IsInsideBuffer(point)) {
          const MeasureType value = std::abs(m_Interpolator->Evaluate(point));
          sample->PushBack(value);

          m_MeanValue += value;
          m_RMSEValue += value * value;
          m_MaximalValue = std::max(value, m_MaximalValue);

          m_NumberOfPixelsCounted++;
        }
      }

      if (!m_NumberOfPixelsCounted) {
        itkExceptionMacro(<< "All the points mapped to outside of the moving image");
      }

      typedef typename itk::Statistics::Histogram<MeasureType, itk::Statistics::DenseFrequencyContainer2> HistogramType;
      typename HistogramType::SizeType size(1);
      size.Fill(m_HistogramSize);

      typedef itk::Statistics::SampleToHistogramFilter<ListSampleType, HistogramType> SampleToHistogramFilterType;
      SampleToHistogramFilterType::Pointer sampleToHistogram = SampleToHistogramFilterType::New();
      sampleToHistogram->SetInput(sample);
      sampleToHistogram->SetAutoMinimumMaximum(true);
      sampleToHistogram->SetHistogramSize(size);
      try {
        sampleToHistogram->Update();
      }
      catch (itk::ExceptionObject& excep) {
        itkExceptionMacro(<< excep);
      }

      m_MeanValue = m_MeanValue / m_NumberOfPixelsCounted;
      m_RMSEValue = std::sqrt(m_RMSEValue / m_NumberOfPixelsCounted);
      m_QuantileValue = sampleToHistogram->GetOutput()->Quantile(0, m_LevelOfQuantile);
    }

  protected:
    PointSetToImageMetrics() {}
    virtual ~PointSetToImageMetrics() {}

  private:
    PointSetToImageMetrics(const Self &);  // purposely not implemented
    void operator=(const Self &);          // purposely not implemented

    typename InterpolatorType::Pointer m_Interpolator;

    double m_LevelOfQuantile = 0.95;
    size_t m_HistogramSize = 1000;
    MeasureType m_MeanValue;
    MeasureType m_RMSEValue;
    MeasureType m_QuantileValue;
    MeasureType m_MaximalValue;
    InfoType m_Info;
    size_t m_NumberOfPixelsCounted;

    typename FixedPointSetType::ConstPointer m_FixedPointSet;
    typename MovingImageType::ConstPointer m_MovingImage;
  };
}
