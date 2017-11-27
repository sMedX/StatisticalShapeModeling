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
  template< typename TPointSet, typename TImage >
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
    itkStaticConstMacro(Dimension, unsigned int, TPointSet::PointDimension);
    static_assert(TPointSet::PointDimension == TImage::ImageDimension, "Invalid dimensions of the input data.");

    /** Types transferred from the base class */
    typedef TPointSet                                             PointSetType;
    typedef typename PointSetType::PointType                      PointType;
    typedef TImage                                                ImageType;
    typedef typename PointSetType::PointsContainerConstIterator   PointIterator;
    typedef itk::LinearInterpolateImageFunction<TImage, double>   InterpolatorType;

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

      auto points = PointSetType::New();
      PointType point;

      for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
        for (size_t i = 0; i < Dimension; ++i) {
          point[i] = surface->GetPoints()->GetPoint(n)[i];
        }
        points->SetPoint(n, point);
      }

      this->SetPointSet(points);
    }

    template< typename MeshType>
    void SetPointSetAsMesh(const MeshType * mesh)
    {
      auto points = PointSetType::New();
      points->SetPoints(const_cast <typename PointSetType::PointsContainer*> (mesh->GetPoints()));
      this->SetPointSet(points);
    }

    /** Get/Set the fixed point set.  */
    itkSetConstObjectMacro(PointSet, PointSetType);
    itkGetConstObjectMacro(PointSet, PointSetType);

    /** Get/Set the moving image.  */
    itkSetConstObjectMacro(Image, ImageType);
    itkGetConstObjectMacro(Image, ImageType);

    /*Get/Set values to compute quantile. */
    itkSetMacro(LevelOfQuantile, double);
    itkGetMacro(LevelOfQuantile, double);

    itkSetMacro(HistogramSize, size_t);
    itkGetMacro(HistogramSize, size_t);

    /*Get metrics values. */
    itkGetMacro(MeanValue, double);
    itkGetMacro(RMSEValue, double);
    itkGetMacro(QuantileValue, double);
    itkGetMacro(MaximalValue, double);

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
      if (!m_PointSet) {
        itkExceptionMacro(<< "Fixed point set has not been assigned");
      }

      if (!m_Image) {
        itkExceptionMacro(<< "Moving image has not been assigned");
      }
      if (m_Image->GetSource()) {
        m_Image->GetSource()->Update();
      }

      m_Interpolator = InterpolatorType::New();
      m_Interpolator->SetInputImage(m_Image);

      m_MeanValue = 0;
      m_RMSEValue = 0;
      m_MaximalValue = 0;

      typedef itk::Vector<double, 1> VectorType;
      typedef itk::Statistics::ListSample<VectorType> ListSampleType;
      ListSampleType::Pointer sample = ListSampleType::New();

      m_NumberOfPixelsCounted = 0;

      for (PointIterator it = m_PointSet->GetPoints()->Begin(); it != m_PointSet->GetPoints()->End(); ++it) {
        PointType point = it.Value();

        if (m_Interpolator->IsInsideBuffer(point)) {
          const double value = std::abs(m_Interpolator->Evaluate(point));
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

      typedef typename itk::Statistics::Histogram<double, itk::Statistics::DenseFrequencyContainer2> HistogramType;
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
    size_t m_NumberOfPixelsCounted;

    double m_MeanValue;
    double m_RMSEValue;
    double m_QuantileValue;
    double m_MaximalValue;
    InfoType m_Info;

    typename PointSetType::ConstPointer m_PointSet;
    typename ImageType::ConstPointer m_Image;
  };
}
