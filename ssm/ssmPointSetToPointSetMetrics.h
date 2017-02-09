#pragma once
#include <boost/filesystem.hpp>

#include <itkSingleValuedCostFunction.h>
#include <itkKdTreeGenerator.h>
#include <itkListSample.h>
#include <itkPointSet.h>

namespace ssm
{
  template< typename TFixedPointSet, typename TMovingPointSet = TFixedPointSet >
  class PointSetToPointSetMetrics : public itk::Object
  {
  public:
    /** Standard class typedefs. */
    typedef PointSetToPointSetMetrics               Self;
    typedef itk::Object                             Superclass;
    typedef itk::SmartPointer< Self >               Pointer;
    typedef itk::SmartPointer< const Self >         ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(PointSetToPointSetMetrics, itk::Object);

    /** Types transferred from the base class */
    typedef typename double                                MeasureType;
    typedef typename TFixedPointSet                        FixedPointSetType;
    typedef typename TMovingPointSet                       MovingPointSetType;
    typedef typename FixedPointSetType::ConstPointer       FixedPointSetConstPointer;
    typedef typename MovingPointSetType::ConstPointer      MovingPointSetConstPointer;

    /**  Type of the parameters. */
    typedef itk::Statistics::ListSample<typename MovingPointSetType::PointType> ListSampleType;
    typedef itk::Statistics::KdTreeGenerator<ListSampleType> TreeGeneratorType;
    typedef typename TreeGeneratorType::KdTreeType TreeType;
    typedef typename TreeType::InstanceIdentifierVectorType NeighborhoodType;

    typedef std::pair<std::string, std::string> PairType;
    typedef std::vector<PairType> InfoType;

    /** Get/Set the Fixed Image.  */
    itkSetConstObjectMacro(FixedPointSet, FixedPointSetType);
    itkGetConstObjectMacro(FixedPointSet, FixedPointSetType);

    /** Get/Set the Moving Image.  */
    itkSetConstObjectMacro(MovingPointSet, MovingPointSetType);
    itkGetConstObjectMacro(MovingPointSet, MovingPointSetType);

    /* Compute symmetric metrics*/
    itkSetMacro(Symmetric, bool);
    itkGetMacro(Symmetric, bool);

    /*Get/Set values to compute quantile. */
    itkSetMacro(LevelOfQuantile, double);
    itkGetMacro(LevelOfQuantile, double);

    itkSetMacro(HistogramSize, size_t);
    itkGetMacro(HistogramSize, size_t);

    itkSetMacro(Info, InfoType);

    /*Get metrics values. */
    itkGetMacro(MeanValue, MeasureType);
    itkGetMacro(RMSEValue, MeasureType);
    itkGetMacro(QuantileValue, MeasureType);
    itkGetMacro(MaximalValue, MeasureType);

    void PrintReport(std::ostream& os) const
    {
      std::string indent = "    ";
      os << "Metric values:" << std::endl;
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

    /** Compute metrics. */
    void Compute()
    {
      std::vector<MeasureType> fixedToMovingMetrics;
      this->ComputeMetrics(fixedToMovingMetrics, m_FixedPointSet, m_MovingPointSet);

      m_MeanValue = fixedToMovingMetrics[0];
      m_RMSEValue = fixedToMovingMetrics[1];
      m_QuantileValue = fixedToMovingMetrics[2];
      m_MaximalValue = fixedToMovingMetrics[3];

      if (m_Symmetric) {
        std::vector<MeasureType> movingToFixedMetrics;
        this->ComputeMetrics(movingToFixedMetrics, m_MovingPointSet, m_FixedPointSet);

        m_MeanValue = (fixedToMovingMetrics[0] + movingToFixedMetrics[0])/2;
        m_RMSEValue = (fixedToMovingMetrics[1] + movingToFixedMetrics[1])/2;
        m_QuantileValue = (fixedToMovingMetrics[2] + movingToFixedMetrics[2])/2;
        m_MaximalValue = (fixedToMovingMetrics[3] + movingToFixedMetrics[3])/2;
      }
    }


  protected:
    PointSetToPointSetMetrics() {}
    virtual ~PointSetToPointSetMetrics() {}

  private:
    PointSetToPointSetMetrics(const Self &);
    void operator=(const Self &);

    FixedPointSetConstPointer m_FixedPointSet;
    MovingPointSetConstPointer m_MovingPointSet;

    double m_Radius = 10;
    size_t m_BucketSize = 16;
    size_t m_HistogramSize = 1000;
    double m_LevelOfQuantile = 0.95;

    bool m_Symmetric = true;
    MeasureType m_MeanValue;
    MeasureType m_RMSEValue;
    MeasureType m_QuantileValue;
    MeasureType m_MaximalValue;
    InfoType m_Info;

    typename ListSampleType::Pointer m_ListOfPoints;
    typename TreeGeneratorType::Pointer m_TreeGenerator;
    typename TreeType::ConstPointer m_Tree;

    void ComputeMetrics(std::vector<MeasureType> & metrics, typename FixedPointSetType::ConstPointer fixedPointSet, typename MovingPointSetType::ConstPointer movingPointSet)
    {
      this->InitializeKdTree(movingPointSet);

      MeasureType mean = itk::NumericTraits<MeasureType>::Zero;
      MeasureType rmse = itk::NumericTraits<MeasureType>::Zero;
      MeasureType maximal = itk::NumericTraits<MeasureType>::Zero;

      typename FixedPointSetType::PointsContainer::ConstPointer fixedContainer = fixedPointSet->GetPoints();
      typename MovingPointSetType::PointsContainer::ConstPointer  movingContainer = movingPointSet->GetPoints();

      typedef itk::Vector<MeasureType, 1> VectorType;
      typedef itk::Statistics::ListSample<VectorType> ListSampleType;
      ListSampleType::Pointer sample = ListSampleType::New();

      for (typename FixedPointSetType::PointsContainerConstIterator fixed = fixedContainer->Begin(); fixed != fixedContainer->End(); ++fixed) {
        typename FixedPointSetType::PointType fixedPoint = fixed.Value();

        NeighborhoodType neighbors;
        std::vector<double> distancies;
        m_Tree->Search(fixedPoint, m_Radius, neighbors, distancies);
        if (distancies.size() == 0) {
          continue;
        }

        MeasureType distance = *std::min_element(distancies.begin(), distancies.end());
        sample->PushBack(distance);
        mean += distance;
        rmse += distance * distance;
        maximal = std::max(maximal, distance);
      }

      mean = mean / fixedPointSet->GetNumberOfPoints();
      rmse = std::sqrt(rmse / fixedPointSet->GetNumberOfPoints());

      // compute quantile
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
      MeasureType quantile = sampleToHistogram->GetOutput()->Quantile(0, m_LevelOfQuantile);

      metrics.resize(0);
      metrics.push_back(mean);
      metrics.push_back(rmse);
      metrics.push_back(quantile);
      metrics.push_back(maximal);
    }

    void InitializeKdTree(typename MovingPointSetType::ConstPointer pointSet)
    {
      m_ListOfPoints = ListSampleType::New();
      for (typename MovingPointSetType::PointsContainerConstIterator it = pointSet->GetPoints()->Begin(); it != pointSet->GetPoints()->End(); ++it) {
        m_ListOfPoints->PushBack(it.Value());
      }

      m_TreeGenerator = TreeGeneratorType::New();
      m_TreeGenerator->SetSample(m_ListOfPoints);
      m_TreeGenerator->SetBucketSize(m_BucketSize);
      try {
        m_TreeGenerator->Update();
      }
      catch (itk::ExceptionObject& excep) {
        itkExceptionMacro(<< excep);
      }

      m_Tree = m_TreeGenerator->GetOutput();
    }
  };
}
