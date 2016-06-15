#ifndef __PointSetToImageMetrics_h
#define __PointSetToImageMetrics_h

#include <itkPointSetToImageMetric.h>
#include <itkCovariantVector.h>
#include <itkPoint.h>

template< typename TFixedPointSet, typename TMovingImage >
class PointSetToImageMetrics :
public itk::PointSetToImageMetric < TFixedPointSet, TMovingImage >
{
  public:

  /** Standard class typedefs. */
  typedef PointSetToImageMetrics                      Self;
  typedef itk::PointSetToImageMetric< TFixedPointSet, TMovingImage > Superclass;
  typedef itk::SmartPointer< Self >                                  Pointer;
  typedef itk::SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetToImageMetrics, Object);

  /** Types transferred from the base class */
  typedef typename Superclass::RealType                RealType;
  typedef typename Superclass::TransformType           TransformType;
  typedef typename Superclass::TransformPointer        TransformPointer;
  typedef typename Superclass::TransformParametersType TransformParametersType;
  typedef typename Superclass::TransformJacobianType   TransformJacobianType;
  typedef typename Superclass::GradientPixelType       GradientPixelType;
  typedef typename Superclass::InputPointType          InputPointType;
  typedef typename Superclass::OutputPointType         OutputPointType;

  typedef typename Superclass::MeasureType               MeasureType;
  typedef typename Superclass::DerivativeType            DerivativeType;
  typedef typename Superclass::FixedPointSetType         FixedPointSetType;
  typedef typename Superclass::MovingImageType           MovingImageType;
  typedef typename Superclass::FixedPointSetConstPointer FixedPointSetConstPointer;
  typedef typename Superclass::MovingImageConstPointer   MovingImageConstPointer;

  typedef typename Superclass::PointIterator     PointIterator;
  typedef typename Superclass::PointDataIterator PointDataIterator;
  typedef itk::LinearInterpolateImageFunction<TMovingImage, double> InterpolatorType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType ParametersType;

  /**  Compute values. */
  void Compute();

  /*Get/Set values to compute quantile. */
  itkSetMacro(LevelOfQuantile, double);
  itkGetMacro(LevelOfQuantile, double);

  itkSetMacro(HistogramSize, unsigned int);
  itkGetMacro(HistogramSize, unsigned int);

  /*Get metrics values. */
  itkGetMacro(MeanValue, MeasureType);
  itkGetMacro(RMSEValue, MeasureType);
  itkGetMacro(QuantileValue, MeasureType);
  itkGetMacro(MaximalValue, MeasureType);

  void PrintReport(std::ostream& os) const;

  virtual MeasureType GetValue(const ParametersType& p) const
  {
    itkExceptionMacro(<< "Not implemented");
  }

  virtual void GetDerivative(const ParametersType& p, DerivativeType& dp) const
  {
    itkExceptionMacro(<< "Not implemented");
  }

  virtual void GetValueAndDerivative(const ParametersType& p, MeasureType& v, DerivativeType& dp) const
  {
    itkExceptionMacro(<< "Not implemented");
  }

protected:
  PointSetToImageMetrics();
  virtual ~PointSetToImageMetrics() {}

private:
  PointSetToImageMetrics(const Self &); //purposely not implemented
  void operator=(const Self &);                   //purposely not implemented

  typename InterpolatorType::Pointer m_Interpolator;

  double m_LevelOfQuantile;
  unsigned int m_HistogramSize;
  RealType m_FixedValue;
  MeasureType m_MeanValue = itk::NumericTraits<MeasureType>::Zero;
  MeasureType m_RMSEValue = itk::NumericTraits<MeasureType>::Zero;
  MeasureType m_QuantileValue = itk::NumericTraits<MeasureType>::Zero;
  MeasureType m_MaximalValue = itk::NumericTraits<MeasureType>::Zero;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "PointSetToImageMetrics.hxx"
#endif

#endif
