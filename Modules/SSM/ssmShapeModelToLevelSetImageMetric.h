#pragma once

#include <itkMacro.h>
#include <itkImageBase.h>
#include <itkTransform.h>
#include <itkInterpolateImageFunction.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkCovariantVector.h>
#include <itkLogger.h>

namespace ssm
{
template <typename TShapeModel,  typename TTransform, typename TImage>
class ShapeModelToLevelSetImageMetric:public itk::SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef ShapeModelToLevelSetImageMetric  Self;
  typedef itk::SingleValuedCostFunction    Superclass;
  typedef itk::SmartPointer<Self>          Pointer;
  typedef itk::SmartPointer<const Self>    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShapeModelToLevelSetImageMetric, itk::Object);

  /** Type used for representing point components  */
  typedef Superclass::ParametersValueType CoordinateRepresentationType;

  /**  Type of the moving Image. */
  typedef TImage                           ImageType;
  typedef typename TImage::PixelType       ImagePixelType;
  typedef typename ImageType::ConstPointer ImageConstPointer;

  /**  Type of the fixed Image. */
  typedef typename TShapeModel::RepresenterType RepresenterType;
  typedef typename RepresenterType::DatasetType DatasetType;
  typedef typename itk::PointSet<typename DatasetType::PixelType, DatasetType::PointDimension> PointSetType;
  typedef typename PointSetType::PointsContainer PointsContainerType;
  typedef typename PointSetType::ConstPointer PointSetConstPointer;
  typedef typename PointSetType::PointsContainer::ConstIterator PointIteratorType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);
  itkStaticConstMacro(PointDimension, unsigned int, TShapeModel::RepresenterType::DatasetType::PointDimension);

  /**  Type of the Transform Base class */
  typedef TTransform                              TransformType;
  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the Interpolator Base class */
  typedef itk::InterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;

  /** Gaussian filter to compute the gradient of the Moving Image */
  typedef typename itk::NumericTraits< ImagePixelType >::RealType RealType;
  typedef itk::CovariantVector<RealType, ImageDimension> GradientPixelType;
  typedef itk::Image<GradientPixelType, ImageDimension> GradientImageType;
  typedef typename itk::GradientRecursiveGaussianImageFilter<TImage, GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  // Set logger
  itkSetObjectMacro(Logger, itk::Logger);

  /** Get/Set the shape model.  */
  itkSetConstObjectMacro(ShapeModel, TShapeModel);
  itkGetConstObjectMacro(ShapeModel, TShapeModel);

  /** Get/Set the Moving Image.  */
  itkSetConstObjectMacro(LevelSetImage, ImageType);
  itkGetConstObjectMacro(LevelSetImage, ImageType);

  /** Connect the Transform. */
  itkSetObjectMacro(Transform, TransformType);

  /** Get a pointer to the Transform.  */
  itkGetModifiableObjectMacro(Transform, TransformType);

  /** Connect the Interpolator. */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the Interpolator.  */
  itkGetModifiableObjectMacro(Interpolator, InterpolatorType);

  /** Get Gradient Image. */
  itkGetModifiableObjectMacro(GradientImage, GradientImageType);

  /** Get the number of pixels considered in the computation. */
  itkGetConstReferenceMacro(NumberOfSamplesCounted, itk::SizeValueType);

  itkSetMacro(RegularizationParameter, double);
  itkGetMacro(RegularizationParameter, double);

  itkSetMacro(Degree, size_t);
  itkGetMacro(Degree, size_t);

  itkSetMacro(NumberOfThreads, size_t);
  itkGetMacro(NumberOfThreads, size_t);
  itkGetMacro(MaximalNumberOfThreads, size_t);

  itkSetMacro(NumberOfUsedPoints, size_t);
  itkGetMacro(NumberOfUsedPoints, size_t);

  itkGetMacro(IsInitialized, bool);

  /** Set/Get the flag for computing the image gradient */
  itkSetMacro(ComputeGradient, bool);
  itkGetConstReferenceMacro(ComputeGradient, bool);

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE {return m_NumberOfParameters;}

  /** Initialize the Metric by making sure that all the components are present and plugged together correctly     */
  virtual void Initialize(void)
  throw(itk::ExceptionObject);

  /** Get the derivatives of the match measure. */
  void GetDerivative(const TransformParametersType & parameters, DerivativeType & Derivative) const ITK_OVERRIDE;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative(const TransformParametersType & parameters, MeasureType & Value, DerivativeType & Derivative) const ITK_OVERRIDE;

  /**Compute penalty.  */
  void CalculateValueAndDerivativePenalty(const TransformParametersType & parameters, MeasureType & value, DerivativeType  & derivative) const;

  void PrintReport() const;

protected:
  ShapeModelToLevelSetImageMetric();
  virtual ~ShapeModelToLevelSetImageMetric() {}

  mutable itk::SizeValueType m_NumberOfSamplesCounted;
  typename PointsContainerType::Pointer m_PointsContainer;
  ImageConstPointer m_LevelSetImage;
  mutable TransformPointer m_Transform;
  InterpolatorPointer m_Interpolator;
  bool m_ComputeGradient;
  typename GradientImageType::Pointer m_GradientImage;
  typename TShapeModel::ConstPointer m_ShapeModel;
  double m_RegularizationParameter;
  unsigned int m_Degree;
  mutable size_t m_NumberOfEvaluations;
  bool m_IsInitialized;
  size_t m_NumberOfUsedPoints;
  size_t m_NumberOfParameters;

  // multi threading members and methods
  struct PerThreadData
  {
    itk::SizeValueType    m_NumberOfSamplesCounted;
    MeasureType           m_Value;
    DerivativeType        m_Derivative;
    TransformJacobianType m_ModelJacobian;
    TransformJacobianType m_SpatialJacobian;
    TransformJacobianType m_Jacobian;
    TransformJacobianType m_JacobianCache;
    PointIteratorType     m_Begin;
    PointIteratorType     m_End;
  };
  size_t m_MaximalNumberOfThreads;
  size_t m_NumberOfThreads;
  mutable std::vector<PerThreadData> m_Threads;

  inline void GetValueAndDerivativeThreadProcessSample(PerThreadData & data) const;
  void MultiThreadingInitialize();
  void SparsePoints();

  itk::Logger::Pointer m_Logger;
  mutable std::ostringstream m_Message;

private:
  ShapeModelToLevelSetImageMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmShapeModelToLevelSetImageMetric.hxx"
#endif
