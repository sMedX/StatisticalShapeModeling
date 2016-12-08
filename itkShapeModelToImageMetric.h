#pragma once

#include <itkMacro.h>
#include <itkImageBase.h>
#include <itkTransform.h>
#include <itkInterpolateImageFunction.h>
#include <itkGradientRecursiveGaussianImageFilter.h>

namespace itk
{
template <typename TShapeModel,  typename TImage>
class ShapeModelToImageMetric:public itk::SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef ShapeModelToImageMetric          Self;
  typedef itk::SingleValuedCostFunction    Superclass;
  typedef itk::SmartPointer<Self>          Pointer;
  typedef itk::SmartPointer<const Self>    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShapeModelToImageMetric, itk::Object);

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
  typedef itk::Transform< CoordinateRepresentationType, itkGetStaticConstMacro(ImageDimension), itkGetStaticConstMacro(PointDimension) > TransformType;
  typedef typename TransformType::Pointer         TransformPointer;
  typedef typename TransformType::InputPointType  InputPointType;
  typedef typename TransformType::OutputPointType OutputPointType;
  typedef typename TransformType::ParametersType  TransformParametersType;
  typedef typename TransformType::JacobianType    TransformJacobianType;

  /**  Type of the Interpolator Base class */
  typedef itk::InterpolateImageFunction<ImageType, CoordinateRepresentationType> InterpolatorType;

  /** Gaussian filter to compute the gradient of the Moving Image */
  typedef typename itk::NumericTraits< ImagePixelType >::RealType RealType;
  typedef itk::CovariantVector<RealType, itkGetStaticConstMacro(ImageDimension)> GradientPixelType;
  typedef itk::Image<GradientPixelType, itkGetStaticConstMacro(ImageDimension)> GradientImageType;
  typedef itk::SmartPointer<GradientImageType> GradientImagePointer;
  typedef typename itk::GradientRecursiveGaussianImageFilter<TImage, GradientImageType> GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  /**  Type of the measure. */
  typedef Superclass::MeasureType MeasureType;

  /**  Type of the derivative. */
  typedef Superclass::DerivativeType DerivativeType;

  /**  Type of the parameters. */
  typedef Superclass::ParametersType ParametersType;

  /** Get/Set the shape model.  */
  itkSetConstObjectMacro(ShapeModel, TShapeModel);
  itkGetConstObjectMacro(ShapeModel, TShapeModel);

  /** Get/Set the Moving Image.  */
  itkSetConstObjectMacro(Image, ImageType);
  itkGetConstObjectMacro(Image, ImageType);

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
  itkGetConstReferenceMacro(NumberOfPixelsCounted, itk::SizeValueType);

  /** Set the parameters defining the Transform. */
  void SetTransformParameters(const ParametersType & parameters) const;

  itkSetMacro(RegularizationParameter, MeasureType)
  itkGetMacro(RegularizationParameter, MeasureType)

  /** Set/Get the flag for computing the image gradient.
   *  When ON the metric derivative is computed using the Jacobian of the
   *  transformation and the image gradient. When OFF the metric derivative
   *  is computed by finite differences. Mode ON results in higher speed
   *  at the price of large memory footprint. Mode OFF results in small
   *  memory footprint at the price of large computation time */
  itkSetMacro(ComputeGradient, bool);
  itkGetConstReferenceMacro(ComputeGradient, bool);

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE { return m_Transform->GetNumberOfParameters(); }

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
  void CalculateValuePenalty(const TransformParametersType& parameters, MeasureType & value) const;
  void CalculateDerivativePenalty(const TransformParametersType& parameters, DerivativeType & derivative) const;

protected:
  ShapeModelToImageMetric();
  virtual ~ShapeModelToImageMetric() {}
  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

  mutable itk::SizeValueType m_NumberOfPixelsCounted = 0;
  typename PointsContainerType::Pointer m_PointsContainer;
  ImageConstPointer m_Image;
  mutable TransformPointer m_Transform;
  InterpolatorPointer m_Interpolator;
  bool m_ComputeGradient;
  GradientImagePointer m_GradientImage;
  typename TShapeModel::ConstPointer m_ShapeModel;
  MeasureType m_RegularizationParameter = 0.1;
  unsigned int m_NumberOfParameters;
  unsigned int m_NumberOfComponents;

private:
  ShapeModelToImageMetric(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShapeModelToImageMetric.hxx"
#endif
