#ifndef __ShapeModelToImageRegistration_h
#define __ShapeModelToImageRegistration_h

#include <itkImage.h>
#include <itkImageToMeshFilter.h>
#include <itkLBFGSOptimizer.h>
#include <itkCompositeTransform.h>
#include <itkPointSet.h>
#include <itkStatisticalModel.h>
#include <itkStatisticalShapeModelTransform.h>
#include <itkPointSetToImageRegistrationMethod.h>

template <typename TInputImage, typename TOutputMesh>
class ShapeModelToImageRegistration : public itk::ImageToMeshFilter<TInputImage, TOutputMesh>
{
public:
  // Standard typedefs
  typedef ShapeModelToImageRegistration Self;
  typedef itk::ImageToMeshFilter<TInputImage, TOutputMesh> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::LBFGSOptimizer OptimizerType;
  typedef typename itk::Optimizer::ScalesType ScalesType;
  typedef itk::CompositeTransform<double, TInputImage::ImageDimension> TransformType;
  typedef itk::StatisticalShapeModelTransform<TOutputMesh, double, TInputImage::ImageDimension> ModelTransformType;
  typedef itk::PointSet<float, TInputImage::ImageDimension> PointSetType;
  typedef itk::Image<unsigned char, TInputImage::ImageDimension> BinaryImageType;
  typedef itk::Image<float, TInputImage::ImageDimension> PotentialImageType;
  typedef itk::StatisticalModel<TOutputMesh> ModelType;
  typedef itk::PointSetToImageMetric<PointSetType, PotentialImageType> MetricType;

  itkNewMacro(Self);
  itkTypeMacro(ShapeModelToImageRegistration, itk::ImageToMeshFilter);

  //Set/Get shape model
  itkSetConstObjectMacro(ShapeModel, ModelType);
  itkGetConstObjectMacro(ShapeModel, ModelType);

  //Get PotentialImage
  itkGetConstObjectMacro(Optimizer, OptimizerType);
  itkGetConstObjectMacro(PotentialImage, PotentialImageType);
  itkGetConstObjectMacro(Transform, TransformType);

  itkSetMacro(NumberOfIterations, unsigned int);
  itkGetMacro(NumberOfIterations, unsigned int);

  itkSetMacro(GradientConvergenceTolerance, double);
  itkGetMacro(GradientConvergenceTolerance, double);

  itkSetMacro(LineSearchAccuracy, double);
  itkGetMacro(LineSearchAccuracy, double);

  itkSetMacro(DefaultStepLength, double);
  itkGetMacro(DefaultStepLength, double);

  itkSetMacro(RegularizationParameter, double);
  itkGetMacro(RegularizationParameter, double);

  //Set/Get scaling
  itkSetMacro(ModelScale, double);
  itkGetMacro(ModelScale, double);

  itkSetMacro(TranslationScale, double);
  itkGetMacro(TranslationScale, double);

  itkSetMacro(ScalingScale, double);
  itkGetMacro(ScalingScale, double);

  itkSetMacro(RotationScale, double);
  itkGetMacro(RotationScale, double);

  void PrintReport(std::ostream& os) const;
  TOutputMesh* GetOutput() { return  dynamic_cast<TOutputMesh*>(this->ProcessObject::GetOutput(0)); }
  TOutputMesh* GetDeformedOutput() { return  dynamic_cast<TOutputMesh*>(this->ProcessObject::GetOutput(1)); }

protected:
  ShapeModelToImageRegistration();
  ~ShapeModelToImageRegistration() {}

  virtual void GenerateData() override;
  void InitializeTransform();
  void ComputePotentialImage();
  void PointSetToImageRegistration();
  void GenerateOutputData();

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);
  static_assert(ImageDimension == 3U, "Invalid dimension of input image. Dimension 3 is supported.");

  OptimizerType::Pointer m_Optimizer;
  typename ModelType::ConstPointer m_ShapeModel;
  typename TInputImage::ConstPointer m_Image;
  typename PotentialImageType::ConstPointer m_PotentialImage;
  typename TOutputMesh::Pointer m_Surface;
  typename PointSetType::Pointer m_PointSet;
  typename TransformType::Pointer m_Transform;
  typename MetricType::Pointer m_Metric;
  typename ModelTransformType::Pointer m_ModelTransform;

  int m_ModelTransformIndex;

  unsigned int m_NumberOfIterations = 100;
  double m_GradientConvergenceTolerance = 1e-07;
  double m_LineSearchAccuracy = 0.1;
  double m_DefaultStepLength = 0.1;
  double m_RegularizationParameter = 0.1;

  ScalesType m_Scales;
  double m_ModelScale = 3;
  double m_TranslationScale = 1;
  double m_ScalingScale = 1;
  double m_RotationScale = 0.1;
  double m_Margin = 5;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "ShapeModelToImageRegistration.hxx"
#endif

#endif
