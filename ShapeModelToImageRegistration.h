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
#include <utils/itkPenalizingMeanSquaresPointSetToImageMetric.h>

template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
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
  typedef itk::PenalizingMeanSquaresPointSetToImageMetric<PointSetType, PotentialImageType> MetricType;

  itkNewMacro(Self);
  itkTypeMacro(ShapeModelToImageRegistration, itk::ImageToMeshFilter);

  //Set/Get shape model
  itkSetConstObjectMacro(ShapeModel, ModelType);
  itkGetConstObjectMacro(ShapeModel, ModelType);

  //Get PotentialImage
  itkGetConstObjectMacro(Optimizer, OptimizerType);
  itkGetConstObjectMacro(PotentialImage, PotentialImageType);

  itkSetConstObjectMacro(InputTransform, TInputTransformType);
  itkGetConstObjectMacro(InputTransform, TInputTransformType);

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

  void PrintReport(std::ostream& os) const;
  TOutputMesh* GetOutput() { return  dynamic_cast<TOutputMesh*>(this->ProcessObject::GetOutput(0)); }
  TOutputMesh* GetMovedOutput() { return  dynamic_cast<TOutputMesh*>(this->ProcessObject::GetOutput(1)); }

protected:
  ShapeModelToImageRegistration();
  ~ShapeModelToImageRegistration() {}

  virtual void GenerateData() override;
  void InitializeTransform();
  void ComputePotentialImage();
  void GenerateOutputData();

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);
  static_assert(ImageDimension == 3U, "Invalid dimension of input image. Dimension 3 is supported.");

  OptimizerType::Pointer m_Optimizer;
  typename ModelType::ConstPointer m_ShapeModel;
  typename PotentialImageType::ConstPointer m_PotentialImage;
  typename TOutputMesh::Pointer m_Surface;
  typename PointSetType::Pointer m_PointSet;
  typename MetricType::Pointer m_Metric;

  typename TransformType::Pointer m_Transform;
  typename TInputTransformType::ConstPointer m_InputTransform;
  typename ModelTransformType::Pointer m_ModelTransform;

  int m_ModelTransformIndex;

  unsigned int m_NumberOfIterations = 100;
  double m_GradientConvergenceTolerance = 1e-07;
  double m_LineSearchAccuracy = 0.1;
  double m_DefaultStepLength = 0.1;
  double m_RegularizationParameter = 0.1;

  ScalesType m_Scales;
  double m_ModelScale = 3;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "ShapeModelToImageRegistration.hxx"
#endif

#endif
