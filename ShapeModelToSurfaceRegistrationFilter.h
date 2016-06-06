#ifndef __ShapeModelToSurfaceRegistrationFilter_h
#define __ShapeModelToSurfaceRegistrationFilter_h

#include <itkImage.h>
#include <itkMeshToMeshFilter.h>
#include <itkLBFGSOptimizer.h>
#include <itkCompositeTransform.h>
#include <itkPointSet.h>
#include <itkStatisticalModel.h>
#include <itkStatisticalShapeModelTransform.h>
#include <itkPointSetToImageRegistrationMethod.h>
#include <utils/itkPenalizingMeanSquaresPointSetToImageMetric.h>

template <typename TInputMesh, typename TOutputMesh = TInputMesh>
class ShapeModelToSurfaceRegistrationFilter : public itk::MeshToMeshFilter<TInputMesh, TOutputMesh>
{
public:
  // Standard typedefs
  typedef ShapeModelToSurfaceRegistrationFilter Self;
  typedef itk::MeshToMeshFilter<TInputMesh, TOutputMesh> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::LBFGSOptimizer OptimizerType;
  typedef typename itk::Optimizer::ScalesType ScalesType;
  typedef itk::StatisticalShapeModelTransform<TOutputMesh, double, TInputMesh::PointDimension> ModelTransformType;
  typedef itk::PointSet<float, TInputMesh::PointDimension> PointSetType;
  typedef itk::Image<unsigned char, TInputMesh::PointDimension> BinaryImageType;
  typedef itk::Image<float, TInputMesh::PointDimension> LevelsetImageType;
  typedef itk::StatisticalModel<TOutputMesh> ModelType;
  typedef itk::PenalizingMeanSquaresPointSetToImageMetric<PointSetType, LevelsetImageType> MetricType;

  itkNewMacro(Self);
  itkTypeMacro(ShapeModelToSurfaceRegistrationFilter, itk::ImageToMeshFilter);

  //Set/Get shape model
  itkSetConstObjectMacro(ShapeModel, ModelType);
  itkGetConstObjectMacro(ShapeModel, ModelType);

  //Get PotentialImage
  itkGetConstObjectMacro(Optimizer, OptimizerType);
  itkGetConstObjectMacro(LevelsetImage, LevelsetImageType);

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

  itkSetMacro(Margin, double);
  itkGetMacro(Margin, double);

  itkSetMacro(Spacing, typename LevelsetImageType::SpacingType);
  itkGetMacro(Spacing, typename LevelsetImageType::SpacingType);

  //Set/Get scaling
  itkSetMacro(ModelScale, double);
  itkGetMacro(ModelScale, double);

  void PrintReport(std::ostream& os) const;

protected:
  ShapeModelToSurfaceRegistrationFilter();
  ~ShapeModelToSurfaceRegistrationFilter() {}

  virtual void GenerateData() override;
  void InitializeTransform();
  void GenerateOutputData();

  itkStaticConstMacro(Dimension, unsigned int, TInputMesh::PointDimension);
  static_assert(Dimension == 3U, "Invalid dimension of input image. Dimension 3 is supported.");

  OptimizerType::Pointer m_Optimizer;
  typename ModelType::ConstPointer m_ShapeModel;
  typename LevelsetImageType::ConstPointer m_LevelsetImage;
  typename TInputMesh::Pointer m_Surface;
  typename PointSetType::Pointer m_PointSet;
  typename MetricType::Pointer m_Metric;
  typename ModelTransformType::Pointer m_Transform;

  int m_ModelTransformIndex;

  unsigned int m_NumberOfIterations = 100;
  double m_GradientConvergenceTolerance = 1e-07;
  double m_LineSearchAccuracy = 0.1;
  double m_DefaultStepLength = 0.1;
  double m_RegularizationParameter = 0.1;

  ScalesType m_Scales;
  typename LevelsetImageType::SpacingType m_Spacing;
  double m_ModelScale = 3;
  double m_Margin = 0.25;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "ShapeModelToSurfaceRegistrationFilter.hxx"
#endif

#endif
