#ifndef __SurfaceToImageRegistrationFilter_h
#define __SurfaceToImageRegistrationFilter_h

#include <itkImage.h>
#include <itkMeshToMeshFilter.h>
#include <itkLBFGSOptimizer.h>
#include <itkCompositeTransform.h>
#include <itkPointSet.h>
#include <itkPointSetToImageRegistrationMethod.h>

template <typename TInputMesh, typename TOutputMesh=TInputMesh>
class SurfaceToImageRegistrationFilter : public itk::MeshToMeshFilter<TInputMesh, TOutputMesh>
{
public:
  // Standard typedefs
  typedef SurfaceToImageRegistrationFilter Self;
  typedef itk::MeshToMeshFilter<TInputMesh, TOutputMesh> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;
  typedef itk::LBFGSOptimizer OptimizerType;
  typedef typename itk::Optimizer::ScalesType ScalesType;
  typedef itk::CompositeTransform<double, TInputMesh::PointDimension> TransformType;
  typedef itk::PointSet<float, TInputMesh::PointDimension> PointSetType;
  typedef itk::Image<unsigned char, TInputMesh::PointDimension> BinaryImageType;
  typedef itk::Image<float, TInputMesh::PointDimension> PotentialImageType;
  typedef itk::PointSetToImageMetric<PointSetType, PotentialImageType> MetricType;

  itkNewMacro(Self);
  itkTypeMacro(surfaceToImageRegistration, itk::MeshToMeshFilter);

  //Set/Get shape model
  itkGetConstObjectMacro(Mask, BinaryImageType);

  //Set/Get PotentialImage
  itkGetConstObjectMacro(PotentialImage, PotentialImageType);
  itkSetConstObjectMacro(PotentialImage, PotentialImageType);

  itkGetConstObjectMacro(Optimizer, OptimizerType);
  itkGetConstObjectMacro(Transform, TransformType);

  itkSetMacro(NumberOfIterations, unsigned int);
  itkGetMacro(NumberOfIterations, unsigned int);

  itkSetMacro(GradientConvergenceTolerance, double);
  itkGetMacro(GradientConvergenceTolerance, double);

  itkSetMacro(LineSearchAccuracy, double);
  itkGetMacro(LineSearchAccuracy, double);

  itkSetMacro(DefaultStepLength, double);
  itkGetMacro(DefaultStepLength, double);

  itkSetMacro(TranslationScale, double);
  itkGetMacro(TranslationScale, double);

  itkSetMacro(ScalingScale, double);
  itkGetMacro(ScalingScale, double);

  itkSetMacro(RotationScale, double);
  itkGetMacro(RotationScale, double);

  void PrintReport(std::ostream& os) const;

protected:
  SurfaceToImageRegistrationFilter();
  ~SurfaceToImageRegistrationFilter() {}

  virtual void GenerateData() override;
  void InitializeTransform();
  void ComputeLabelImage();
  void GenerateOutputData();

  itkStaticConstMacro(PointDimension, unsigned int, TInputMesh::PointDimension);
  static_assert(PointDimension == 3U, "Invalid dimension of input image. Dimension 3 is supported.");

  OptimizerType::Pointer m_Optimizer;
  typename TInputMesh::Pointer m_Surface;
  typename PointSetType::Pointer m_PointSet;
  typename BinaryImageType::ConstPointer m_Mask;
  typename PotentialImageType::ConstPointer m_PotentialImage;
  typename TransformType::Pointer m_Transform;
  typename MetricType::Pointer m_Metric;

  unsigned int m_NumberOfIterations = 100;
  double m_GradientConvergenceTolerance = 1e-07;
  double m_LineSearchAccuracy = 0.1;
  double m_DefaultStepLength = 0.1;

  ScalesType m_Scales;
  double m_TranslationScale = 1;
  double m_ScalingScale = 1;
  double m_RotationScale = 0.1;
  double m_Margin = 5;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "SurfaceToImageRegistrationFilter.hxx"
#endif

#endif
