#pragma once

#include <itkImage.h>
#include <itkMeshToMeshFilter.h>
#include <itkLBFGSOptimizer.h>
#include <itkCompositeTransform.h>
#include <itkPointSet.h>

#include "ssmTransformInitializer.h"

namespace ssm
{
  template <typename TInputMesh, typename TOutputMesh = TInputMesh>
  class SurfaceToImageRegistrationMethod : public itk::MeshToMeshFilter < TInputMesh, TOutputMesh >
  {
  public:
    // Standard typedefs
    typedef SurfaceToImageRegistrationMethod Self;
    typedef itk::MeshToMeshFilter<TInputMesh, TOutputMesh> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::LBFGSOptimizer OptimizerType;
    typedef typename itk::Optimizer::ScalesType ScalesType;
    typedef itk::Transform<double, TInputMesh::PointDimension> TransformType;
    typedef itk::CompositeTransform<double, TInputMesh::PointDimension> CompositeTransformType;
    typedef itk::PointSet<float, TInputMesh::PointDimension> PointSetType;
    typedef itk::Image<unsigned char, TInputMesh::PointDimension> BinaryImageType;
    typedef itk::Image<float, TInputMesh::PointDimension> LevelsetImageType;
    typedef itk::PointSetToImageMetric<PointSetType, LevelsetImageType> MetricType;
    typedef itk::LinearInterpolateImageFunction<LevelsetImageType, double> InterpolatorType;

    itkNewMacro(Self);
    itkTypeMacro(SurfaceToImageRegistrationMethod, itk::MeshToMeshFilter);

    // set type of transform
    itkSetEnumMacro(TypeOfTransform, EnumTransform);
    itkGetEnumMacro(TypeOfTransform, EnumTransform);
    void SetTypeOfTransform(size_t transform) { m_TypeOfTransform = static_cast<EnumTransform>(transform); }

    //Set/Get PotentialImage
    itkGetConstObjectMacro(LevelsetImage, LevelsetImageType);
    itkSetConstObjectMacro(LevelsetImage, LevelsetImageType);

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

    void PrintReport(std::ostream& os) const;

  protected:
    SurfaceToImageRegistrationMethod();
    ~SurfaceToImageRegistrationMethod() {}

    virtual void GenerateData() override;
    void InitializeTransform();
    void ComputeLabelImage();
    void GenerateOutputData();

    itkStaticConstMacro(PointDimension, unsigned int, TInputMesh::PointDimension);
    static_assert(PointDimension == 3U, "Invalid dimension of input mesh. Dimension 3 is supported.");

    OptimizerType::Pointer m_Optimizer;
    typename TInputMesh::Pointer m_Surface;
    typename PointSetType::Pointer m_PointSet;
    typename BinaryImageType::ConstPointer m_Mask;
    typename LevelsetImageType::ConstPointer m_LevelsetImage;
    typename MetricType::Pointer m_Metric;
    typename TransformType::Pointer m_Transform;
    typename CompositeTransformType::Pointer m_CompositeTransform;
    typename InterpolatorType::Pointer m_Interpolator;

    EnumTransform m_TypeOfTransform = EnumTransform::Euler3D;
    size_t m_NumberOfIterations = 500;
    double m_LineSearchAccuracy = 0.1;
    double m_DefaultStepLength = 0.1;
    double m_GradientConvergenceTolerance = 1e-07;
    ScalesType m_Scales;
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmSurfaceToImageRegistrationMethod.hxx"
#endif
