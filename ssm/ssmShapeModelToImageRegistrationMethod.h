#pragma once

#include <itkTimeProbe.h>
#include <itkImage.h>
#include <itkStatisticalModel.h>
#include <itkTransform.h>
#include <itkCompositeTransform.h>
#include <itkStatisticalShapeModelTransform.h>
#include <itkLBFGSOptimizer.h>

#include "ssmShapeModelToImageMetric.h"
#include "ssmTransformInitializer.h"

namespace ssm
{
  template <typename TShapeModel, typename TOutputMesh>
  class ShapeModelToImageRegistrationMethod : public itk::ProcessObject
  {
  public:
    /** Standard class typedefs. */
    typedef ShapeModelToImageRegistrationMethod       Self;
    typedef itk::ProcessObject                        Superclass;
    typedef itk::SmartPointer< Self >                 Pointer;
    typedef itk::SmartPointer< const Self >           ConstPointer;

    /** Constants for the image dimensions */
    itkStaticConstMacro(ImageDimension, unsigned int, 3);
    itkStaticConstMacro(PointDimension, unsigned int, TShapeModel::RepresenterType::DatasetType::PointDimension);
    static_assert(PointDimension == 3U, "Invalid dimension of input shape model. Dimension 3 is supported.");
    static_assert(ImageDimension == 3U, "Invalid dimension of input shape model. Dimension 3 is supported.");

    /** Typedefs of the shape model. */
    typedef TShapeModel ShapeModelType;
    typedef TOutputMesh OutputMeshType;
    typedef typename ShapeModelType::RepresenterType RepresenterType;
    typedef typename RepresenterType::DatasetType DatasetType;
    typedef typename itk::PointSet<typename DatasetType::PixelType, DatasetType::PointDimension> PointSetType;

    /** Typedefs of the input image. */
    typedef itk::Image<float, itkGetStaticConstMacro(ImageDimension)> LevelSetImageType;

    /**  Type of the Transform Base class */
    typedef double CoordinateRepresentationType;
    typedef itk::LBFGSOptimizer OptimizerType;
    typedef typename itk::Optimizer::ScalesType ScalesType;
    typedef typename itk::StatisticalShapeModelTransform<DatasetType, double, PointDimension> ShapeTransformType;
    typedef itk::ShapeModelToImageMetric<TShapeModel, LevelSetImageType> MetricType;
    typedef ssm::TransformInitializer<CoordinateRepresentationType> TransformInitializerType;
    typedef itk::Transform<CoordinateRepresentationType, PointDimension> TransformType;
    typedef itk::CompositeTransform<CoordinateRepresentationType, PointDimension> CompositeTransformType;

    /** Smart Pointer type to a DataObject. */
    typedef typename itk::DataObject::Pointer DataObjectPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(ShapeModelToImageRegistrationMethod, itk::ProcessObject);

    /** Returns the transform resulting from the registration process  */
    const OutputMeshType* GetOutput() const;

    /** Method to return the latest modified time of this object or any of its cached ivars */
    virtual itk::ModifiedTimeType GetMTime() const ITK_OVERRIDE;

    // Set/Get shape model and parameters
    itkSetConstObjectMacro(ShapeModel, ShapeModelType);
    itkGetConstObjectMacro(ShapeModel, ShapeModelType);

    // Set/Get level set image
    itkGetConstObjectMacro(LevelSetImage, LevelSetImageType);
    itkSetConstObjectMacro(LevelSetImage, LevelSetImageType);

    // Set/Get transform and optimizer
    itkSetObjectMacro(ShapeTransform, ShapeTransformType);
    itkGetObjectMacro(ShapeTransform, ShapeTransformType);

    itkSetMacro(Scales, ScalesType);
    itkGetMacro(Scales, ScalesType);

    itkGetConstObjectMacro(Optimizer, OptimizerType);

    // Set/Get optimizer parameters
    itkSetMacro(NumberOfIterations, unsigned int);
    itkGetMacro(NumberOfIterations, unsigned int);

    itkSetMacro(RegularizationParameter, double);
    itkGetMacro(RegularizationParameter, double);

    // Set/Get model scales
    itkSetMacro(ModelScale, double);
    itkGetMacro(ModelScale, double);

    itkSetMacro(Degree, unsigned int);
    itkGetMacro(Degree, unsigned int);

    itkSetObjectMacro(TransformInitializer, TransformInitializerType);
    itkGetObjectMacro(TransformInitializer, TransformInitializerType);

    void PrintReport(std::ostream& os);
    itk::TimeProbe::CountType GetElapsedTime() { return m_Clock.GetTotal(); }

  protected:
    ShapeModelToImageRegistrationMethod();
    ~ShapeModelToImageRegistrationMethod() {}

    void GenerateData();
    void InitializeTransform();
    void InitializeOptimizer();
    void InitializeMetric();
    void GenerateOutputData();

    typename ShapeModelType::ConstPointer m_ShapeModel;
    typename LevelSetImageType::ConstPointer m_LevelSetImage;
    typename OptimizerType::Pointer m_Optimizer;
    typename OutputMeshType::Pointer m_OutputMesh;
    typename MetricType::Pointer m_Metric;
    OptimizerType::ParametersType m_InitialParameters;
    size_t m_NumberOfSpatialParameters;
    size_t m_NumberOfShapeModelParameters;

    typename TransformType::Pointer m_Transform;
    typename TransformType::Pointer m_SpatialTransform;
    typename ShapeTransformType::Pointer m_ShapeTransform;
    TransformInitializerType::Pointer m_TransformInitializer;
    ScalesType m_SpatialScales;
    ScalesType m_Scales;
    double m_ModelScale = 3;
    unsigned int m_NumberOfIterations = 500;
    double m_LineSearchAccuracy = 0.1;
    double m_DefaultStepLength = 0.1;
    double m_GradientConvergenceTolerance = 1e-05;
    double m_RegularizationParameter = 0.1;
    unsigned int m_Degree = 2;
    itk::TimeProbe m_Clock;

  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmShapeModelToImageRegistrationMethod.hxx"
#endif
