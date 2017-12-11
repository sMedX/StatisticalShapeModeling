#pragma once

#include <itkImage.h>
#include <itkTimeProbe.h>
#include <itkStatisticalModel.h>
#include <itkTransform.h>
#include <itkLBFGSOptimizer.h>
#include <itkLogger.h>
#include <itkSingleValuedCostFunction.h>

#include "ssmShapeModelMultiTransform.h"
#include "ssmShapeModelToLevelSetImageMetric.h"

namespace ssm
{
  template <typename TShapeModel, typename TOutputMesh>
  class ShapeModelRegistrationMethodBase : public itk::ProcessObject
  {
  public:
    /** Standard class typedefs. */
    typedef ShapeModelRegistrationMethodBase       Self;
    typedef itk::ProcessObject                     Superclass;
    typedef itk::SmartPointer< Self >              Pointer;
    typedef itk::SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkTypeMacro(ShapeModelRegistrationMethodBase, Superclass);

    /** Constants for the image dimensions */
    itkStaticConstMacro(PointDimension, unsigned int, 3U);
    static_assert(TShapeModel::RepresenterType::DatasetType::PointDimension == PointDimension, "Invalid dimension of input shape model. Dimension 3 is supported.");

    /** Typedefs of the shape model. */
    typedef TShapeModel ShapeModelType;
    typedef TOutputMesh OutputMeshType;
    typedef typename ShapeModelType::RepresenterType::DatasetType DatasetType;
    typedef typename itk::PointSet<typename DatasetType::PixelType, PointDimension> PointSetType;

    /**  Type of the Transform Base class */
    typedef double CoordinateRepresentationType;
    typedef itk::LBFGSOptimizer OptimizerType;
    typedef itk::SingleValuedCostFunction MetricType;
    typedef typename itk::Optimizer::ScalesType ScalesType;
    typedef itk::Transform<double, PointDimension, PointDimension> SpatialTransformType;
    typedef ShapeModelMultiTransform<DatasetType, double> ShapeModelMultiTransformType;

    /** Returns the transform resulting from the registration process  */
    const OutputMeshType* GetOutput() const;

    /** Method to return the latest modified time of this object or any of its cached vars */
    virtual itk::ModifiedTimeType GetMTime() const ITK_OVERRIDE;

    // Set logger
    itkSetObjectMacro(Logger, itk::Logger);

    // Set/Get shape model and parameters
    itkSetConstObjectMacro(ShapeModel, ShapeModelType);
    itkGetConstObjectMacro(ShapeModel, ShapeModelType);

    // Set/Get transform and optimizer
    itkSetObjectMacro(SpatialTransform, SpatialTransformType);
    itkGetObjectMacro(SpatialTransform, SpatialTransformType);
    itkGetObjectMacro(Transform, ShapeModelMultiTransformType);

    itkGetMacro(Scales, ScalesType);

    itkSetMacro(SpatialScales, ScalesType);
    itkGetMacro(SpatialScales, ScalesType);

    itkSetMacro(ModelScale, double);
    itkGetMacro(ModelScale, double);

    itkSetMacro(NumberOfIterations, size_t);
    itkGetMacro(NumberOfIterations, size_t);

    itkSetObjectMacro(Optimizer, OptimizerType);
    itkGetObjectMacro(Optimizer, OptimizerType);

    itkSetObjectMacro(Metric, MetricType);
    itkGetObjectMacro(Metric, MetricType);

    itkGetMacro(ElapsedTime, double);

    void SetNumberOfEpochs(const unsigned int & input) { m_NumberOfEpochs = input; m_UseNumberOfEpochs = true; };
    itkGetMacro(NumberOfEpochs, unsigned int);

    void SetStep(const unsigned int & input) { m_Step = input; m_UseNumberOfEpochs = false; };
    itkGetMacro(Step, unsigned int);

    void PrintReport();

  protected:
    ShapeModelRegistrationMethodBase();
    ~ShapeModelRegistrationMethodBase() {}

    void Initialize();
    void InitializeTransform();
    void InitializeOptimizer();
    void ComputeMultiStageData();
    void PerformRegistration();
    void GenerateOutputData();
    void ExceptionHandler(itk::ExceptionObject & excep, const size_t & line);

    typename ShapeModelType::ConstPointer m_ShapeModel;

    typename OptimizerType::Pointer m_Optimizer;
    typename OutputMeshType::Pointer m_OutputMesh;
    OptimizerType::ParametersType m_InitialParameters;

    typename MetricType::Pointer m_Metric;

    typename ShapeModelMultiTransformType::Pointer m_Transform;
    typename SpatialTransformType::Pointer m_SpatialTransform;
    ScalesType m_SpatialScales;
    ScalesType m_Scales;
    ScalesType m_InverseScales;
    double m_ModelScale;

    bool m_UseNumberOfEpochs;
    unsigned int m_NumberOfEpochs;
    unsigned int m_Step;
    itk::Array<unsigned int> m_NumberOfUsedComponents;
    itk::Array<double> m_CostFunctionValues;

    size_t m_NumberOfIterations;

    itk::TimeProbe m_Time;
    double m_ElapsedTime;
    itk::Logger::Pointer m_Logger;
    std::ostringstream m_Message;
  };
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmShapeModelRegistrationMethodBase.hxx"
#endif
