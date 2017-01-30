#pragma once

#include <itkTransformMeshFilter.h>
#include <itkLinearInterpolateImageFunction.h>

#include "ssmShapeModelToImageRegistrationMethod.h"

namespace ssm
{
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::ShapeModelToImageRegistrationMethod()
  {
    this->SetNumberOfRequiredInputs(0);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, TOutputMesh::New());

    m_ShapeModel = nullptr;
    m_LevelSetImage = nullptr;
    m_Metric = nullptr;
    m_TransformInitializer = nullptr;
    m_SpatialTransform = nullptr;
    m_ShapeTransform = nullptr;
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  const typename ShapeModelToImageRegistrationMethod< TShapeModel, TOutputMesh>::OutputMeshType *
  ShapeModelToImageRegistrationMethod< TShapeModel, TOutputMesh>::GetOutput() const
  {
    return static_cast< const TOutputMesh*>(this->ProcessObject::GetOutput(0));
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  itk::ModifiedTimeType ShapeModelToImageRegistrationMethod< TShapeModel, TOutputMesh>::GetMTime() const
  {
    itk::ModifiedTimeType mtime = Superclass::GetMTime();
    itk::ModifiedTimeType m;

    if (m_ShapeModel) {
      m = m_ShapeModel->GetMTime();
      mtime = (m > mtime ? m : mtime);
    }

    if (m_LevelSetImage) {
      m = m_LevelSetImage->GetMTime();
      mtime = (m > mtime ? m : mtime);
    }

    if (m_ShapeTransform) {
      m = m_ShapeTransform->GetMTime();
      mtime = (m > mtime ? m : mtime);
    }

    return mtime;
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::InitializeTransform()
  {
    // initialize spatial transform
    if (m_TransformInitializer != nullptr) {
      try {
        m_TransformInitializer->Update();
      }
      catch (itk::ExceptionObject& excep) {
        itkExceptionMacro(<< excep);
      }
      m_SpatialTransform = m_TransformInitializer->GetTransform();
      m_NumberOfSpatialParameters = m_SpatialTransform->GetNumberOfParameters();
      m_SpatialScales = m_TransformInitializer->GetScales();
    }
    else {
      m_NumberOfSpatialParameters = 0;
      m_SpatialScales.SetSize(m_NumberOfSpatialParameters);
    }

    // initialize shape model transform
    m_ShapeTransform = ShapeTransformType::New();
    m_ShapeTransform->SetStatisticalModel(m_ShapeModel);
    m_ShapeTransform->SetIdentity();
    m_NumberOfShapeModelParameters = m_ShapeTransform->GetNumberOfParameters();

    // initialize composite transform
    if (m_SpatialTransform != nullptr) {
      CompositeTransformType::Pointer transform = CompositeTransformType::New();
      transform->AddTransform(m_SpatialTransform);
      transform->AddTransform(m_ShapeTransform);
      m_Transform = transform;
    }
    else {
      m_Transform = m_ShapeTransform;
    }

    // initialize scales
    m_Scales.set_size(m_Transform->GetNumberOfParameters());
    size_t count = 0;

    for (size_t i = 0; i < m_NumberOfShapeModelParameters; ++i, ++count) {
      m_Scales[count] = 1 / m_ModelScale;
    }

    for (size_t i = 0; i < m_NumberOfSpatialParameters; ++i, ++count) {
      m_Scales[count] = 1 / m_SpatialScales[i];
    }
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::InitializeMetric()
  {
    m_Metric = MetricType::New();
    m_Metric->SetShapeModel(m_ShapeModel);
    m_Metric->SetImage(m_LevelSetImage);
    m_Metric->SetSpatialTransform(m_SpatialTransform);
    m_Metric->SetRegularizationParameter(m_RegularizationParameter);
    m_Metric->SetDegree(m_Degree);
    try {
      m_Metric->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      itkExceptionMacro(<< excep);
    }
  }
  //----------------------------------------------------------------------------
  // initialize optimizer
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::InitializeOptimizer()
  {
    if (m_NumberOfIterations < 1) {
      itkExceptionMacro(<< "number of iterations is zero");
    }

    m_Optimizer = itk::LBFGSOptimizer::New();
    m_Optimizer->SetCostFunction(m_Metric);
    m_Optimizer->SetScales(m_Scales);
    m_Optimizer->SetMaximumNumberOfFunctionEvaluations(m_NumberOfIterations);
    m_Optimizer->SetGradientConvergenceTolerance(m_GradientConvergenceTolerance);
    m_Optimizer->SetLineSearchAccuracy(m_LineSearchAccuracy);
    m_Optimizer->SetDefaultStepLength(m_DefaultStepLength);
    m_Optimizer->MinimizeOn();
    m_Optimizer->SetGlobalWarningDisplay(false);
    m_Optimizer->SetTrace(false);
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::GenerateData()
  {
    m_Clock.Start();

    // initialize transform, metric, optimizer and multi-stage data
    this->InitializeTransform();
    this->InitializeMetric();
    this->InitializeOptimizer();

    // run optimization
    m_Optimizer->SetInitialPosition(m_Transform->GetParameters());
    try {
      m_Optimizer->StartOptimization();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
    }

    // set results to transform
    m_Transform->SetParameters(m_Optimizer->GetCurrentPosition());

    // generate output data
    this->GenerateOutputData();

    m_Clock.Stop();
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::GenerateOutputData()
  {
    // compute transformed mesh
    typedef typename itk::TransformMeshFilter<DatasetType, TOutputMesh, TransformType> TransformMeshFilterType;
    typename TransformMeshFilterType::Pointer transform = TransformMeshFilterType::New();
    transform->SetInput(m_ShapeModel->GetRepresenter()->GetReference());
    transform->SetTransform(m_Transform);
    try {
      transform->Update();
    }
    catch (itk::ExceptionObject& excep) {
      itkExceptionMacro(<< excep);
    }

    // connect to the output
    this->SetNthOutput(0, transform->GetOutput());
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelToImageRegistrationMethod<TShapeModel, TOutputMesh>::PrintReport(std::ostream& os)
  {
    os << "shape model info" << std::endl;
    os << "number of components       " << m_ShapeModel->GetNumberOfPrincipalComponents() << std::endl;
    os << "number of points           " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
    os << "number of cells            " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
    os << std::endl;

    os << "metric info" << std::endl;
    os << "name of class              " << m_Metric->GetNameOfClass() << std::endl;
    os << "number of threads          " << m_Metric->GetNumberOfThreads() << " / " << m_Metric->GetMaximalNumberOfThreads() << std::endl;
    os << "regularization parameter   " << m_RegularizationParameter << std::endl;
    os << "degree                     " << m_Degree << std::endl;
    os << std::endl;

    os << "optimizer info" << std::endl;
    os << "stop condition description " << m_Optimizer->GetStopConditionDescription() << std::endl;
    os << "cost function value        " << m_Optimizer->GetValue() << std::endl;
    os << "number of evaluations      " << m_Optimizer->GetMaximumNumberOfFunctionEvaluations() << std::endl;
    os << "gradient tolerance         " << m_Optimizer->GetGradientConvergenceTolerance() << std::endl;
    os << "elapsed time, sec          " << m_Clock.GetTotal() << std::endl;
    os << std::endl;

    os << m_Transform->GetTransformTypeAsString() << ", " << m_Transform->GetTransformCategory() << std::endl;
    os << "number of parameters " << m_Transform->GetNumberOfParameters() << std::endl;
    os << std::endl;
  }
}
