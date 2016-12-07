#pragma once

#include <itkTransformMeshFilter.h>
#include <itkLinearInterpolateImageFunction.h>

#include "itkShapeModelRegistrationMethod.h"

namespace itk
{
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::ShapeModelRegistrationMethod()
  {
    this->SetNumberOfRequiredInputs(0);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, TOutputMesh::New());

    m_ShapeModel = nullptr;
    m_LevelSetImage = nullptr;
    m_ShapeTransform = nullptr;
    m_Metric = nullptr;
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  const typename ShapeModelRegistrationMethod< TShapeModel, TOutputMesh>::OutputMeshType *
  ShapeModelRegistrationMethod< TShapeModel, TOutputMesh>::GetOutput() const
  {
    return static_cast< const TOutputMesh*>(this->ProcessObject::GetOutput(0));
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  itk::ModifiedTimeType ShapeModelRegistrationMethod< TShapeModel, TOutputMesh>::GetMTime() const
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::InitializeTransform()
  {
    int numberOfOptimizationParameters = m_ShapeModel->GetNumberOfPrincipalComponents();

    // initialize transforms of the model
    m_Scales.set_size(numberOfOptimizationParameters);
    m_Scales.Fill(1.0 / m_ModelScale);

    // shape model transform
    m_ShapeTransform = ShapeTransformType::New();
    m_ShapeTransform->SetStatisticalModel(m_ShapeModel);
    m_ShapeTransform->SetIdentity();
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::InitializeMetric()
  {
    m_Metric = MetricType::New();
    m_Metric->SetShapeModel(m_ShapeModel);
    m_Metric->SetImage(m_LevelSetImage);
    m_Metric->SetTransform(m_ShapeTransform);
    m_Metric->SetRegularizationParameter(m_RegularizationParameter);
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::InitializeOptimizer()
  {
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::GenerateData()
  {
    // initialize transform, metric, optimizer and multi-stage data
    this->InitializeTransform();
    this->InitializeMetric();
    this->InitializeOptimizer();

    // run optimization
    m_Optimizer->SetInitialPosition(m_ShapeTransform->GetParameters());
    try {
      m_Optimizer->StartOptimization();
    }
    catch (itk::ExceptionObject& excep) {
      itkExceptionMacro(<< excep);
    }

    // generate output data
    this->GenerateOutputData();
  }
  //----------------------------------------------------------------------------
  template <typename TShapeModel, typename TOutputMesh>
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::GenerateOutputData()
  {
    // compute transformed mesh
    typedef itk::TransformMeshFilter<DatasetType, TOutputMesh, ShapeTransformType> TransformFilterType;
    typename TransformFilterType::Pointer transform = TransformFilterType::New();
    transform->SetInput(m_ShapeModel->GetRepresenter()->GetReference());
    transform->SetTransform(m_ShapeTransform);
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::PrintReport()
  {
    std::cout << "shape model info" << std::endl;
    std::cout << "number of components       " << m_ShapeModel->GetNumberOfPrincipalComponents() << std::endl;
    std::cout << "number of points           " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
    std::cout << "number of cells            " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
    std::cout << std::endl;

    std::cout << "metric info" << std::endl;
    std::cout << "name of class              " << m_Metric->GetNameOfClass() << std::endl;
    std::cout << "regularization parameter   " << m_RegularizationParameter << std::endl;
    std::cout << std::endl;

    if (m_NumberOfIterations > 0) {
      std::cout << "optimizer info" << std::endl;
      std::cout << "stop condition description " << m_Optimizer->GetStopConditionDescription() << std::endl;
      std::cout << "cost function value        " << m_Optimizer->GetValue() << std::endl;
      std::cout << std::endl;
    }

    std::cout << m_ShapeTransform->GetTransformTypeAsString() << ", " << m_ShapeTransform->GetTransformCategory() << std::endl;
    std::cout << "number of parameters " << m_ShapeTransform->GetNumberOfParameters() << std::endl;
    std::cout << std::endl;
  }
}
