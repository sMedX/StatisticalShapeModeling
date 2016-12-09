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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::InitializeOptimizer()
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::GenerateData()
  {
    m_Clock.Start();

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

    m_Clock.Stop();
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
  void ShapeModelRegistrationMethod<TShapeModel, TOutputMesh>::PrintReport(std::ostream& os)
  {
    os << "shape model info" << std::endl;
    os << "number of components       " << m_ShapeModel->GetNumberOfPrincipalComponents() << std::endl;
    os << "number of points           " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
    os << "number of cells            " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
    os << std::endl;

    os << "metric info" << std::endl;
    os << "name of class              " << m_Metric->GetNameOfClass() << std::endl;
    os << "regularization parameter   " << m_RegularizationParameter << std::endl;
    os << "degree                     " << m_Degree << std::endl;
    os << std::endl;

    os << "optimizer info" << std::endl;
    os << "stop condition description " << m_Optimizer->GetStopConditionDescription() << std::endl;
    os << "cost function value        " << m_Optimizer->GetValue() << std::endl;
    os << "elapsed time, sec          " << m_Clock.GetTotal() << std::endl;
    os << std::endl;

    os << m_ShapeTransform->GetTransformTypeAsString() << ", " << m_ShapeTransform->GetTransformCategory() << std::endl;
    os << "number of parameters " << m_ShapeTransform->GetNumberOfParameters() << std::endl;
    os << std::endl;
  }
}
