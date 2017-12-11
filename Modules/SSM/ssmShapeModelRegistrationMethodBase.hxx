#pragma once

#include <itkTransformMeshFilter.h>

#include "ssmShapeModelRegistrationMethodBase.h"

namespace ssm
{
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TOutputMesh>
ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::ShapeModelRegistrationMethodBase()
{
  this->SetNumberOfRequiredInputs(0);
  this->SetNumberOfRequiredOutputs(1);
  this->SetNthOutput(0, TOutputMesh::New());

  m_ShapeModel = nullptr;

  m_SpatialTransform = nullptr;
  m_Transform = nullptr;

  m_ModelScale = 3;
  m_NumberOfIterations = 500;

  m_UseNumberOfEpochs = true;
  m_NumberOfEpochs = 1;
  m_Step = 0;

  m_Optimizer = OptimizerType::New();
  m_Optimizer->SetMaximumNumberOfFunctionEvaluations(m_NumberOfIterations);
  m_Optimizer->SetLineSearchAccuracy(0.1);
  m_Optimizer->SetDefaultStepLength(0.5);
  m_Optimizer->SetGradientConvergenceTolerance(1.0e-05);
  m_Optimizer->SetGlobalWarningDisplay(false);
  m_Optimizer->SetTrace(false);
  m_Optimizer->MinimizeOn();

  m_ElapsedTime = 0;
  m_Logger = itk::Logger::New();
}
//----------------------------------------------------------------------------
//
template <typename TShapeModel, typename TOutputMesh>
const typename ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh >::OutputMeshType *
ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::GetOutput() const
{
  return static_cast<const OutputMeshType*>(this->ProcessObject::GetOutput(0));
}
//----------------------------------------------------------------------------
//
template <typename TShapeModel, typename TOutputMesh>
itk::ModifiedTimeType ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::GetMTime() const
{
  itk::ModifiedTimeType mtime = Superclass::GetMTime();
  itk::ModifiedTimeType m;

  if (m_ShapeModel) {
    m = m_ShapeModel->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  if (m_SpatialTransform) {
    m = m_SpatialTransform->GetMTime();
    mtime = (m > mtime ? m : mtime);
  }

  return mtime;
}
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::Initialize()
{
  if (!m_ShapeModel) {
    itkExceptionMacro(<< "Shape model is not initialized.");
  }

  if (!m_SpatialTransform) {
    itkExceptionMacro(<< "Spatial transform does not initialized.");
  }

  if (m_SpatialTransform->GetNumberOfParameters() != m_SpatialScales.size()) {
    itkExceptionMacro(<< "the number of parameters of spatial transform and size of spatial scales are not equal");
  }
}
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::InitializeTransform()
{
  // shape model multi transform
  m_Transform = ShapeModelMultiTransformType::New();
  m_Transform->SetShapeModel(m_ShapeModel);
  m_Transform->SetSpatialTransform(m_SpatialTransform);

  // define scales and bounds
  m_Scales.set_size(m_Transform->GetNumberOfParameters());

  size_t count = 0;

  for (size_t n = 0; n < m_ShapeModel->GetNumberOfPrincipalComponents(); ++n, ++count) {
    m_Scales[count] = m_ModelScale;
  }

  for (size_t n = 0; n < m_SpatialScales.Size(); ++n, ++count) {
    m_Scales[count] = m_SpatialScales[n];
  }

  // define inverse scales
  m_InverseScales.set_size(m_Scales.Size());

  for (size_t n = 0; n < m_Scales.Size(); ++n) {
    m_InverseScales[n] = 1 / m_Scales[n];
  }

  m_InitialParameters = m_Transform->GetParameters();
}
//----------------------------------------------------------------------------
// initialize optimizer
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::InitializeOptimizer()
{
  m_Optimizer->SetCostFunction(m_Metric);
  m_Optimizer->SetInitialPosition(m_Transform->GetParameters());
  m_Optimizer->SetScales(m_InverseScales);
  m_Optimizer->SetMaximumNumberOfFunctionEvaluations(m_NumberOfIterations);
}
//----------------------------------------------------------------------------
// perform registration
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::PerformRegistration()
{
  // compute multi-stage data
  this->ComputeMultiStageData();

  // perform multi-stage registration
  for (size_t n = 0; n < m_NumberOfUsedComponents.GetSize(); ++n) {
    m_Transform->SetNumberOfUsedComponents(m_NumberOfUsedComponents[n]);
    m_Optimizer->SetInitialPosition(m_Transform->GetParameters());
    try {
      m_Optimizer->StartOptimization();
    }
    catch (itk::ExceptionObject& excep) {
      m_Message.str("");
      m_Message << excep << std::endl;
      m_Logger->Fatal(m_Message.str());

      m_CostFunctionValues[n] = m_Optimizer->GetValue();
      m_Transform->SetParameters(m_Optimizer->GetCurrentPosition());
      break;
    }

    m_CostFunctionValues[n] = m_Optimizer->GetValue();
    m_Transform->SetParameters(m_Optimizer->GetCurrentPosition());
  }
}
//----------------------------------------------------------------------------
// initialize data for multi-stage optimization
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::ComputeMultiStageData()
{
  if (m_UseNumberOfEpochs) {
    if (m_NumberOfEpochs == 0) {
      m_Step = 1;
      m_NumberOfEpochs = m_ShapeModel->GetNumberOfPrincipalComponents();
    }
    else if (m_NumberOfEpochs > 1) {
      m_Step = itk::Math::Ceil<unsigned int>(double(m_ShapeModel->GetNumberOfPrincipalComponents()) / (m_NumberOfEpochs - 1));
      m_NumberOfEpochs = itk::Math::Ceil<unsigned int>(double(m_ShapeModel->GetNumberOfPrincipalComponents()) / m_Step) + 1;
    }
  }
  else {
    if (m_Step == 0) {
      m_NumberOfEpochs = 1;
    }
    else {
      m_NumberOfEpochs = itk::Math::Ceil<unsigned int>(double(m_ShapeModel->GetNumberOfPrincipalComponents()) / m_Step) + 1;
    }
  }

  m_NumberOfUsedComponents.SetSize(m_NumberOfEpochs);

  for (size_t n = 0; n < (m_NumberOfEpochs - 1); ++n) {
    m_NumberOfUsedComponents[n] = n * m_Step;
  }
  m_NumberOfUsedComponents[m_NumberOfEpochs - 1] = m_ShapeModel->GetNumberOfPrincipalComponents();

  m_CostFunctionValues.SetSize(m_NumberOfEpochs);
  m_CostFunctionValues.Fill(NAN);
}
//----------------------------------------------------------------------------
//
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::GenerateOutputData()
{
  // compute transformed mesh
  typedef itk::TransformMeshFilter<DatasetType, TOutputMesh, ShapeModelMultiTransformType> TransformFilterType;
  auto transform = TransformFilterType::New();
  transform->SetInput(m_ShapeModel->GetRepresenter()->GetReference());
  transform->SetTransform(m_Transform);
  try {
    transform->Update();
  }
  catch (itk::ExceptionObject & excep) {
    this->ExceptionHandler(excep, __LINE__);
  }

  // connect to the output
  this->SetNthOutput(0, transform->GetOutput());

  m_Time.Stop();
  m_ElapsedTime = m_Time.GetTotal();
}
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::PrintReport()
{
  std::cout << this->GetNameOfClass() << std::endl;
  std::cout << std::endl;

  std::cout << "shape model info" << std::endl;
  std::cout << "number of components       " << m_ShapeModel->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points           " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << "number of cells            " << m_ShapeModel->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
  std::cout << std::endl;

  if (m_NumberOfIterations > 0) {
    std::cout << "stop condition description " << m_Optimizer->GetStopConditionDescription() << std::endl;
    std::cout << "number of used components  " << m_NumberOfUsedComponents << std::endl;
    std::cout << "cost function values       " << m_CostFunctionValues << std::endl;
    std::cout << std::endl;
  }

  std::cout << m_Transform->GetTransformTypeAsString() << ", " << m_Transform->GetTransformCategory() << std::endl;
  std::cout << "number of parameters " << m_Transform->GetNumberOfParameters() << std::endl;
  std::cout << std::endl;
  std::cout << "initial position" << std::endl << m_InitialParameters << std::endl;
  std::cout << "final position" << std::endl << m_Transform->GetParameters() << std::endl;
  std::cout << std::endl;
  std::cout << "Elapsed time " << m_ElapsedTime << std::endl;
  std::cout << std::endl;
}
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TOutputMesh>
void ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>::ExceptionHandler(itk::ExceptionObject & excep, const size_t & line)
{
  std::cout << "Error: " << this->GetNameOfClass() << std::endl;
  std::cout << "File: " << __FILE__ << ", line " << line << std::endl;
  std::cout << excep << std::endl;
}
//----------------------------------------------------------------------------
}
