#pragma once

#include <itkCastImageFilter.h>

#include "ssmShapeModelToImageRegistrationMethod.h"
#include "ssmBinaryImageToLevelSetImageFilter.h"

namespace ssm
{
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::ShapeModelToImageRegistrationMethod()
{
  m_Image = nullptr;
  m_Metric = MetricType::New();

  m_LevelSetImage = nullptr;
  m_ComputeLevelSetImage = true;
}
//----------------------------------------------------------------------------
//
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::Initialize()
{
  if (!m_Image) {
    itkExceptionMacro(<< "Image is not initialized.");
  }

  Superclass::Initialize();

  this->ComputeLevelSet();
  this->InitializeTransform();
  this->InitializeMetric();
  this->InitializeOptimizer();
}
//----------------------------------------------------------------------------
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::ComputeLevelSet()
{
  // compute level set image
  if (m_ComputeLevelSetImage) {
    typedef ssm::BinaryImageToLevelSetImageFilter<InputImageType, FloatImageType> BinaryImageToLevelSetImageType;
    BinaryImageToLevelSetImageType::Pointer levelset = BinaryImageToLevelSetImageType::New();
    levelset->SetInput(m_Image);
    levelset->Update();
    m_LevelSetImage = levelset->GetOutput();
  }
  else {
    typedef itk::CastImageFilter<InputImageType, LevelSetImageType> CastFilterType;
    auto cast = CastFilterType::New();
    cast->SetInput(m_Image);
    cast->Update();
    m_LevelSetImage = cast->GetOutput();
  }
}
//----------------------------------------------------------------------------
// initialize metric
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::InitializeMetric()
{
  m_Metric->SetShapeModel(this->m_ShapeModel);
  m_Metric->SetLevelSetImage(this->m_LevelSetImage);
  m_Metric->SetTransform(this->m_Transform);
  m_Metric->SetLogger(this->m_Logger);
  try {
    m_Metric->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    this->ExceptionHandler(excep, __LINE__);
  }
  Superclass::SetMetric(m_Metric);
}
//----------------------------------------------------------------------------
//
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::GenerateData()
{
  this->m_Time.Start();

  // initialize the interconnects between components
  try {
    this->Initialize();
  }
  catch (itk::ExceptionObject & excep) {
    throw excep;
  }

  // perform registration
  this->PerformRegistration();

  // generate output data and print report
  this->GenerateOutputData();
}
//----------------------------------------------------------------------------
// print report
template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistrationMethod<TShapeModel, TInputImage, TOutputMesh>::PrintReport()
{
  m_Metric->PrintReport();
  Superclass::PrintReport();
}
}
