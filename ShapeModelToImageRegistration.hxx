#ifndef __ShapeModelToImageRegistration_hxx
#define __ShapeModelToImageRegistration_hxx

#include <itkLBFGSBOptimizer.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>

#include "ShapeModelToImageRegistration.h"

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::ShapeModelToImageRegistration()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(2);

  this->SetNthOutput(0, TOutputMesh::New());
  this->SetNthOutput(1, TOutputMesh::New());

  m_InputTransform = nullptr;
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::GenerateData()
{
  m_Surface = m_ShapeModel->GetRepresenter()->GetReference();

  //define point set
  m_PointSet = PointSetType::New();
  m_PointSet->SetPoints(m_Surface->GetPoints());

  typename PointSetType::PointDataContainer::Pointer pointData = PointSetType::PointDataContainer::New();
  pointData->Reserve(m_PointSet->GetNumberOfPoints());

  for (auto it = pointData->Begin(); it != pointData->End(); ++it) {
    it->Value() = 0;
  }
  m_PointSet->SetPointData(pointData);

  this->ComputePotentialImage();
  this->InitializeTransform();

  m_Metric = MetricType::New();
  m_Metric->SetRegularizationParameter(m_RegularizationParameter);
  m_Metric->SetNumberOfModelComponents(m_ShapeModel->GetNumberOfPrincipalComponents());

  //define interpolator
  typedef itk::LinearInterpolateImageFunction<FloatImageType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  //define optimizer
  m_Optimizer = OptimizerType::New();
  m_Optimizer->SetMaximumNumberOfFunctionEvaluations(m_NumberOfIterations);
  m_Optimizer->SetScales(m_Scales);
  m_Optimizer->SetGradientConvergenceTolerance(m_GradientConvergenceTolerance);
  m_Optimizer->SetLineSearchAccuracy(m_LineSearchAccuracy);
  m_Optimizer->SetDefaultStepLength(m_DefaultStepLength);
  m_Optimizer->MinimizeOn();

  //run registration with a shape model
  //point set to image registration method
  typedef itk::PointSetToImageRegistrationMethod<PointSetType, FloatImageType> RegistrationFilterType;
  typename RegistrationFilterType::Pointer m_Registration = RegistrationFilterType::New();
  m_Registration->SetInitialTransformParameters(m_Transform->GetParameters());
  m_Registration->SetMetric(m_Metric);
  m_Registration->SetInterpolator(interpolator);
  m_Registration->SetOptimizer(m_Optimizer);
  m_Registration->SetTransform(m_Transform);
  m_Registration->SetFixedPointSet(m_PointSet);
  m_Registration->SetMovingImage(m_PotentialImage);

  try {
    m_Registration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  //generate output data
  this->GenerateOutputData();
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::ComputePotentialImage()
{
  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, PotentialImageType> DistanceFilterType;
  DistanceFilterType::Pointer distanceMapFilter = DistanceFilterType::New();
  distanceMapFilter->SetUseImageSpacing(true);
  distanceMapFilter->SetInput(this->GetInput());
  distanceMapFilter->Update();
  m_PotentialImage = distanceMapFilter->GetOutput();
}

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::InitializeTransform()
{
  //TInputTransformType::Pointer inputTransform = const_cast<TInputTransformType *>(this->GetInputTransform());
  m_InputTransform = this->GetInputTransform();

  int numberOfInputParameters = m_InputTransform->GetNumberOfParameters();
  int numberOfModelComponents = m_ShapeModel->GetNumberOfPrincipalComponents();
  int numberOfOptimizationParameters = numberOfInputParameters + numberOfModelComponents;

  //initialize transforms of the model
  m_Scales.set_size(numberOfOptimizationParameters);
  m_Scales.Fill(1);

  //shape model transform
  m_ModelTransform = ModelTransformType::New();
  m_ModelTransform->SetStatisticalModel(m_ShapeModel);
  m_ModelTransform->SetIdentity();

  for (int i = 0; i < numberOfModelComponents; ++i) {
    m_Scales[numberOfInputParameters+i] = 1.0 / m_ModelScale;
  }

  // composite transform OutputPoint = m_InputTransform( m_ModelTransform(InputPoint) )
  m_Transform = TransformType::New();

  //first transform in queue
  for (int n = 0; n < m_InputTransform->GetNumberOfTransforms(); ++n) {
    m_Transform->AddTransform(m_InputTransform->GetNthTransform(n));
    m_Transform->SetNthTransformToOptimize(0, false);
  }

  //second transform in queue
  m_Transform->AddTransform(m_ModelTransform);
  m_Transform->SetNthTransformToOptimize(0, true);
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::GenerateOutputData()
{
  //compute output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, ModelTransformType> ModelTransformFilterType;
  typename ModelTransformFilterType::Pointer modelTransformMesh = ModelTransformFilterType::New();
  modelTransformMesh->SetInput(m_Surface);
  modelTransformMesh->SetTransform(m_ModelTransform);

  try {
    modelTransformMesh->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  this->GraftNthOutput(0, modelTransformMesh->GetOutput());

  //compute moved output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, TransformType> TransformFilterType;
  typename TransformFilterType::Pointer transformMesh = TransformFilterType::New();
  transformMesh->SetInput(m_Surface);
  transformMesh->SetTransform(m_Transform);

  try {
    transformMesh->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  this->GraftNthOutput(1, transformMesh->GetOutput());
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh, typename TInputTransformType>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh, TInputTransformType>::PrintReport(std::ostream& os) const
{
  os << this->GetNameOfClass() << std::endl;
  os << std::endl;

  if (m_NumberOfIterations > 0) {
    os << "      stop condition description " << m_Optimizer->GetStopConditionDescription() << std::endl;
    os << "            line search accuracy " << m_Optimizer->GetLineSearchAccuracy() << std::endl;
    os << "  gradient convergence tolerance " << m_Optimizer->GetGradientConvergenceTolerance() << std::endl;
    os << "             default step length " << m_Optimizer->GetDefaultStepLength() << std::endl;
    os << "  number of function evaluations " << m_Optimizer->GetMaximumNumberOfFunctionEvaluations() << std::endl;
    os << "             cost function value " << m_Optimizer->GetValue() << std::endl;
    os << std::endl;
  }

  os << m_Transform->GetTransformTypeAsString() << std::endl;
  os << "scales" << std::endl;
  os << m_Scales << std::endl;

  for (int n = 0; n < m_Transform->GetNumberOfTransforms(); ++n) {
    os << std::endl;
    os << n << ") " << m_Transform->GetNthTransformConstPointer(n)->GetTransformTypeAsString() << ", number of parameters ";
    os << m_Transform->GetNthTransformConstPointer(n)->GetNumberOfParameters() << std::endl;
    os << "transform category " << m_Transform->GetNthTransformConstPointer(n)->GetTransformCategory() << std::endl;
    os << m_Transform->GetNthTransformConstPointer(n)->GetParameters() << std::endl;
  }
}

//----------------------------------------------------------------------------
#endif
