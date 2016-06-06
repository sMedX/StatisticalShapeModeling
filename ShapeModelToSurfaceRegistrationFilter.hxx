#ifndef __ShapeModelToSurfaceRegistrationFilter_hxx
#define __ShapeModelToSurfaceRegistrationFilter_hxx

#include <itkLBFGSBOptimizer.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkBoundingBox.h>

#include "SurfaceToLevelSetImageFilter.h"
#include "ShapeModelToSurfaceRegistrationFilter.h"

//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
ShapeModelToSurfaceRegistrationFilter<TInputMesh, TOutputMesh>::ShapeModelToSurfaceRegistrationFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);

  m_Spacing.Fill(1);
}
//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void ShapeModelToSurfaceRegistrationFilter<TInputMesh, TOutputMesh>::GenerateData()
{
  // compute level set image
  m_Surface = const_cast <TInputMesh*> (this->GetInput());

  typedef SurfaceToLevelSetImageFilter<TInputMesh, LevelsetImageType> SurfaceToPotentialImageFilterType;
  SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
  levelset->SetMargin(m_Margin);
  levelset->SetSpacing(m_Spacing);
  levelset->SetInput(m_Surface);

  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  m_LevelsetImage = levelset->GetOutput();

  //define point set
  m_PointSet = PointSetType::New();
  m_PointSet->SetPoints(m_ShapeModel->GetRepresenter()->GetReference()->GetPoints());

  typename PointSetType::PointDataContainer::Pointer pointData = PointSetType::PointDataContainer::New();
  pointData->Reserve(m_PointSet->GetNumberOfPoints());

  for (auto it = pointData->Begin(); it != pointData->End(); ++it) {
    it->Value() = 0;
  }
  m_PointSet->SetPointData(pointData);

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
  m_Registration->SetMovingImage(m_LevelsetImage);

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
template <typename TInputMesh, typename TOutputMesh>
void ShapeModelToSurfaceRegistrationFilter<TInputMesh, TOutputMesh>::InitializeTransform()
{
  int numberOfOptimizationParameters = m_ShapeModel->GetNumberOfPrincipalComponents();

  //initialize transforms of the model
  m_Scales.set_size(numberOfOptimizationParameters);
  m_Scales.Fill(1.0 / m_ModelScale);

  //shape model transform
  m_Transform = ModelTransformType::New();
  m_Transform->SetStatisticalModel(m_ShapeModel);
  m_Transform->SetIdentity();
}
//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void ShapeModelToSurfaceRegistrationFilter<TInputMesh, TOutputMesh>::GenerateOutputData()
{
  //compute output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, ModelTransformType> ModelTransformFilterType;
  typename ModelTransformFilterType::Pointer transformModelSurface = ModelTransformFilterType::New();
  transformModelSurface->SetInput(m_ShapeModel->GetRepresenter()->GetReference());
  transformModelSurface->SetTransform(m_Transform);

  try {
    transformModelSurface->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  this->GraftOutput(transformModelSurface->GetOutput());
}
//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void ShapeModelToSurfaceRegistrationFilter<TInputMesh, TOutputMesh>::PrintReport(std::ostream& os) const
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

  os << std::endl;
  os << m_Transform->GetTransformTypeAsString() << ", number of parameters ";
  os << m_Transform->GetNumberOfParameters() << std::endl;
  os << "transform category " << m_Transform->GetTransformCategory() << std::endl;
  os << m_Transform->GetParameters() << std::endl;
}

//----------------------------------------------------------------------------
#endif
