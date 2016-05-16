#ifndef __PointSetToImageRegistration_hxx
#define __PointSetToImageRegistration_hxx

#include <itkEuler3DTransform.h>
#include <itkScaleTransform.h>
#include <itkTransformMeshFilter.h>
#include <itkImageMomentsCalculator.h>
#include <itkLBFGSBOptimizer.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkBoundingBox.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresPointSetToImageMetric.h>
#include "PointSetToImageRegistration.h"

//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
PointSetToImageRegistration<TInputMesh, TOutputMesh>::PointSetToImageRegistration()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(1);
  this->SetNthOutput(0, TOutputMesh::New());
}
//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void PointSetToImageRegistration<TInputMesh, TOutputMesh>::GenerateData()
{
  m_Surface = const_cast <TInputMesh*> (this->GetInput());

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

  //define metric
  typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, PotentialImageType> MetricType;
  m_Metric = MetricType::New();

  //run point set to image registration
  typedef itk::PointSetToImageRegistrationMethod<PointSetType, FloatImageType> RegistrationFilterType;
  typename RegistrationFilterType::Pointer registration = RegistrationFilterType::New();
  registration->SetInitialTransformParameters(m_Transform->GetParameters());
  registration->SetMetric(m_Metric);
  registration->SetInterpolator(interpolator);
  registration->SetOptimizer(m_Optimizer);
  registration->SetTransform(m_Transform);
  registration->SetFixedPointSet(m_PointSet);
  registration->SetMovingImage(m_PotentialImage);

  try {
    registration->Update();
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
void PointSetToImageRegistration<TInputMesh, TOutputMesh>::ComputePotentialImage()
{
  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, PotentialImageType> DistanceFilterType;
  DistanceFilterType::Pointer distanceMapFilter = DistanceFilterType::New();
  distanceMapFilter->SetUseImageSpacing(true);
  distanceMapFilter->SetInput(m_Label);
  distanceMapFilter->Update();

  m_PotentialImage = distanceMapFilter->GetOutput();
}

//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void PointSetToImageRegistration<TInputMesh, TOutputMesh>::InitializeTransform()
{
  // Compute a bounding box around the input mesh
  typedef  itk::BoundingBox<int, TOutputMesh::PointDimension, float, typename TOutputMesh::PointsContainer> BoundingBoxType;
  typename BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
  boundingBox->SetPoints(m_Surface->GetPoints());
  boundingBox->ComputeBoundingBox();

  BinaryImageType::SpacingType spacing = m_Label->GetSpacing();
  BinaryImageType::PointType origin = boundingBox->GetMinimum();
  BinaryImageType::SizeType size;

  for (int n = 0; n < PointDimension; ++n) {
    origin[n] -= m_Margin;
    size[n] = (boundingBox->GetMaximum()[n] - boundingBox->GetMinimum()[n] + 2 * m_Margin) / spacing[n];
  }

  typedef itk::TriangleMeshToBinaryImageFilter<TOutputMesh, BinaryImageType> ShapeToBinaryImageFilterType;
  typename ShapeToBinaryImageFilterType::Pointer shapeToImage = ShapeToBinaryImageFilterType::New();
  shapeToImage->SetInput(m_Surface);
  shapeToImage->SetSize(size);
  shapeToImage->SetOrigin(origin);
  shapeToImage->SetSpacing(spacing);
  shapeToImage->SetDirection(m_Label->GetDirection());
  shapeToImage->SetOutsideValue(0);
  shapeToImage->SetInsideValue(1);

  try {
    shapeToImage->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  // moment calculators
  typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
  typedef typename ImageCalculatorType::VectorType VectorType;

  ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
  movingCalculator->SetImage(shapeToImage->GetOutput());
  movingCalculator->Compute();

  ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
  fixedCalculator->SetImage(m_Label);
  fixedCalculator->Compute();

  typename TransformType::InputPointType center;
  typename TransformType::OutputVectorType translation;

  for (unsigned int i = 0; i < PointDimension; i++) {
    center[i] = movingCalculator->GetCenterOfGravity()[i];
    translation[i] = fixedCalculator->GetCenterOfGravity()[i] - movingCalculator->GetCenterOfGravity()[i];
  }

  //number of optimization parameters
  int numberOfRotationComponents = 3;
  int numberOfTranslationComponents = 3;
  int numberOfScalingComponents = 3;
  int numberOfOptimizationParameters = 15;

  //initialize transforms of the model
  m_Scales.set_size(numberOfOptimizationParameters);
  size_t count = 0;

  //second Euler 3D transform
  typedef itk::Euler3DTransform<double> Euler3DTransformType;
  Euler3DTransformType::Pointer rigid3DTransform2 = Euler3DTransformType::New();
  rigid3DTransform2->SetIdentity();
  rigid3DTransform2->SetCenter(center);
  rigid3DTransform2->SetTranslation(translation);

  for (int i = 0; i < numberOfRotationComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_RotationScale;
  }

  for (int i = 0; i < numberOfTranslationComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_TranslationScale;
  }

  //scale transform
  typedef itk::ScaleTransform<double, PointDimension> ScaleTransformType;
  typename ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();
  scaleTransform->SetIdentity();
  scaleTransform->SetCenter(center);

  for (int i = 0; i < numberOfScalingComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_ScalingScale;
  }

  //first Euler 3D transform
  typedef itk::Euler3DTransform<double> VersorTransformType;
  VersorTransformType::Pointer rigid3DTransform1 = VersorTransformType::New();
  rigid3DTransform1->SetIdentity();
  rigid3DTransform1->SetCenter(center);

  for (int i = 0; i < numberOfRotationComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_RotationScale;
  }

  for (int i = 0; i < numberOfTranslationComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_TranslationScale;
  }

  // composite transform OutputPoint = rigid3DTransform2( scaleTransform ( rigid3DTransform1(InputPoint)))
  m_Transform = TransformType::New();
  m_Transform->AddTransform(rigid3DTransform2);   //transform 0
  m_Transform->AddTransform(scaleTransform);      //transform 1
  m_Transform->AddTransform(rigid3DTransform1);   //transform 2

  m_Transform->SetAllTransformsToOptimizeOn();
}
//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void PointSetToImageRegistration<TInputMesh, TOutputMesh>::GenerateOutputData()
{
  //compute moved output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, TransformType> TransformFilterType;
  typename TransformFilterType::Pointer transform = TransformFilterType::New();
  transform->SetInput(m_Surface);
  transform->SetTransform(m_Transform);

  try {
    transform->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  this->GraftNthOutput(0, transform->GetOutput());
}

//----------------------------------------------------------------------------
template <typename TInputMesh, typename TOutputMesh>
void PointSetToImageRegistration<TInputMesh, TOutputMesh>::PrintReport(std::ostream& os) const
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
