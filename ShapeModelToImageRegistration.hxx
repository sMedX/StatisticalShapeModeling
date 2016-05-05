#ifndef __ShapeModelToImageRegistration_hxx
#define __ShapeModelToImageRegistration_hxx

#include <itkEuler3DTransform.h>
#include <itkVersorTransform.h>
#include <itkVersorRigid3DTransform.h>
#include <itkScaleTransform.h>
#include <itkCompositeTransform.h>
#include <itkImageMomentsCalculator.h>
#include <itkLBFGSBOptimizer.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAbsImageFilter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkBoundingBox.h>
#include <itkPointSetToImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresPointSetToImageMetric.h>
#include <utils/itkPenalizingMeanSquaresPointSetToImageMetric.h>

#include "ShapeModelToImageRegistration.h"

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
ShapeModelToImageRegistration<TInputImage, TOutputMesh>::ShapeModelToImageRegistration()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(2);

  this->SetNthOutput(0, TOutputMesh::New());
  this->SetNthOutput(1, TOutputMesh::New());
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::GenerateData()
{
  m_Image = this->GetInput();
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
  m_Transform->SetAllTransformsToOptimizeOn();

  //perform geometrical registration
  m_Transform->SetNthTransformToOptimize(m_ModelTransformIndex, false);

  typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, FloatImageType> MetricType;
  typename MetricType::Pointer metric = MetricType::New();
  m_Metric = metric;

  this->PointSetToImageRegistration();

  //perform registration with a shape model
  m_Transform->SetNthTransformToOptimize(m_ModelTransformIndex, true);

  typedef itk::PenalizingMeanSquaresPointSetToImageMetric<PointSetType, FloatImageType> PenalizingMetricType;
  typename PenalizingMetricType::Pointer penalizingMetric = PenalizingMetricType::New();
  penalizingMetric->SetRegularizationParameter(m_RegularizationParameter);
  penalizingMetric->SetNumberOfModelComponents(m_ShapeModel->GetNumberOfPrincipalComponents());
  m_Metric = penalizingMetric;

  this->PointSetToImageRegistration();

  //generate output data
  this->GenerateOutputData();
}

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::PointSetToImageRegistration()
{
  if (m_NumberOfIterations < 1) {
    return;
  }

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
}
//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::ComputePotentialImage()
{
  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, FloatImageType> DistanceFilterType;
  DistanceFilterType::Pointer distanceMapFilter = DistanceFilterType::New();
  distanceMapFilter->SetUseImageSpacing(true);
  distanceMapFilter->SetInput(m_Image);

  typedef itk::AbsImageFilter<FloatImageType, FloatImageType> AbsImageFilterType;
  AbsImageFilterType::Pointer absImageFilter = AbsImageFilterType::New();
  absImageFilter->SetInput(distanceMapFilter->GetOutput());
  absImageFilter->Update();
  m_PotentialImage = absImageFilter->GetOutput();
}

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::InitializeTransform()
{
  // Compute a bounding box around the reference shape
  typedef  itk::BoundingBox<int, TOutputMesh::PointDimension, float, typename TOutputMesh::PointsContainer> BoundingBoxType;
  typename BoundingBoxType::Pointer boundingBox = BoundingBoxType::New();
  boundingBox->SetPoints(m_Surface->GetPoints());
  boundingBox->ComputeBoundingBox();

  BinaryImageType::SpacingType spacing = m_Image->GetSpacing();
  BinaryImageType::PointType origin = boundingBox->GetMinimum();
  BinaryImageType::SizeType size;

  for (int n = 0; n < ImageDimension; ++n) {
    origin[n] -= m_Margin;
    size[n] = (boundingBox->GetMaximum()[n] - boundingBox->GetMinimum()[n] + 2 * m_Margin) / spacing[n];
  }

  typedef itk::TriangleMeshToBinaryImageFilter<TOutputMesh, BinaryImageType> ShapeToBinaryImageFilterType;
  typename ShapeToBinaryImageFilterType::Pointer shapeToImage = ShapeToBinaryImageFilterType::New();
  shapeToImage->SetInput(m_Surface);
  shapeToImage->SetSize(size);
  shapeToImage->SetOrigin(origin);
  shapeToImage->SetSpacing(spacing);
  shapeToImage->SetDirection(m_Image->GetDirection());
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
  fixedCalculator->SetImage(m_Image);
  fixedCalculator->Compute();

  typename TransformType::InputPointType center;
  typename TransformType::OutputVectorType translation;

  for (unsigned int i = 0; i < ImageDimension; i++) {
    center[i] = movingCalculator->GetCenterOfGravity()[i];
    translation[i] = fixedCalculator->GetCenterOfGravity()[i] - movingCalculator->GetCenterOfGravity()[i];
  }

  //number of optimization parameters
  int numberOfRotationComponents = 3;
  int numberOfTranslationComponents = 3;
  int numberOfScalingComponents = 3;
  int numberOfModelComponents = m_ShapeModel->GetNumberOfPrincipalComponents();
  int numberOfOptimizationParameters = numberOfModelComponents + 15;

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
  typedef itk::ScaleTransform<double, ImageDimension> ScaleTransformType;
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

  //shape model transform
  m_ModelTransform = ModelTransformType::New();
  m_ModelTransform->SetStatisticalModel(m_ShapeModel);
  m_ModelTransform->SetIdentity();

  for (int i = 0; i < numberOfModelComponents; ++i, ++count) {
    m_Scales[count] = 1.0 / m_ModelScale;
  }

  // composite transform OutputPoint = rigid3DTransform2( scaleTransform ( rigid3DTransform1( m_ModelTransform(InputPoint))))
  m_Transform = TransformType::New();
  m_Transform->AddTransform(rigid3DTransform2);   //transform 0
  m_Transform->AddTransform(scaleTransform);      //transform 1
  m_Transform->AddTransform(rigid3DTransform1);   //transform 2

  m_ModelTransformIndex = 3;
  m_Transform->AddTransform(m_ModelTransform);    //transform 3
}

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::GenerateOutputData()
{
  //compute moved output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, TransformType> TransformFilterType;
  typename TransformFilterType::Pointer transform0 = TransformFilterType::New();
  transform0->SetInput(m_Surface);
  transform0->SetTransform(m_Transform);

  try {
    transform0->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  this->GraftNthOutput(0, transform0->GetOutput());

  //compute deformed output
  typedef itk::TransformMeshFilter<TOutputMesh, TOutputMesh, ModelTransformType> ModelTransformFilterType;
  typename ModelTransformFilterType::Pointer transform1 = ModelTransformFilterType::New();
  transform1->SetInput(m_Surface);
  transform1->SetTransform(m_ModelTransform);

  try {
    transform1->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

 this->GraftNthOutput(1, transform1->GetOutput());
}

//----------------------------------------------------------------------------
template <typename TInputImage, typename TOutputMesh>
void ShapeModelToImageRegistration<TInputImage, TOutputMesh>::PrintReport(std::ostream& os) const
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
