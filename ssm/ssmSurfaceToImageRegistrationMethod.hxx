#pragma once

#include <itkTranslationTransform.h>
#include <itkEuler3DTransform.h>
#include <itkScaleTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkTransformMeshFilter.h>
#include <itkImageMomentsCalculator.h>
#include <itkLBFGSBOptimizer.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkBoundingBox.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresPointSetToImageMetric.h>
#include <itkBinaryThresholdImageFilter.h>

#include "ssmSurfaceToImageRegistrationMethod.h"

namespace ssm
{
  //----------------------------------------------------------------------------
  template <typename TInputMesh, typename TOutputMesh>
  SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::SurfaceToImageRegistrationMethod()
  {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, TOutputMesh::New());
  }
  //----------------------------------------------------------------------------
  template <typename TInputMesh, typename TOutputMesh>
  void SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::GenerateData()
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

    this->InitializeTransform();

    // setup interpolator
    m_Interpolator = InterpolatorType::New();

    // setup metric
    typedef itk::MeanSquaresPointSetToImageMetric<PointSetType, LevelsetImageType> MetricType;
    m_Metric = MetricType::New();
    m_Metric->SetMovingImage(m_LevelsetImage);
    m_Metric->SetFixedPointSet(m_PointSet);
    m_Metric->SetTransform(m_Transform);
    m_Metric->SetInterpolator(m_Interpolator);
    try {
      m_Metric->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      itkExceptionMacro(<< excep);
    }

    // setup optimizer
    m_Optimizer = OptimizerType::New();
    m_Optimizer->SetCostFunction(m_Metric);
    m_Optimizer->SetInitialPosition(m_Transform->GetParameters());
    m_Optimizer->SetMaximumNumberOfFunctionEvaluations(m_NumberOfIterations);
    m_Optimizer->SetScales(m_Scales);
    m_Optimizer->SetGradientConvergenceTolerance(m_GradientConvergenceTolerance);
    m_Optimizer->SetLineSearchAccuracy(m_LineSearchAccuracy);
    m_Optimizer->SetDefaultStepLength(m_DefaultStepLength);
    m_Optimizer->MinimizeOn();
    try {
      m_Optimizer->StartOptimization();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      itkExceptionMacro(<< excep);
    }

    this->GenerateOutputData();
  }
  //----------------------------------------------------------------------------
  template <typename TInputMesh, typename TOutputMesh>
  void SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::ComputeLabelImage()
  {
    typedef itk::BinaryThresholdImageFilter <LevelsetImageType, BinaryImageType> BinaryThresholdImageFilterType;
    BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
    threshold->SetInput(m_LevelsetImage);
    threshold->SetLowerThreshold(std::numeric_limits<LevelsetImageType::PixelType>::lowest());
    threshold->SetUpperThreshold(0);
    threshold->SetInsideValue(1);
    threshold->SetOutsideValue(0);
    threshold->Update();

    m_Mask = threshold->GetOutput();
  }

  //----------------------------------------------------------------------------
  template <typename TInputMesh, typename TOutputMesh>
  void SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::InitializeTransform()
  {
    // compute label of the input level set image to initialize transform
    this->ComputeLabelImage();

    // Compute a bounding box of the input mesh
    typename TOutputMesh::BoundingBoxType::ConstPointer boundingBox = m_Surface->GetBoundingBox();

    BinaryImageType::SpacingType spacing = m_Mask->GetSpacing();
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
    shapeToImage->SetDirection(m_Mask->GetDirection());
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
    fixedCalculator->SetImage(m_Mask);
    fixedCalculator->Compute();

    typename TransformType::InputPointType center;
    typename TransformType::OutputVectorType translation;

    for (size_t i = 0; i < PointDimension; i++) {
      center[i] = movingCalculator->GetCenterOfGravity()[i];
      translation[i] = fixedCalculator->GetCenterOfGravity()[i] - movingCalculator->GetCenterOfGravity()[i];
    }

    // number of optimization parameters
    size_t numberOfTranslationComponents = 3;
    size_t numberOfRotationComponents = 3;
    size_t numberOfScalingComponents = 3;

    switch (m_TypeOfTransform) {
    case EnumTransformType::Translation: {
      // Translation transform
      typedef itk::TranslationTransform<double, PointDimension> TranslationTransformType;
      TranslationTransformType::Pointer translationTransform = TranslationTransformType::New();
      translationTransform->Translate(translation);
      m_Transform = translationTransform;

      //define scales
      m_Scales.set_size(m_Transform->GetNumberOfParameters());

      for (size_t i = 0; i < numberOfTranslationComponents; ++i) {
        m_Scales[i] = 1.0 / m_TranslationScale;
      }

      break;
    }
    case EnumTransformType::Euler3D: {
      // Euler3DTransform
      typedef itk::Euler3DTransform<double> Euler3DTransformType;
      Euler3DTransformType::Pointer Euler3DTransform = Euler3DTransformType::New();
      Euler3DTransform->SetIdentity();
      Euler3DTransform->SetCenter(center);
      Euler3DTransform->SetTranslation(translation);
      m_Transform = Euler3DTransform;

      //define scales
      m_Scales.set_size(m_Transform->GetNumberOfParameters());
      size_t count = 0;

      for (size_t i = 0; i < numberOfRotationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_RotationScale;
      }

      for (size_t i = 0; i < numberOfTranslationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_TranslationScale;
      }

      break;
    }
    case EnumTransformType::Similarity: {
      // Similarity3DTransform
      typedef itk::Similarity3DTransform<double> Similarity3DTransformType;
      Similarity3DTransformType::Pointer similarity3DTransform = Similarity3DTransformType::New();
      similarity3DTransform->SetIdentity();
      similarity3DTransform->SetCenter(center);
      similarity3DTransform->SetTranslation(translation);
      m_Transform = similarity3DTransform;

      //define scales
      m_Scales.set_size(m_Transform->GetNumberOfParameters());
      size_t count = 0;

      for (size_t i = 0; i < numberOfRotationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_RotationScale;
      }

      for (size_t i = 0; i < numberOfTranslationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_TranslationScale;
      }

      m_Scales[count] = 1.0 / m_ScalingScale;

      break;
    }
    case EnumTransformType::Affine: {
      //initialize transforms of the model
      int numberOfOptimizationParameters = 15;
      m_Scales.set_size(numberOfOptimizationParameters);
      size_t count = 0;

      //second Euler 3D transform
      typedef itk::Euler3DTransform<double> Euler3DTransformType;
      Euler3DTransformType::Pointer Euler3DTransform2 = Euler3DTransformType::New();
      Euler3DTransform2->SetIdentity();
      Euler3DTransform2->SetCenter(center);
      Euler3DTransform2->SetTranslation(translation);

      for (size_t i = 0; i < numberOfRotationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_RotationScale;
      }

      for (size_t i = 0; i < numberOfTranslationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_TranslationScale;
      }

      //scale transform
      typedef itk::ScaleTransform<double, PointDimension> ScaleTransformType;
      typename ScaleTransformType::Pointer scaleTransform = ScaleTransformType::New();
      scaleTransform->SetIdentity();
      scaleTransform->SetCenter(center);

      for (size_t i = 0; i < numberOfScalingComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_ScalingScale;
      }

      //first Euler 3D transform
      Euler3DTransformType::Pointer Euler3DTransform1 = Euler3DTransformType::New();
      Euler3DTransform1->SetIdentity();
      Euler3DTransform1->SetCenter(center);

      for (size_t i = 0; i < numberOfRotationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_RotationScale;
      }

      for (size_t i = 0; i < numberOfTranslationComponents; ++i, ++count) {
        m_Scales[count] = 1.0 / m_TranslationScale;
      }

      // composite transform OutputPoint = rigid3DTransform2( scaleTransform ( rigid3DTransform1(InputPoint)))
      m_CompositeTransform = CompositeTransformType::New();
      m_CompositeTransform->AddTransform(Euler3DTransform2);   //transform 0
      m_CompositeTransform->AddTransform(scaleTransform);      //transform 1
      m_CompositeTransform->AddTransform(Euler3DTransform1);   //transform 2
      m_CompositeTransform->SetAllTransformsToOptimizeOn();
      m_Transform = m_CompositeTransform;
    }
    default:
      throw itk::ExceptionObject("Invalid transform type");
    }
  }
  //----------------------------------------------------------------------------
  template <typename TInputMesh, typename TOutputMesh>
  void SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::GenerateOutputData()
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
  void SurfaceToImageRegistrationMethod<TInputMesh, TOutputMesh>::PrintReport(std::ostream& os) const
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
    os << "scales " << m_Scales << std::endl;

    if (m_TypeOfTransform != EnumTransformType::Affine) {
      os << std::endl;
      os << "1) " << m_Transform->GetTransformTypeAsString() << ", number of parameters ";
      os << m_Transform->GetNumberOfParameters() << std::endl;
      os << "transform category " << m_Transform->GetTransformCategory() << std::endl;
      os << m_Transform->GetParameters() << std::endl;
    }
    else {
      for (int n = 0; n < m_CompositeTransform->GetNumberOfTransforms(); ++n) {
        os << std::endl;
        os << n + 1 << ") " << m_CompositeTransform->GetNthTransformConstPointer(n)->GetTransformTypeAsString() << ", number of parameters ";
        os << m_CompositeTransform->GetNthTransformConstPointer(n)->GetNumberOfParameters() << std::endl;
        os << "transform category " << m_CompositeTransform->GetNthTransformConstPointer(n)->GetTransformCategory() << std::endl;
        os << m_CompositeTransform->GetNthTransformConstPointer(n)->GetParameters() << std::endl;
      }
    }
    os << std::endl;
  }
}
