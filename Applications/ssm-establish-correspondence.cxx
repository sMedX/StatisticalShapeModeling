﻿#include <boost/lexical_cast.hpp>

#include <itkImageMomentsCalculator.h>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmShapeModelToImageRegistrationMethod.h"
#include "ssmMeshToLevelSetImageFilter.h"
#include "ssmInitializeSpatialTransform.h"
#include "ssmCorrespondenceOptions.h"

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer model, const ssm::CorrespondenceOptions & options);

int main(int argc, char** argv)
{
  ssm::CorrespondenceOptions options;
  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  // read options from config file
  if (!options.ParseConfigFile()) {
    return EXIT_FAILURE;
  }
  options.PrintConfig();

  //----------------------------------------------------------------------------
  // read the reference surface
  MeshType::Pointer reference = MeshType::New();
  if (!readMesh<MeshType>(reference, options.GetReferenceFileName())) {
    return EXIT_FAILURE;
  }
  printMeshInfo<MeshType>(reference, options.GetReferenceFileName());

  //----------------------------------------------------------------------------
  // read list of files
  StringVector listOfInputFiles;
  try {
    listOfInputFiles = readListFromFile(options.GetInputList());
  }
  catch (ifstream::failure & e) {
    cerr << "could not read the data-list: " << e.what() << endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;

  for (const auto & fileName : listOfInputFiles) {
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }

    vectorOfSurfaces.push_back(surface);
    printMeshInfo<MeshType>(surface, fileName);
  }

  size_t numberOfSurfaces = vectorOfSurfaces.size();

  //----------------------------------------------------------------------------
  // perform establishing of correspondences
  itk::TimeProbe clock;
  clock.Start();

  StringVector listOfOutputFiles;

  for (size_t stage = 0; stage < options.GetNumberOfStages(); ++stage) {
    std::cout << "establish correspondence stage (" << stage + 1 << " / " << options.GetNumberOfStages() << ")" << std::endl;

    // build GP model for the reference surface
    StatisticalModelType::Pointer model = buildGPModel(reference, options.GetParameters()[0], options.GetScale(), options.GetNumberOfComponents()[0]);

    // initialize reference by zero
    for (MeshType::PointsContainerIterator iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
      iter.Value().Fill(0);
    }

    // preform GP model to surface registration
    for (size_t count = 0; count < numberOfSurfaces; ++count) {
      std::cout << "surface (" << count+1 << " / " << numberOfSurfaces << ")" << std::endl;

      MeshType::Pointer surface = vectorOfSurfaces[count];
      MeshType::Pointer output = shapeModelToSurfaceRegistration(surface, model, options);

      // add output to reference
      for (MeshType::PointsContainerIterator iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
        MeshType::PointType point = output->GetPoint(iter.Index());
        for (size_t d = 0; d < MeshType::PointDimension; ++d) {
          iter.Value()[d] += point[d] / numberOfSurfaces;
        }
      }

      // compute metrics and write output surface to file
      if (stage + 1 == options.GetNumberOfStages()) {
        // compute metrics
        typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilter;
        MeshToLevelSetImageFilter::Pointer levelset = MeshToLevelSetImageFilter::New();
        levelset->SetInput(surface);
        try {
          levelset->Update();
        }
        catch (itk::ExceptionObject& excep) {
          std::cerr << excep << std::endl;
          return EXIT_FAILURE;
        }

        typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
        typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
        PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
        metrics->SetPointSetAsMesh<MeshType>(output);
        metrics->SetImage(levelset->GetOutput());
        metrics->Compute();
        metrics->PrintReport();

        // print report to *.csv file
        std::cout << "print report to the file: " << options.GetReportFileName() << std::endl;
        std::cout << std::endl;
        metrics->PrintReportToFile(options.GetReportFileName(), getBaseNameFromPath(listOfInputFiles[count]));

        // write output surface to file
        auto fileName = options.FormatOutput(listOfInputFiles[count]);
        listOfOutputFiles.push_back(fileName);

        printMeshInfo<MeshType>(output, fileName);
        if (!writeMesh<MeshType>(output, fileName)) {
          EXIT_FAILURE;
        }
      }
    }

    // write the reference surface to the file
    std::string fileName = addFileNameSuffix(options.GetReferenceFileName(), "-" + std::to_string(stage + 1));
    if (!writeMesh<MeshType>(reference, fileName)) {
      EXIT_FAILURE;
    }
  }

  // write list of output files
  try {
    writeListToFile(options.GetOutputList(), listOfOutputFiles);
  }
  catch (std::ofstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  clock.Stop();
  std::cout << "elapsed time, sec " << clock.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer initialModel, const ssm::CorrespondenceOptions & options)
{
  // compute level set image
  typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilter;
  MeshToLevelSetImageFilter::Pointer levelset = MeshToLevelSetImageFilter::New();
  levelset->SetInput(surface);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  // shape model to image registration
  StatisticalModelType::Pointer model = initialModel;
  size_t numberOfStages = options.GetParameters().size();

  for (size_t stage = 0; stage < numberOfStages; ++stage) {
    std::cout << "registration stage (" << stage+1 << " / " << numberOfStages << ")" << std::endl;

    // if stage > 0 update GP model
    if (stage > 0) {
      model = buildGPModel(surface, options.GetParameters()[stage], options.GetScale(), options.GetNumberOfComponents()[stage]);
    }

    // initialize spatial transform
    MeshType::BoundingBoxType::ConstPointer boundingBox = model->DrawMean()->GetBoundingBox();
    BinaryImageType::SpacingType spacing(1);
    BinaryImageType::PointType origin = boundingBox->GetMinimum();
    BinaryImageType::SizeType size;
    for (size_t i = 0; i < Dimension; ++i) {
      size[i] = (boundingBox->GetMaximum()[i] - boundingBox->GetMinimum()[i]) / spacing[i];
    }

    typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> ShapeToBinaryImageFilterType;
    auto shapeToImage = ShapeToBinaryImageFilterType::New();
    shapeToImage->SetInput(model->DrawMean());
    shapeToImage->SetSize(size);
    shapeToImage->SetOrigin(origin);
    shapeToImage->SetSpacing(spacing);
    shapeToImage->SetOutsideValue(0);
    shapeToImage->SetInsideValue(1);
    try {
      shapeToImage->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      throw;
    }

    // moment calculators
    typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
    auto movingCalculator = ImageCalculatorType::New();
    movingCalculator->SetImage(shapeToImage->GetOutput());
    movingCalculator->Compute();

    auto fixedCalculator = ImageCalculatorType::New();
    fixedCalculator->SetImage(levelset->GetMask());
    fixedCalculator->Compute();

    typedef ImageCalculatorType::VectorType VectorType;
    VectorType center = movingCalculator->GetCenterOfGravity();
    VectorType translation = fixedCalculator->GetCenterOfGravity() - movingCalculator->GetCenterOfGravity();

    typedef ssm::InitializeSpatialTransform<double> TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTransformType(options.GetTransform());
    initializer->SetCenter(center);
    initializer->SetTranslation(translation);
    try {
      initializer->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      throw;
    }
    initializer->PrintReport();

    // perform registration
    typedef ssm::ShapeModelToImageRegistrationMethod<StatisticalModelType, FloatImageType, MeshType> ShapeModelRegistrationMethod;
    auto shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetComputeLevelSetImage(false);
    shapeModelToSurfaceRegistration->SetNumberOfIterations(options.GetNumberOfIterations());
    shapeModelToSurfaceRegistration->SetSpatialTransform(initializer->GetTransform());
    shapeModelToSurfaceRegistration->SetSpatialScales(initializer->GetScales());
    shapeModelToSurfaceRegistration->GetMetric()->SetRegularizationParameter(options.GetRegularization()[stage]);
    shapeModelToSurfaceRegistration->GetMetric()->SetDegree(2);
    try {
      shapeModelToSurfaceRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }
    shapeModelToSurfaceRegistration->PrintReport();

    surface = const_cast<MeshType*> (shapeModelToSurfaceRegistration->GetOutput());
  }

  return surface;
}

//==============================================================================
// Build Gaussian process model
StatisticalModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions)
{
  itk::TimeProbe clock;
  clock.Start();

  // create kernel
  typedef DataTypeShape::PointType PointType;
  auto gaussianKernel = new GaussianKernel<PointType>(parameters);

  typedef std::shared_ptr<const statismo::ScalarValuedKernel<PointType>> MatrixPointerType;
  MatrixPointerType kernel(gaussianKernel);

  typedef std::shared_ptr<statismo::MatrixValuedKernel<PointType>> KernelPointerType;
  KernelPointerType unscaledKernel(new statismo::UncorrelatedMatrixValuedKernel<PointType>(kernel.get(), Dimension));
  KernelPointerType modelBuildingKernel(new statismo::ScaledKernel<PointType>(unscaledKernel.get(), scale));

  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  typedef RepresenterType::DatasetPointerType DatasetPointerType;
  RepresenterType::Pointer representer = RepresenterType::New();
  representer->SetReference(surface);

  // build model
  typedef itk::LowRankGPModelBuilder<MeshType> ModelBuilderType;
  ModelBuilderType::Pointer modelBuilder = ModelBuilderType::New();
  modelBuilder->SetRepresenter(representer);

  // build and save model to file
  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model;
  try {
    model = modelBuilder->BuildNewModel(representer->IdentitySample(), *modelBuildingKernel.get(), numberOfBasisFunctions);
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  clock.Stop();
  std::cout << "Gaussian process model" << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << "parameters           " << "(" << parameters << ", " << scale << ")" << std::endl;
  std::cout << "elapsed time, sec    " << clock.GetTotal() << std::endl;
  std::cout << std::endl;

  return model;
}
