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

ShapeModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, ShapeModelType::Pointer model, const ssm::CorrespondenceOptions & options);

int main(int argc, char** argv)
{
  // parse options
  ssm::CorrespondenceOptions options;

  if (!options.ParseOptions(argc, argv)) {
    return EXIT_FAILURE;
  }

  if (!checkFileName(options.GetReportFileName())) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // read the reference surface
  auto reference = MeshType::New();
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
  catch (std::ifstream::failure & e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;

  for (const auto & fileName : listOfInputFiles) {
    auto surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }

    vectorOfSurfaces.push_back(surface);
    printMeshInfo<MeshType>(surface, fileName);
  }

  size_t numberOfSurfaces = vectorOfSurfaces.size();

  //----------------------------------------------------------------------------
  // perform establishing of correspondences
  itk::TimeProbe tp;
  tp.Start();

  StringVector listOfOutputFiles;

  for (size_t stage = 0; stage < options.GetNumberOfStages(); ++stage) {
    std::cout << "establish correspondence stage (" << stage + 1 << " / " << options.GetNumberOfStages() << ")" << std::endl;

    // build GP model for the reference surface
    auto model = buildGPModel(reference, options.GetGPModelParameters()[0], options.GetGPModelScale(), options.GetGPModelNumberOfComponents()[0]);

    // initialize reference by zero
    for (auto iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
      iter.Value().Fill(0);
    }

    // preform GP model to surface registration
    for (size_t count = 0; count < numberOfSurfaces; ++count) {
      std::cout << "surface (" << count+1 << " / " << numberOfSurfaces << ")" << std::endl;

      auto surface = vectorOfSurfaces[count];
      auto output = shapeModelToSurfaceRegistration(surface, model, options);

      // add output to reference
      for (auto iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
        auto point = output->GetPoint(iter.Index());
        for (size_t d = 0; d < MeshType::PointDimension; ++d) {
          iter.Value()[d] += point[d] / numberOfSurfaces;
        }
      }

      // compute metrics and write output surface to file
      if (stage + 1 == options.GetNumberOfStages()) {
        // compute metrics
        typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilter;
        auto levelset = MeshToLevelSetImageFilter::New();
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
        auto metrics = PointSetToImageMetricsType::New();
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
    auto fileName = addFileNameSuffix(options.GetReferenceFileName(), "-" + std::to_string(stage + 1));
    if (!writeMesh<MeshType>(reference, fileName)) {
      EXIT_FAILURE;
    }
  }

  // write list of output files
  try {
    writeListToFile(options.GetOutputList(), listOfOutputFiles);
  }
  catch (std::ofstream::failure & e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  tp.Stop();
  std::cout << "elapsed time, sec " << tp.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, ShapeModelType::Pointer initialModel, const ssm::CorrespondenceOptions & options)
{
  // compute level set image
  typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilter;
  auto levelset = MeshToLevelSetImageFilter::New();
  levelset->SetInput(surface);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  // shape model to image registration
  auto model = initialModel;
  size_t numberOfStages = options.GetGPModelParameters().size();

  for (size_t stage = 0; stage < numberOfStages; ++stage) {
    std::cout << "registration stage (" << stage+1 << " / " << numberOfStages << ")" << std::endl;

    // if stage > 0 update GP model
    if (stage > 0) {
      model = buildGPModel(surface, options.GetGPModelParameters()[stage], options.GetGPModelScale(), options.GetGPModelNumberOfComponents()[stage]);
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
    auto initializer = TransformInitializerType::New();
    initializer->SetTransformType(options.GetTransform());
    initializer->SetCenter(center);
    initializer->SetTranslation(translation);
    try {
      initializer->Initialize();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }
    initializer->PrintReport();

    // perform registration
    typedef ssm::ShapeModelToImageRegistrationMethod<ShapeModelType, FloatImageType, MeshType> ShapeModelRegistrationMethod;
    auto shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetComputeLevelSetImage(false);
    shapeModelToSurfaceRegistration->SetNumberOfIterations(options.GetNumberOfIterations());
    shapeModelToSurfaceRegistration->SetSpatialTransform(initializer->GetTransform());
    shapeModelToSurfaceRegistration->SetSpatialScales(initializer->GetScales());
    shapeModelToSurfaceRegistration->GetMetric()->SetRegularizationParameter(options.GetGPModelRegularization()[stage]);
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
ShapeModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions)
{
  itk::TimeProbe tp;
  tp.Start();

  // create kernel
  typedef std::shared_ptr<const statismo::ScalarValuedKernel<PointType>> MatrixPointerType;
  MatrixPointerType kernel(new GaussianKernel<PointType>(parameters));

  typedef std::shared_ptr<statismo::MatrixValuedKernel<PointType>> KernelPointerType;
  KernelPointerType unscaledKernel(new statismo::UncorrelatedMatrixValuedKernel<PointType>(kernel.get(), Dimension));
  KernelPointerType modelBuildingKernel(new statismo::ScaledKernel<PointType>(unscaledKernel.get(), scale));

  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  auto representer = RepresenterType::New();
  representer->SetReference(surface);

  // build model
  typedef itk::LowRankGPModelBuilder<MeshType> ModelBuilderType;
  auto modelBuilder = ModelBuilderType::New();
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

  tp.Stop();
  std::cout << "Gaussian process model" << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << "parameters           " << "(" << parameters << ", " << scale << ")" << std::endl;
  std::cout << "elapsed time, sec    " << tp.GetTotal() << std::endl;
  std::cout << std::endl;

  return model;
}
