#include <boost/program_options.hpp>
#include <itkImageMomentsCalculator.h>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmShapeModelToImageRegistrationMethod.h"
#include "ssm/ssmSurfaceToLevelSetImageFilter.h"

struct RegistrationOptions
{
  size_t transform = 1;
  std::vector<double> parameters;
  std::vector<double> regularization;
  double scale = 100;
  unsigned int iterations = 500;
  unsigned int degree;
  double margin = 0.10;
  double spacing = 1.0;
};

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer model, RegistrationOptions & options);

int main(int argc, char** argv)
{
  RegistrationOptions options;

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();
  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  std::string referenceFile;
  parser->GetCommandLineArgument("-reference", referenceFile);

  parser->GetCommandLineArgument("-transform", options.transform);
  parser->GetCommandLineArgument("-iterations", options.iterations);
  parser->GetCommandLineArgument("-regularization", options.regularization);
  parser->GetCommandLineArgument("-parameters", options.parameters);
  parser->GetCommandLineArgument("-scale", options.scale);
  parser->GetCommandLineArgument("-degree", options.degree);

  unsigned int numberOfStages = 1;
  parser->GetCommandLineArgument("-stages", numberOfStages);

  unsigned int numberOfComponents = 100;
  parser->GetCommandLineArgument("-components", numberOfComponents);

  std::cout << std::endl;
  std::cout << " shape model to image registration" << std::endl;
  std::cout << "   output surface file " << outputFile << std::endl;
  std::cout << "      number of stages " << numberOfStages << std::endl;
  std::cout << "  number of iterations " << options.iterations << std::endl;
  std::cout << "                degree " << options.degree << std::endl;

  for (int n = options.regularization.size(); n < options.parameters.size(); ++n) {
    options.regularization.push_back(options.regularization.back());
  }

  std::cout << "        regularization ";
  for (auto value : options.regularization) {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  std::cout << "            parameters ";
  for (auto value : options.parameters) {
    std::cout << value << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read the reference surface
  MeshType::Pointer reference = MeshType::New();
  if (!readMesh<MeshType>(reference, referenceFile)) {
    return EXIT_FAILURE;
  }
  std::cout << "reference " << referenceFile << std::endl;
  std::cout << "number of cells " << reference->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << reference->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read list of files
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "could not read the data-list: " << e.what() << endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (const auto & fileName : listOfFiles) {
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    std::cout << "surface " << fileName << std::endl;
    std::cout << "number of cells   " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points  " << surface->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;
  }

  size_t numberOfSurfaces = vectorOfSurfaces.size();

  //----------------------------------------------------------------------------
  // perform establishing of correspondences
  itk::TimeProbe clock;
  clock.Start();

  for (size_t stage = 0; stage < numberOfStages; ++stage) {
    std::cout << "establish correspondence stage (" << stage+1 << " / " << numberOfStages << ")" << std::endl;

    // build GP model for the reference surface
    StatisticalModelType::Pointer model = BuildGPModel(reference, options.parameters[0], options.scale, numberOfComponents);

    // initialize reference to zero
    typedef itk::VectorContainer<long unsigned int,MeshType::PointType> VContainer;
    
    for (VContainer::Iterator iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
      iter.Value().Fill(0);
    }

    // preform GP model to surface registration
    for (size_t count = 0; count < numberOfSurfaces; ++count) {
      std::cout << "surface (" << count+1 << " / " << numberOfSurfaces << ")" << std::endl;

      MeshType::Pointer surface = vectorOfSurfaces[count];
      MeshType::Pointer output = shapeModelToSurfaceRegistration(surface, model, options);

      // add output to reference
      for (VContainer::Iterator iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
        MeshType::PointType point = output->GetPoint(iter.Index());
        for (size_t d = 0; d < MeshType::PointDimension; ++d) {
          iter.Value()[d] += point[d] / numberOfSurfaces;
        }
      }

      // compute metrics and write output surface to file
      if (stage + 1 == numberOfStages) {
        // compute metrics
        typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
        SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
        levelset->SetMargin(options.margin);
        levelset->SetSpacing(options.spacing);
        levelset->SetInput(surface);
        try {
          levelset->Update();
        }
        catch (itk::ExceptionObject& excep) {
          std::cerr << excep << std::endl;
          return EXIT_FAILURE;
        }

        typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
        PointSetType::Pointer points = PointSetType::New();
        points->SetPoints(output->GetPoints());

        typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
        PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
        metrics->SetFixedPointSet(points);
        metrics->SetMovingImage(levelset->GetOutput());
        metrics->Compute();
        metrics->PrintReport(std::cout);

        // print report to *.csv file
        if (parser->ArgumentExists("-report")) {
          std::string fileName;
          parser->GetCommandLineArgument("-report", fileName);
          std::cout << "print report to file " << fileName << std::endl;
          metrics->PrintReportToFile(fileName, getBaseNameFromPath(vectorOfFiles[count]));
        }

        // write output surface to file
        typedef boost::filesystem::path fp;
        fp path = fp(outputFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(outputFile).filename().string());
        std::string fileName = path.string();
        std::cout << "write surface to file " << fileName << std::endl;
        if (!writeMesh<MeshType>(output, fileName)) {
          EXIT_FAILURE;
        }
      }
    }

    // write reference to file
    if (parser->ArgumentExists("-output-reference")) {
      std::string fileName;
      parser->GetCommandLineArgument("-output-reference", fileName);
      fileName = addFileNameSuffix(fileName, "-" + std::to_string(stage + 1));
      if (!writeMesh<MeshType>(reference, fileName)) {
        EXIT_FAILURE;
      }
    }
  }

  clock.Stop();
  std::cout << "elapsed time, sec " << clock.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer initialModel, RegistrationOptions & options)
{
  // compute level set image
  typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
  SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
  levelset->SetMargin(options.margin);
  levelset->SetSpacing(options.spacing);
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
  size_t numberOfStages = options.parameters.size();

  for (size_t stage = 0; stage < numberOfStages; ++stage) {
    std::cout << "registration stage (" << stage+1 << " / " << numberOfStages << ")" << std::endl;

    // if stage > 0 update GP model
    if (stage > 0) {
      model = BuildGPModel(surface, options.parameters[stage], options.scale, model->GetNumberOfPrincipalComponents());
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
    ShapeToBinaryImageFilterType::Pointer shapeToImage = ShapeToBinaryImageFilterType::New();
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
    ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
    movingCalculator->SetImage(shapeToImage->GetOutput());
    movingCalculator->Compute();

    ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
    fixedCalculator->SetImage(levelset->GetMask());
    fixedCalculator->Compute();

    typedef ImageCalculatorType::VectorType VectorType;
    VectorType center = movingCalculator->GetCenterOfGravity();
    VectorType translation = fixedCalculator->GetCenterOfGravity() - movingCalculator->GetCenterOfGravity();

    typedef ssm::TransformInitializer<double> TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTypeOfTransform(options.transform);
    initializer->SetCenter(center);
    initializer->SetTranslation(translation);
    try {
      initializer->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      throw;
    }
    initializer->PrintReport(std::cout);

    // perform registration
    typedef ssm::ShapeModelToImageRegistrationMethod<StatisticalModelType, MeshType> ShapeModelRegistrationMethod;
    ShapeModelRegistrationMethod::Pointer shapeModelToSurfaceRegistration;
    shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetNumberOfIterations(options.iterations);
    shapeModelToSurfaceRegistration->SetRegularizationParameter(options.regularization[stage]);
    shapeModelToSurfaceRegistration->SetDegree(options.degree);
    shapeModelToSurfaceRegistration->SetTransformInitializer(initializer);
    try {
      shapeModelToSurfaceRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }
    shapeModelToSurfaceRegistration->PrintReport(std::cout);

    surface = const_cast<MeshType*> (shapeModelToSurfaceRegistration->GetOutput());
  }

  return surface;
}

//==============================================================================
// Build Gaussian process model
StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions)
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
