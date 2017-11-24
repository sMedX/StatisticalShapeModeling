#include <fstream>

#include <boost/program_options.hpp>

#include <itkImageMomentsCalculator.h>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>
#include <utils/statismo-build-models-utils.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmShapeModelToImageRegistrationMethod.h"
#include "ssmSurfaceToLevelSetImageFilter.h"

struct ProgramOptions
{
  bool help;
  std::string listFile;
  std::string surfaceFile;
  std::string reportFile;
  std::string referenceFile;
  std::string outputReferenceFile;

  std::vector<size_t> components;
  double scale = 100;
  std::vector<double> parameters;
  std::vector<double> regularization;
  size_t stages = 1;
  size_t transform = 1;
  size_t iterations = 500;
  size_t degree=2;
  double margin = 0.10;
  double spacing = 1.0;
};

typedef boost::filesystem::path fp;
namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer model, ProgramOptions & options);

int main(int argc, char** argv)
{
  ProgramOptions options;
  po::options_description description = initializeProgramOptions(options);
  po::variables_map vm;
  try {
    po::parsed_options parsedOptions = po::command_line_parser(argc, argv).options(description).run();
    po::store(parsedOptions, vm);
    po::notify(vm);
  }
  catch (po::error& e) {
    cerr << "An exception occurred while parsing the command line:" << endl;
    cerr << e.what() << endl << endl;
    cout << description << endl;
    return EXIT_FAILURE;
  }
  if (options.help == true) {
    cout << description << endl;
    return EXIT_SUCCESS;
  }

  if (options.parameters.size() == 0 || options.components.size() == 0 || options.regularization.size() == 0 ) {
    std::cout << "parameters, components and regularization factors mus be specified" << std::endl;
    return EXIT_FAILURE;
  }

  for (size_t n = options.components.size(); n < options.parameters.size(); ++n) {
    options.components.push_back(options.components.back());
  }

  for (size_t n = options.regularization.size(); n < options.parameters.size(); ++n) {
    options.regularization.push_back(options.regularization.back());
  }

  std::cout << "            parameters ";
  for (auto value : options.parameters) {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  std::cout << "            components ";
  for (auto value : options.components) {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  std::cout << "        regularization ";
  for (auto value : options.regularization) {
    std::cout << value << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read the reference surface
  MeshType::Pointer reference = MeshType::New();
  if (!readMesh<MeshType>(reference, options.referenceFile)) {
    return EXIT_FAILURE;
  }
  std::cout << "reference surface " << options.referenceFile << std::endl;
  std::cout << "number of cells " << reference->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << reference->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read list of files
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(options.listFile);
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

  for (size_t stage = 0; stage < options.stages; ++stage) {
    std::cout << "establish correspondence stage (" << stage + 1 << " / " << options.stages << ")" << std::endl;

    // build GP model for the reference surface
    StatisticalModelType::Pointer model = buildGPModel(reference, options.parameters[0], options.scale, options.components[0]);

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
      if (stage + 1 == options.stages) {
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
        if (options.reportFile != "") {
          std::cout << "print report to file " << options.reportFile << std::endl;
          metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(vectorOfFiles[count]));
        }

        // write output surface to file
        typedef boost::filesystem::path fp;
        fp path = fp(options.surfaceFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(options.surfaceFile).filename().string());
        std::string fileName = path.string();
        std::cout << "write surface to file " << fileName << std::endl;
        if (!writeMesh<MeshType>(output, fileName)) {
          EXIT_FAILURE;
        }
      }
    }

    // write the reference surface to the file
    if (options.outputReferenceFile != "") {
      std::string fileName = addFileNameSuffix(options.outputReferenceFile, "-" + std::to_string(stage + 1));
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
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer initialModel, ProgramOptions & options)
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
      model = buildGPModel(surface, options.parameters[stage], options.scale, options.components[stage]);
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

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("list,l", po::value<std::string>(&options.listFile), "The path to the file with list of surfaces to establish correspondence.")
    ("reference,r", po::value<std::string>(&options.referenceFile), "The path to the input reference surface.")
    ("surface,s", po::value<std::string>(&options.surfaceFile), "The path for the output surfaces.")
    ;

  po::options_description input("Optional input options");
  input.add_options()
    ("components", po::value<std::vector<size_t>>(&options.components)->multitoken(), "The number of components to build GP shape model.")
    ("scale", po::value<double>(&options.scale)->default_value(options.scale), "The scaling factor to scale the kernel.")
    ("parameters", po::value<std::vector<double>>(&options.parameters)->multitoken(), "The parameters to build GP shape model.")
    ("regularization", po::value<std::vector<double>>(&options.regularization)->multitoken(), "The regularization factor.")
    ("stages", po::value<size_t>(&options.stages)->default_value(options.stages), "The number of stages to establish correspondence.")
    ("transform", po::value<size_t>(&options.transform)->default_value(options.transform), "The type of the used spatial transform.")
    ("iterations", po::value<size_t>(&options.iterations)->default_value(options.iterations), "The number of iterations.")
    ("degree", po::value<size_t>(&options.degree)->default_value(options.degree), "The degree of residuals to compute shape model to image metric.")
    ;

  po::options_description output("Optional output options");
  output.add_options()
    ("output-reference", po::value<std::string>(&options.outputReferenceFile), "The path for the output updated reference surfaces.")
    ;

  po::options_description report("Optional report options");
  report.add_options()
    ("report,r", po::value<std::string>(&options.reportFile), "The path for the file to print report.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(input).add(output).add(report).add(help);

  return description;
}
