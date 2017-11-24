#include <fstream>

#include <itkImageMomentsCalculator.h>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmShapeModelToImageRegistrationMethod.h"
#include "ssmSurfaceToLevelSetImageFilter.h"
#include "ssmCorrespondenceOptions.h"

/*
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
*/

//typedef boost::filesystem::path fp;
//namespace po = boost::program_options;
//po::options_description initializeProgramOptions(ProgramOptions& poParameters);

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer buildGPModel(MeshType::Pointer surface, double parameters, double scale, size_t numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer model, ssm::CorrespondenceOptions & options);

template<typename T>
std::vector<T> to_array(const std::string& str)
{
  std::stringstream stream(str);

  std::string item;
  std::vector<T> result;

  while (std::getline(stream, item, ' ')) {
      try {
        result.push_back(std::stod(item));
      }
      catch (...) {
        std::cout << "error " << item << std::endl;
        return result;
      }
  }

  return result;
}

int main(int argc, char** argv)
{
  /*
  boost::property_tree::ptree pt;
  boost::property_tree::read_ini("config.ini", pt);
  ssm::printTree(pt, std::cout, 0);

  auto x = to_array<double>(pt.get<std::string>("CORRESPONDENCE.parameters"));
  std::cout << "size " << x.size() << std::endl;

  for (auto i : x) {
    std::cout << i << ' ';
  }
  std::cout << '\n';


  return 0;
  */
  ssm::CorrespondenceOptions options;
  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  // read options from config file
  if (!options.ReadConfigFile()) {
    return EXIT_FAILURE;
  }
  options.PrintConfig();
  /*
  auto parameters = options.GetParameters();
  auto components = options.GetNumberOfComponents();
  auto regularization = options.GetRegularization();

  if (parameters.size() == 0 || components.size() == 0 || regularization.size() == 0 ) {
    std::cout << "parameters, components and regularization factors must be specified" << std::endl;
    return EXIT_FAILURE;
  }

  for (size_t n = components.size(); n < parameters.size(); ++n) {
    components.push_back(components.back());
  }

  for (size_t n = regularization.size(); n < parameters.size(); ++n) {
    regularization.push_back(regularization.back());
  }

  std::cout << "    parameters ";
  for (auto value : parameters) {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  std::cout << "    components ";
  for (auto value : components) {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  std::cout << "regularization ";
  for (auto value : regularization) {
    std::cout << value << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;
  */
  //----------------------------------------------------------------------------
  // read the reference surface
  auto reference = MeshType::New();
  if (!readMesh<MeshType>(reference, options.GetReferenceFileName())) {
    return EXIT_FAILURE;
  }
  printMeshInfo<MeshType>(reference, options.GetReferenceFileName());

  //----------------------------------------------------------------------------
  // read list of files
  StringList listOfInputFiles;
  try {
    listOfInputFiles = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (const auto & fileName : listOfInputFiles) {
    auto surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    printMeshInfo<MeshType>(surface, fileName);
  }

  size_t numberOfSurfaces = vectorOfSurfaces.size();

  //----------------------------------------------------------------------------
  // perform establishing of correspondences
  itk::TimeProbe time;
  time.Start();

  StringList listOfOutputFiles;

  for (size_t stage = 0; stage < options.GetNumberOfStages(); ++stage) {
    std::cout << "correspondence stage (" << stage + 1 << " / " << options.GetNumberOfStages() << ")" << std::endl;

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
        typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
        auto levelset = SurfaceToLevelSetImageFilter::New();
        levelset->SetMargin(0.10);
        levelset->SetSpacing(1);
        levelset->SetInput(surface);
        try {
          levelset->Update();
        }
        catch (itk::ExceptionObject& excep) {
          std::cerr << excep << std::endl;
          return EXIT_FAILURE;
        }

        typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
        auto points = PointSetType::New();
        points->SetPoints(output->GetPoints());

        typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
        PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
        metrics->SetFixedPointSet(points);
        metrics->SetMovingImage(levelset->GetOutput());
        metrics->Compute();
        metrics->PrintReport(std::cout);

        // print report to *.csv file
        std::cout << "print report to the file: " << options.GetReportFileName() << std::endl;
        std::cout << std::endl;
        metrics->PrintReportToFile(options.GetReportFileName(), getBaseNameFromPath(vectorOfFiles[count]));

        // write output surface to file
        auto fileName = options.FormatOutput(vectorOfFiles[count]);
        listOfOutputFiles.push_back(fileName);

        std::cout << "output file " << fileName << std::endl;
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

  // write list of files
  try {
    writeListToFile(options.GetOutputList(), listOfOutputFiles);
  }
  catch (std::ofstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  time.Stop();
  std::cout << "elapsed time, sec " << time.GetTotal() << std::endl;

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer initialModel, ssm::CorrespondenceOptions & options)
{
  // compute level set image
  typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
  auto levelset = SurfaceToLevelSetImageFilter::New();
  levelset->SetMargin(0.1);
  levelset->SetSpacing(1);
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
    std::cout << "registration stage (" << stage + 1 << " / " << numberOfStages << ")" << std::endl;

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

    auto center = movingCalculator->GetCenterOfGravity();
    auto translation = fixedCalculator->GetCenterOfGravity() - movingCalculator->GetCenterOfGravity();

    typedef ssm::TransformInitializer<double> TransformInitializerType;
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    initializer->SetTypeOfTransform(options.GetTransform());
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
    auto shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetNumberOfIterations(options.GetNumberOfIterations());
    shapeModelToSurfaceRegistration->SetRegularizationParameter(options.GetRegularization()[stage]);
    shapeModelToSurfaceRegistration->SetDegree(2);
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
  itk::TimeProbe time;
  time.Start();

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

  time.Stop();
  std::cout << "Gaussian process model" << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << "parameters           " << "(" << parameters << ", " << scale << ")" << std::endl;
  std::cout << "elapsed time, sec    " << time.GetTotal() << std::endl;
  std::cout << std::endl;

  return model;
}
/*
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
*/