#include <boost/program_options.hpp>

#include <itkImageMomentsCalculator.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkStandardMeshRepresenter.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmShapeModelToImageRegistrationMethod.h"
#include "ssmMeshToLevelSetImageFilter.h"
#include "ssmInitializeSpatialTransform.h"

struct ProgramOptions
{
  bool help;
  std::string modelFile;
  std::string surfaceFile;
  std::string outputFile;
  std::string levelsetFile;
  std::string reportFile;
  std::string transformFile;
  size_t transform = 2;
  size_t iterations = 500;
  size_t degree = 2;
  double regularization = 0.1;
};

namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);

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

  //----------------------------------------------------------------------------
  // read surface
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, options.surfaceFile)) {
    return EXIT_FAILURE;
  }
  std::cout << "surface " << surface << std::endl;
  std::cout << "number of cells   " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points  " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read statistical shape model
  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  RepresenterType::Pointer representer = RepresenterType::New();

  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model = StatisticalModelType::New();
  try {
    model->Load(representer, options.modelFile.c_str());
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "model " << options.modelFile << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // compute level set image
  double margin = 0.10;
  BinaryImageType::SpacingType spacing(1);

  MeshType::BoundingBoxType::ConstPointer boundingBox = model->DrawMean()->GetBoundingBox();
  BinaryImageType::PointType origin = boundingBox->GetMinimum();
  BinaryImageType::SizeType size;
  for (size_t n = 0; n < Dimension; ++n) {
    size[n] = (boundingBox->GetMaximum()[n] - boundingBox->GetMinimum()[n]) / spacing[n];
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
    return EXIT_FAILURE;
  }

  // initialize spatial transform
  typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
  ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
  movingCalculator->SetImage(shapeToImage->GetOutput());
  movingCalculator->Compute();

  ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
  fixedCalculator->SetImage(shapeToImage->GetOutput());
  fixedCalculator->Compute();

  typedef ImageCalculatorType::VectorType VectorType;
  VectorType center = movingCalculator->GetCenterOfGravity();
  VectorType translation = fixedCalculator->GetCenterOfGravity() - movingCalculator->GetCenterOfGravity();

  typedef ssm::InitializeSpatialTransform<double> TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransformType(options.transform);
  initializer->SetCenter(center);
  initializer->SetTranslation(translation);
  try {
    initializer->Initialize();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializer->PrintReport();

  //----------------------------------------------------------------------------
  // perform shape model to image registration
  typedef ssm::ShapeModelToImageRegistrationMethod<StatisticalModelType, BinaryImageType, MeshType> RegistrationMethodType;
  auto shapeModelToSurfaceRegistration = RegistrationMethodType::New();
  shapeModelToSurfaceRegistration->SetShapeModel(model);
  shapeModelToSurfaceRegistration->SetImage(shapeToImage->GetOutput());
  shapeModelToSurfaceRegistration->SetNumberOfIterations(options.iterations);
  shapeModelToSurfaceRegistration->SetSpatialTransform(initializer->GetTransform());
  shapeModelToSurfaceRegistration->SetSpatialScales(initializer->GetScales());
  shapeModelToSurfaceRegistration->GetMetric()->SetRegularizationParameter(options.regularization);
  shapeModelToSurfaceRegistration->GetMetric()->SetDegree(options.degree);
  try {
    shapeModelToSurfaceRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  shapeModelToSurfaceRegistration->PrintReport();

  typedef std::pair<std::string, std::string> PairType;
  std::vector<PairType> info;
  info.push_back(PairType("Metric", std::to_string(shapeModelToSurfaceRegistration->GetOptimizer()->GetValue())));

  //----------------------------------------------------------------------------
  // compute metrics
  typedef itk::PointSet<MeshType::PointType, MeshType::PointDimension> PointSetType;
  typedef ssm::PointSetToImageMetrics<PointSetType, RegistrationMethodType::LevelSetImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetPointSetAsMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput());
  metrics->SetImage(shapeModelToSurfaceRegistration->GetLevelSetImage());
  metrics->SetInfo(info);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (options.reportFile != "") {
    metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(options.surfaceFile));
  }

  // write surface
  if (!writeMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput(), options.outputFile)) {
    return EXIT_FAILURE;
  }

  // write levelset image
  if (options.levelsetFile != "" ) {
    if (!writeImage(shapeModelToSurfaceRegistration->GetLevelSetImage(), options.reportFile)) {
      return EXIT_FAILURE;
    }
  }

  // write transform
  if (options.transformFile != "") {
    std::cout << "write transform to the file: " << options.transformFile << std::endl;
    writeTransform(initializer->GetTransform(), options.transformFile);
  }

  return EXIT_SUCCESS;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("model,m", po::value<std::string>(&options.modelFile), "The path to the input shape model file.")
    ("surface,s", po::value<std::string>(&options.surfaceFile), "The path to the input surface file.")
    ("output,o", po::value<std::string>(&options.outputFile), "The path for the output surface file.")
    ;

  po::options_description input("Optional input options");
  input.add_options()
    ("transform", po::value<size_t>(&options.transform)->default_value(options.transform), "The type of the used spatial transform.")
    ("iterations", po::value<size_t>(&options.iterations)->default_value(options.iterations), "The number of iterations.")
    ("degree", po::value<size_t>(&options.degree)->default_value(options.degree), "The degree of residuals to compute shape model to image metric.")
    ("regularization", po::value<double>(&options.regularization)->default_value(options.regularization), "The regularization factor.")
    ;

  po::options_description output("Optional output options");
  output.add_options()
    ("output-levelset", po::value<std::string>(&options.levelsetFile), "The path for the output level-set image.")
    ("output-transform", po::value<std::string>(&options.transformFile), "The path for the output transform file.")
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
