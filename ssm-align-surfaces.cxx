#include <boost/program_options.hpp>

#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmSurfaceToLevelSetImageFilter.h"
#include "ssm/ssmSurfaceToImageRegistrationMethod.h"
#include "ssm/ssmMeshPropertiesCalculator.h"

typedef std::vector<MeshType::Pointer> MeshVectorType;
typedef itk::Transform<double, MeshType::PointDimension> TransformType;
typedef std::vector<TransformType::ConstPointer> TransformVectorType;
FloatImageType::Pointer computeLevelSetImage(FloatImageType::Pointer levelSetImage, MeshVectorType & vectorOfSurfaces, TransformVectorType & vectorOfTransform, size_t typeOfransform);

struct ProgramOptions
{
  bool help;
  std::string listFile;
  std::string surfaceFile;
  std::string levelsetFile;
  std::string reportFile;
  size_t transform = 1;
  size_t stages = 1;
  size_t iterations = 500;
};

typedef boost::filesystem::path fp;
namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);

int main(int argc, char** argv) {
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

  // read list of files
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(options.listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "Could not read the data-list:" << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  MeshVectorType vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (const auto & fileName : listOfFiles) {
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    std::cout << "surface " << fileName << std::endl;
    std::cout << "number of cells  " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;
  }

  std::cout << "number of surfaces " << vectorOfSurfaces.size() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // compute initial level set image

  // compute size of the initial level set image
  typedef itk::Vector<double, Dimension> VectorType;
  std::vector<VectorType> vectorOfCenters;
  itk::Matrix<double, Dimension, 2> maximalBoundingBox;
  for (size_t i = 0; i < Dimension; ++i) {
    maximalBoundingBox[i][0] = std::numeric_limits<double>::max();
    maximalBoundingBox[i][1] = std::numeric_limits<double>::min();
  }

  VectorType centerOfMaximalBoundingBox;
  centerOfMaximalBoundingBox.Fill(0);

  for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
    MeshType::BoundingBoxType::ConstPointer boundingbox = vectorOfSurfaces[count]->GetBoundingBox();

    // compute mask for surface
    BinaryImageType::SpacingType spacing(1);
    BinaryImageType::SizeType size;

    for (unsigned i = 0; i < Dimension; ++i) {
      spacing[i] = 1;
      size[i] = (boundingbox->GetMaximum()[i] - boundingbox->GetMinimum()[i]) / spacing[i];
    }

    typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> TriangleMeshToBinaryImageFilterType;
    TriangleMeshToBinaryImageFilterType::Pointer surfaceToMask = TriangleMeshToBinaryImageFilterType::New();
    surfaceToMask->SetInput(vectorOfSurfaces[count]);
    surfaceToMask->SetInsideValue(1);
    surfaceToMask->SetSize(size);
    surfaceToMask->SetSpacing(spacing);
    surfaceToMask->SetOrigin(vectorOfSurfaces[count]->GetBoundingBox()->GetMinimum());
    try {
      surfaceToMask->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    // compute center
    typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
    ImageCalculatorType::Pointer calculator = ImageCalculatorType::New();
    calculator->SetImage(surfaceToMask->GetOutput());
    calculator->Compute();
    VectorType center = calculator->GetCenterOfGravity();

    vectorOfCenters.push_back(center);

    for (unsigned int i = 0; i < Dimension; ++i) {
      maximalBoundingBox[i][0] = std::min(maximalBoundingBox[i][0], boundingbox->GetMinimum()[i] - center[i]);
      maximalBoundingBox[i][1] = std::max(maximalBoundingBox[i][1], boundingbox->GetMaximum()[i] - center[i]);
    }
  }

  std::cout << "maximal bounding box, mm" << std::endl;
  std::cout << maximalBoundingBox << std::endl;
  std::cout << "center " << centerOfMaximalBoundingBox << std::endl;

  // allocate zero image
  FloatImageType::SpacingType spacing(1);
  FloatImageType::PointType origin;
  itk::Size<Dimension> size;

  for (size_t i = 0; i < Dimension; ++i) {
    double offset = 0.25;
    double distance = maximalBoundingBox[i][1] - maximalBoundingBox[i][0];

    size[i] = (1 + 2*offset) * distance / spacing[i];
    origin[i] = maximalBoundingBox[i][0] - offset*distance;
  }

  FloatImageType::RegionType region;
  region.SetSize(size);

  FloatImageType::Pointer levelSetImage = FloatImageType::New();
  levelSetImage->SetRegions(region);
  levelSetImage->Allocate();
  levelSetImage->FillBuffer(0);
  levelSetImage->SetSpacing(spacing);
  levelSetImage->SetOrigin(origin);

  std::cout << levelSetImage->GetLargestPossibleRegion() << std::endl;
  std::cout << "spacing " << levelSetImage->GetSpacing() << std::endl;
  std::cout << " origin " << levelSetImage->GetOrigin() << std::endl;

  // compute initial level-set image
  for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
    // transform the n-th surface
    typedef itk::TranslationTransform<double, Dimension> TransformType;
    TransformType::Pointer transform = TransformType::New();
    transform->SetOffset(centerOfMaximalBoundingBox - vectorOfCenters[count]);

    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformFilterType;
    TransformFilterType::Pointer transformSurface = TransformFilterType::New();
    transformSurface->SetInput(vectorOfSurfaces[count]);
    transformSurface->SetTransform(transform);
    try {
      transformSurface->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      return EXIT_FAILURE;
    }
    vectorOfSurfaces[count] = transformSurface->GetOutput();

    // compute level-set image for the n-th surface
    typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
    SurfaceToLevelSetImageFilterType::Pointer surfaceToLevelSetImage = SurfaceToLevelSetImageFilterType::New();
    surfaceToLevelSetImage->SetInput(vectorOfSurfaces[count]);
    surfaceToLevelSetImage->SetOrigin(levelSetImage->GetOrigin());
    surfaceToLevelSetImage->SetSpacing(levelSetImage->GetSpacing());
    surfaceToLevelSetImage->SetSize(levelSetImage->GetLargestPossibleRegion().GetSize());
    try {
      surfaceToLevelSetImage->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    // add the n-th image to the level-set image
    typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
    FilterType::Pointer multiply = FilterType::New();
    multiply->SetInput(surfaceToLevelSetImage->GetOutput());
    multiply->SetConstant(1 / (double)vectorOfSurfaces.size());

    typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
    AddImageFilterType::Pointer add = AddImageFilterType::New();
    add->SetInput1(levelSetImage);
    add->SetInput2(multiply->GetOutput());
    try {
      add->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    levelSetImage = add->GetOutput();
  }

  //----------------------------------------------------------------------------
  // perform alignment of the surfaces
  TransformVectorType vectorOfTransforms;
  vectorOfTransforms.resize(vectorOfSurfaces.size());

  for (size_t stage = 0; stage < options.stages; ++stage) {
    for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
      std::cout << "stage " << stage + 1 << "/" << options.stages << ", surface " << count + 1 << "/" << vectorOfSurfaces.size() << ", " << vectorOfFiles[count] << std::endl;

      // perform surface to image registration
      typedef ssm::SurfaceToImageRegistrationMethod<MeshType> SurfaceToImageRegistrationMethodType;
      SurfaceToImageRegistrationMethodType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationMethodType::New();
      surfaceToImageRegistration->SetInput(vectorOfSurfaces[count]);
      surfaceToImageRegistration->SetNumberOfIterations(options.iterations);
      surfaceToImageRegistration->SetTypeOfTransform(options.transform);
      surfaceToImageRegistration->SetLevelsetImage(levelSetImage);
      try {
        surfaceToImageRegistration->Update();
      }
      catch (itk::ExceptionObject& excep) {
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
      }
      vectorOfTransforms[count] = surfaceToImageRegistration->GetTransform();
      surfaceToImageRegistration->PrintReport(std::cout);
    }

    // update reference level-set image
    levelSetImage = computeLevelSetImage(levelSetImage, vectorOfSurfaces, vectorOfTransforms, options.transform);
  }

  //----------------------------------------------------------------------------
  // write reference level set image
  if (options.levelsetFile != "") {
    std::cout << "output the level-set image " << options.levelsetFile << std::endl;
    std::cout << "   size " << levelSetImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "spacing " << levelSetImage->GetSpacing() << std::endl;
    std::cout << " origin " << levelSetImage->GetOrigin() << std::endl;
    std::cout << std::endl;

    typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
    FilterType::Pointer multiply = FilterType::New();
    multiply->SetInput(levelSetImage);
    multiply->SetConstant(-1);

    if (!writeImage<FloatImageType>(multiply->GetOutput(), options.levelsetFile)) {
      return EXIT_FAILURE;
    }
  }

  // write alignment surfaces
  for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {

    // define full file name for output surface
    fp path = fp(options.surfaceFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(options.surfaceFile).filename().string());
    std::string outputSurfaceFile = path.string();

    std::cout << "output file " << outputSurfaceFile << std::endl;
    if (!writeMesh<MeshType>(vectorOfSurfaces[count], outputSurfaceFile)) {
      EXIT_FAILURE;
    }

    // compute metrics
    typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
    PointSetType::Pointer pointSet = PointSetType::New();
    pointSet->SetPoints(vectorOfSurfaces[count]->GetPoints());

    typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
    PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
    metrics->SetFixedPointSet(pointSet);
    metrics->SetMovingImage(levelSetImage);
    metrics->Compute();
    metrics->PrintReport(std::cout);

    // print report to *.csv file
    if (options.reportFile!="") {
      std::cout << "print report to the file: " << options.reportFile << std::endl;
      metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(outputSurfaceFile));
    }
  }

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
FloatImageType::Pointer computeLevelSetImage(FloatImageType::Pointer levelSetImage, MeshVectorType & vectorOfSurfaces, TransformVectorType & vectorOfTransform, size_t typeOfransform)
{
  double scale = 1;
  size_t index = 6;

  if (typeOfransform > 1) {
    double initialRadius = 0;
    double scaledRadius = 0;

    for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
      MeshType::Pointer surface = vectorOfSurfaces[count];
      TransformType::ConstPointer transform = vectorOfTransform[count];

      typedef ssm::MeshPropertiesCalculator<MeshType> MeshMomentsCalculatorType;
      MeshMomentsCalculatorType::Pointer calculator = MeshMomentsCalculatorType::New();
      calculator->SetMesh(surface);
      calculator->Compute();

      double radius = calculator->GetRadius();
      double scale = transform->GetParameters()[index];

      initialRadius += radius;
      scaledRadius += scale * radius;
    }

    scale = initialRadius / scaledRadius;
  }

  // fill level-set image
  levelSetImage->FillBuffer(0);

  for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
    // modify transform
    TransformType::Pointer transform = const_cast<TransformType*>(vectorOfTransform[count].GetPointer());

    if (typeOfransform > 1) {
      TransformType::ParametersType parameters = transform->GetParameters();
      parameters[index] *= scale;
      transform->SetParameters(parameters);
    }

    // transform surface
    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformFilterType;
    TransformFilterType::Pointer transformMeshFilter = TransformFilterType::New();
    transformMeshFilter->SetInput(vectorOfSurfaces[count]);
    transformMeshFilter->SetTransform(transform);
    try {
      transformMeshFilter->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      throw;
    }
    vectorOfSurfaces[count] = transformMeshFilter->GetOutput();

    // compute level-set image
    typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
    SurfaceToLevelSetImageFilterType::Pointer surfaceToLevelSetImage = SurfaceToLevelSetImageFilterType::New();
    surfaceToLevelSetImage->SetInput(vectorOfSurfaces[count]);
    surfaceToLevelSetImage->SetOrigin(levelSetImage->GetOrigin());
    surfaceToLevelSetImage->SetSpacing(levelSetImage->GetSpacing());
    surfaceToLevelSetImage->SetSize(levelSetImage->GetLargestPossibleRegion().GetSize());
    try {
      surfaceToLevelSetImage->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }

    // add current level set image to update reference image
    typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
    FilterType::Pointer multiply = FilterType::New();
    multiply->SetInput(surfaceToLevelSetImage->GetOutput());
    multiply->SetConstant(1 / (double)vectorOfSurfaces.size());

    typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
    AddImageFilterType::Pointer add = AddImageFilterType::New();
    add->SetInput1(levelSetImage);
    add->SetInput2(multiply->GetOutput());
    try {
      add->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }

    levelSetImage = add->GetOutput();
  }

  return levelSetImage;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("list,l", po::value<std::string>(&options.listFile), "The path to the file with list of surfaces to align.")
    ("output,o", po::value<std::string>(&options.surfaceFile), "The path to the output surface file.")
    ;

  po::options_description input("Optional input options");
  input.add_options()
    ("stages", po::value<size_t>(&options.stages), "The number of stages to align input surfaces.")
    ("transform", po::value<size_t>(&options.transform), "The type of the used spatial transform.")
    ("iterations,t", po::value<size_t>(&options.iterations), "The number of iterations.")
    ;

  po::options_description output("Optional output options");
  output.add_options()
    ("levelset,l", po::value<std::string>(&options.levelsetFile), "The path for output level-set image.")
    ;

  po::options_description report("Optional report options");
  report.add_options()
    ("report,r", po::value<std::string>(&options.reportFile), "The path to the file to print report.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(input).add(output).add(report).add(help);

  return description;
}
