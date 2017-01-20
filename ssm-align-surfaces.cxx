#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmSurfaceToLevelSetImageFilter.h"
#include "ssm/ssmSurfaceToImageRegistrationMethod.h"

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-output", surfaceFile);

  unsigned int numberOfStages = 1;
  parser->GetCommandLineArgument("-stages", numberOfStages);

  unsigned int transform = 1;
  parser->GetCommandLineArgument("-transform", transform);

  unsigned int numberOfIterations = 500;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  std::cout << std::endl;
  std::cout << "parameters" << std::endl;
  std::cout << "    list of files " << listFile << std::endl;
  std::cout << "   output surface " << surfaceFile << std::endl;
  std::cout << " number of stages " << numberOfStages << std::endl;
  std::cout << "type of transform " << transform << std::endl;
  std::cout << "       iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  // read list of files
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "Could not read the data-list:" << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (StringList::const_iterator it = listOfFiles.begin(); it != listOfFiles.end(); ++it) {
    std::string fileName = it->c_str();
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    std::cout << "input surface polydata info " << fileName << std::endl;
    std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
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
  for (unsigned int i = 0; i < Dimension; ++i) {
    maximalBoundingBox[i][0] = std::numeric_limits<double>::max();
    maximalBoundingBox[i][1] = std::numeric_limits<double>::min();
  }

  VectorType centerOfMaximalBoundingBox;
  centerOfMaximalBoundingBox.Fill(0);

  for (int n = 0; n < vectorOfSurfaces.size(); ++n) {
    MeshType::BoundingBoxType::ConstPointer boundingbox = vectorOfSurfaces[n]->GetBoundingBox();

    // compute mask for surface
    BinaryImageType::SpacingType spacing(1);
    BinaryImageType::SizeType size;

    for (unsigned i = 0; i < Dimension; ++i) {
      spacing[i] = 1;
      size[i] = (boundingbox->GetMaximum()[i] - boundingbox->GetMinimum()[i]) / spacing[i];
    }

    typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> TriangleMeshToBinaryImageFilterType;
    TriangleMeshToBinaryImageFilterType::Pointer surfaceToMask = TriangleMeshToBinaryImageFilterType::New();
    surfaceToMask->SetInput(vectorOfSurfaces[n]);
    surfaceToMask->SetInsideValue(1);
    surfaceToMask->SetSize(size);
    surfaceToMask->SetSpacing(spacing);
    surfaceToMask->SetOrigin(vectorOfSurfaces[n]->GetBoundingBox()->GetMinimum());
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

  for (unsigned int i = 0; i < Dimension; ++i) {
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
  for (int n = 0; n < vectorOfSurfaces.size(); ++n) {
    // transform the n-th surface
    typedef itk::TranslationTransform<double, Dimension> TransformType;
    TransformType::Pointer transform = TransformType::New();
    transform->SetOffset(centerOfMaximalBoundingBox - vectorOfCenters[n]);

    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformFilterType;
    TransformFilterType::Pointer transformSurface = TransformFilterType::New();
    transformSurface->SetInput(vectorOfSurfaces[n]);
    transformSurface->SetTransform(transform);
    try {
      transformSurface->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cout << excep << std::endl;
      return EXIT_FAILURE;
    }
    vectorOfSurfaces[n] = transformSurface->GetOutput();

    // compute level-set image for the n-th surface
    typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
    SurfaceToLevelSetImageFilterType::Pointer surfaceToLevelSetImage = SurfaceToLevelSetImageFilterType::New();
    surfaceToLevelSetImage->SetInput(vectorOfSurfaces[n]);
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
  for (int stage = 0; stage < numberOfStages; ++stage) {
    typedef ssm::SurfaceToImageRegistrationMethod<MeshType> SurfaceToImageRegistrationMethodType;
    typedef SurfaceToImageRegistrationMethodType::EnumTransformType EnumTransformType;
    EnumTransformType typeOfTransform = static_cast<EnumTransformType>(transform);

    std::cout << "perform registration" << std::endl;
    std::cout << "stage " << stage + 1 << "/" << numberOfStages << std::endl;
    std::cout << "type of transform " << typeOfTransform << std::endl;
    std::cout << std::endl;

    // allocate image to update reference image
    FloatImageType::Pointer updateLevelSetImage = FloatImageType::New();
    updateLevelSetImage->SetRegions(levelSetImage->GetLargestPossibleRegion());
    updateLevelSetImage->Allocate();
    updateLevelSetImage->FillBuffer(0);
    updateLevelSetImage->CopyInformation(levelSetImage);

    for (int n = 0; n < vectorOfSurfaces.size(); ++n) {
      std::cout << "stage " << stage + 1 << "/" << numberOfStages << ", surface " << n + 1 << "/" << vectorOfSurfaces.size() << ", " << vectorOfFiles[n] << std::endl;

      // perform surface to image registration
      typedef ssm::SurfaceToImageRegistrationMethod<MeshType> SurfaceToImageRegistrationMethodType;
      SurfaceToImageRegistrationMethodType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationMethodType::New();
      surfaceToImageRegistration->SetInput(vectorOfSurfaces[n]);
      surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
      surfaceToImageRegistration->SetTypeOfTransform(typeOfTransform);
      surfaceToImageRegistration->SetLevelsetImage(levelSetImage);
      try {
        surfaceToImageRegistration->Update();
      }
      catch (itk::ExceptionObject& excep) {
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
      }
      vectorOfSurfaces[n] = surfaceToImageRegistration->GetOutput();
      surfaceToImageRegistration->PrintReport(std::cout);

      // compute level set image
      typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
      SurfaceToLevelSetImageFilterType::Pointer surfaceToLevelSetImage = SurfaceToLevelSetImageFilterType::New();
      surfaceToLevelSetImage->SetInput(vectorOfSurfaces[n]);
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

      // add current level set image to update reference image
      typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
      FilterType::Pointer multiply = FilterType::New();
      multiply->SetInput(surfaceToLevelSetImage->GetOutput());
      multiply->SetConstant(1 / (double) vectorOfSurfaces.size());

      typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
      AddImageFilterType::Pointer add = AddImageFilterType::New();
      add->SetInput1(updateLevelSetImage);
      add->SetInput2(multiply->GetOutput());
      try {
        add->Update();
      }
      catch (itk::ExceptionObject& excep) {
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
      }

      updateLevelSetImage = add->GetOutput();
    }

    //update reference image
    levelSetImage = updateLevelSetImage;
  }

  //----------------------------------------------------------------------------
  // write reference level set image
  if (parser->ArgumentExists("-levelset")) {
    std::string fileName;
    parser->GetCommandLineArgument("-levelset", fileName);

    std::cout << "output the level-set image " << fileName << std::endl;
    std::cout << "   size " << levelSetImage->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "spacing " << levelSetImage->GetSpacing() << std::endl;
    std::cout << " origin " << levelSetImage->GetOrigin() << std::endl;
    std::cout << std::endl;

    typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
    FilterType::Pointer multiply = FilterType::New();
    multiply->SetInput(levelSetImage);
    multiply->SetConstant(-1);

    if (!writeImage<FloatImageType>(multiply->GetOutput(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  // write alignment surfaces
  for (int count = 0; count < vectorOfSurfaces.size(); ++count) {

    // define full file name for output surface
    typedef boost::filesystem::path fp;
    fp path = fp(surfaceFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(surfaceFile).filename().string());
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
    if (parser->ArgumentExists("-report")) {
      std::string fileName;
      parser->GetCommandLineArgument("-report", fileName);
      std::cout << "print report to the file: " << fileName << std::endl;
      metrics->PrintReportToFile(fileName, getBaseNameFromPath(outputSurfaceFile));
    }
  }

  return EXIT_SUCCESS;
}

