#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmMeshToLevelSetImageFilter.h"
#include "ssmMeshToImageRegistrationMethod.h"
#include "ssmMeshPropertiesCalculator.h"
#include "ssmAlignmentOptions.h"

typedef std::vector<MeshType::Pointer> MeshVectorType;
typedef itk::Transform<double, MeshType::PointDimension> TransformType;
typedef std::vector<TransformType::ConstPointer> TransformVectorType;
FloatImageType::Pointer computeLevelSetImage(FloatImageType::Pointer levelSetImage, MeshVectorType & vectorOfSurfaces, TransformVectorType & vectorOfTransform, size_t typeOfransform);

typedef boost::filesystem::path fp;

int main(int argc, char** argv) {

  ssm::AlignmentOptions options;
  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  // read options from config file
  if (!options.ParseConfigFile()) {
    return EXIT_FAILURE;
  }
  options.PrintConfig();

  // read list of files
  StringVector listOfInputFiles;
  try {
    listOfInputFiles = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  MeshVectorType vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (const auto & fileName : listOfInputFiles) {
    auto surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    printMeshInfo<MeshType>(surface, fileName);

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
    auto boundingbox = vectorOfSurfaces[count]->GetBoundingBox();

    // compute mask for surface
    BinaryImageType::SpacingType spacing(1);
    BinaryImageType::SizeType size;

    for (unsigned i = 0; i < Dimension; ++i) {
      spacing[i] = 1;
      size[i] = (boundingbox->GetMaximum()[i] - boundingbox->GetMinimum()[i]) / spacing[i];
    }

    typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> TriangleMeshToBinaryImageFilterType;
    auto surfaceToMask = TriangleMeshToBinaryImageFilterType::New();
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
    auto calculator = ImageCalculatorType::New();
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

  auto levelSetImage = FloatImageType::New();
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
    auto transform = TransformType::New();
    transform->SetOffset(centerOfMaximalBoundingBox - vectorOfCenters[count]);

    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformFilterType;
    auto transformSurface = TransformFilterType::New();
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
    typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilterType;
    auto surfaceToLevelSetImage = MeshToLevelSetImageFilterType::New();
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
    auto multiply = FilterType::New();
    multiply->SetInput(surfaceToLevelSetImage->GetOutput());
    multiply->SetConstant(1 / (double)vectorOfSurfaces.size());

    typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
    auto add = AddImageFilterType::New();
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

  for (size_t stage = 0; stage < options.GetNumberOfStages(); ++stage) {
    for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
      std::cout << "stage " << stage + 1 << "/" << options.GetNumberOfStages() << ", surface " << count + 1 << "/" << vectorOfSurfaces.size() << ", " << vectorOfFiles[count] << std::endl;

      // perform surface to image registration
      typedef ssm::MeshToImageRegistrationMethod<MeshType> MeshToImageRegistrationMethodType;
      auto surfaceToImageRegistration = MeshToImageRegistrationMethodType::New();
      surfaceToImageRegistration->SetInput(vectorOfSurfaces[count]);
      surfaceToImageRegistration->SetNumberOfIterations(options.GetNumberOfIterations());
      surfaceToImageRegistration->SetTypeOfTransform(options.GetTransform());
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

    // update reference level set image
    levelSetImage = computeLevelSetImage(levelSetImage, vectorOfSurfaces, vectorOfTransforms, options.GetTransform());
  }

  //----------------------------------------------------------------------------
  // write reference level set image
  typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
  auto multiply = FilterType::New();
  multiply->SetInput(levelSetImage);
  multiply->SetConstant(-1);

  if (!writeImage<FloatImageType>(multiply->GetOutput(), options.GetReferenceFileName())) {
    return EXIT_FAILURE;
  }
  printImageInfo<FloatImageType>(multiply->GetOutput(), options.GetReferenceFileName());

  //----------------------------------------------------------------------------
  // write aligned surfaces
  StringVector listOfOutputFiles;

  for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
    const auto & fileName = options.FormatOutput(vectorOfFiles[count]);
    listOfOutputFiles.push_back(fileName);

    std::cout << "output file " << fileName << std::endl;
    if (!writeMesh<MeshType>(vectorOfSurfaces[count], fileName)) {
      EXIT_FAILURE;
    }

    // compute metrics
    typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
    typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
    auto metrics = PointSetToImageMetricsType::New();
    metrics->SetPointSetAsMesh<MeshType>(vectorOfSurfaces[count]);
    metrics->SetImage(levelSetImage);
    metrics->Compute();
    metrics->PrintReport(std::cout);

    // print report to *.csv file
    std::cout << "print report to the file: " << options.GetReportFile() << std::endl;
    std::cout << std::endl;

    metrics->PrintReportToFile(options.GetReportFile(), getBaseNameFromPath(fileName));
  }

  // write list of files
  try {
    writeListToFile(options.GetOutputList(), listOfOutputFiles);
  }
  catch (std::ofstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
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
      auto surface = vectorOfSurfaces[count];
      auto transform = vectorOfTransform[count];

      typedef ssm::MeshPropertiesCalculator<MeshType> MeshMomentsCalculatorType;
      auto calculator = MeshMomentsCalculatorType::New();
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
    auto transform = const_cast<TransformType*>(vectorOfTransform[count].GetPointer());

    if (typeOfransform > 1) {
      auto parameters = transform->GetParameters();
      parameters[index] *= scale;
      transform->SetParameters(parameters);
    }

    // transform surface
    typedef itk::TransformMeshFilter<MeshType, MeshType, TransformType> TransformFilterType;
    auto transformMeshFilter = TransformFilterType::New();
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
    typedef ssm::MeshToLevelSetImageFilter<MeshType, FloatImageType> MeshToLevelSetImageFilterType;
    auto surfaceToLevelSetImage = MeshToLevelSetImageFilterType::New();
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
    auto multiply = FilterType::New();
    multiply->SetInput(surfaceToLevelSetImage->GetOutput());
    multiply->SetConstant(1 / (double)vectorOfSurfaces.size());

    typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
    auto add = AddImageFilterType::New();
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
