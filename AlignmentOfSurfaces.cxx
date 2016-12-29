#include <itkMesh.h>
#include <itkImage.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"
#include "ssmSurfaceToLevelSetImageFilter.h"
#include "itkSurfaceToImageRegistrationMethod.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;
using fp = boost::filesystem::path;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-output", surfaceFile);

  int numberOfStages = 3;
  parser->GetCommandLineArgument("-stages", numberOfStages);

  int transform = 1;
  parser->GetCommandLineArgument("-transform", transform);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

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
  FloatImageType::SpacingType spacing;
  spacing.Fill(1);

  typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
  SurfaceToLevelSetImageFilterType::Pointer surfaceToLevelSetImage = SurfaceToLevelSetImageFilterType::New();
  surfaceToLevelSetImage->SetMargin(1.0);
  surfaceToLevelSetImage->SetSpacing(spacing);
  surfaceToLevelSetImage->SetInput(vectorOfSurfaces[0]);
  try {
    surfaceToLevelSetImage->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  FloatImageType::Pointer levelSetImage = surfaceToLevelSetImage->GetOutput();

  //----------------------------------------------------------------------------
  // compute reference image 

  // define types
  typedef itk::SurfaceToImageRegistrationMethod<MeshType> SurfaceToImageRegistrationMethodType;
  typedef SurfaceToImageRegistrationMethodType::EnumTransformType EnumTransformType;
  EnumTransformType typeOfTransform;

  for (int stage = 0; stage < numberOfStages; ++stage) {
    switch (stage) {
    case 0:
      typeOfTransform = EnumTransformType::Translation;
      break;
    case 1:
      typeOfTransform = EnumTransformType::Euler3D;
      break;
    default:
      typeOfTransform = static_cast<EnumTransformType>(transform);
    }

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

      surfaceToImageRegistration->PrintReport(std::cout);

      // update surface in the vector
      vectorOfSurfaces[n] = surfaceToImageRegistration->GetOutput();

      // compute level set image
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
  if (parser->ArgumentExists("-reference")) {
    std::string fileName;
    parser->GetCommandLineArgument("-reference", fileName);

    std::cout << "output reference image " << fileName << std::endl;
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
    fp path = fp(surfaceFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(surfaceFile).filename().string());
    std::string outputSurfaceFile = path.string();

    std::cout << "output surface info" << std::endl;
    std::cout << outputSurfaceFile << std::endl;
    std::cout << " number of cells " << vectorOfSurfaces[count]->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << vectorOfSurfaces[count]->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;

    if (!writeMesh<MeshType>(vectorOfSurfaces[count], outputSurfaceFile)) {
      EXIT_FAILURE;
    }

    // compute metrics
    typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
    PointSetType::Pointer pointSet = PointSetType::New();
    pointSet->SetPoints(vectorOfSurfaces[count]->GetPoints());

    typedef PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
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


