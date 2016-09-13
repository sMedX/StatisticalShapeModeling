#include <itkMesh.h>
#include <itkImage.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"
#include "SurfaceToLevelSetImageFilter.h"
#include "SurfaceToImageRegistrationFilter.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  int numberOfStages = 3;
  parser->GetCommandLineArgument("-stage", numberOfStages);

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

  // read surface in to vector
  StringList listOfFiles = getFileList(listFile);
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
  // compute initial reference level set image
  FloatImageType::SpacingType spacing;
  spacing.Fill(1);

  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilterType;
  SurfaceToLevelSetImageFilterType::Pointer levelset = SurfaceToLevelSetImageFilterType::New();
  levelset->SetMargin(0.5);
  levelset->SetSpacing(spacing);
  levelset->SetInput(vectorOfSurfaces[0]);

  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  FloatImageType::Pointer reference = levelset->GetOutput();

  //----------------------------------------------------------------------------
  // compute reference image 

  // define types
  typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
  typedef SurfaceToImageRegistrationFilterType::EnumTransformType EnumTransformType;
  EnumTransformType typeOfTransform;

  for (int stage = 0; stage < numberOfStages; ++stage) {
    switch ( stage ) {
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
    FloatImageType::Pointer updateReference = FloatImageType::New();
    updateReference->SetRegions(reference->GetLargestPossibleRegion());
    updateReference->Allocate();
    updateReference->FillBuffer(0);
    updateReference->CopyInformation(reference);

    for (int n = 0; n < vectorOfSurfaces.size(); ++n) {
      std::cout << "stage " << stage + 1 << "/" << numberOfStages << ", surface " << n + 1 << "/" << vectorOfSurfaces.size() << ", " << vectorOfFiles[n] << std::endl;

      // perform surface to image registration
      SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
      surfaceToImageRegistration->SetInput(vectorOfSurfaces[n]);
      surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
      surfaceToImageRegistration->SetTypeOfTransform(typeOfTransform);
      surfaceToImageRegistration->SetLevelsetImage(reference);

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
      SurfaceToLevelSetImageFilterType::Pointer levelset = SurfaceToLevelSetImageFilterType::New();
      levelset->SetOrigin(reference->GetOrigin());
      levelset->SetSpacing(reference->GetSpacing());
      levelset->SetSize(reference->GetLargestPossibleRegion().GetSize());
      levelset->SetInput(vectorOfSurfaces[n]);

      try {
        levelset->Update();
      }
      catch (itk::ExceptionObject& excep) {
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
      }

      // add current levelset to update reference image
      typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
      FilterType::Pointer multiply = FilterType::New();
      multiply->SetInput(levelset->GetOutput());
      multiply->SetConstant(1 / (double) vectorOfSurfaces.size());
      multiply->Update();

      typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
      AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
      addfilter->SetInput1(updateReference);
      addfilter->SetInput2(multiply->GetOutput());
      addfilter->Update();

      updateReference = addfilter->GetOutput();
    }

    //update reference image
    reference = updateReference;
  }

  //----------------------------------------------------------------------------
  // write reference level set image
  if (parser->ArgumentExists("-potential")) {
    std::string fileName;
    parser->GetCommandLineArgument("-potential", fileName);

    std::cout << "output reference image " << fileName << std::endl;
    std::cout << "   size " << reference->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "spacing " << reference->GetSpacing() << std::endl;
    std::cout << " origin " << reference->GetOrigin() << std::endl;
    std::cout << std::endl;

    if (!writeImage<FloatImageType>(reference, fileName)) {
      return EXIT_FAILURE;
    }
  }

  std::string reportFileName;
  std::ofstream rfile;

  // open file to write report
  if (parser->ArgumentExists("-report")) {
    parser->GetCommandLineArgument("-report", reportFileName);
    std::cout << "write report to the file: " << reportFileName << std::endl;
    rfile.open(reportFileName, std::ofstream::out);
  }

  // write alignment surfaces
  for (int count = 0; count < vectorOfSurfaces.size(); ++count) {

    // define full file name for output surface
    std::string fileName = getDirectoryFromPath(surfaceFile) + getBaseNameFromPath(vectorOfFiles[count]) + "-" + getFileNameFromPath(surfaceFile);

    std::cout << "output surface polydata info" << std::endl;
    std::cout << fileName << std::endl;
    std::cout << " number of cells " << vectorOfSurfaces[count]->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << vectorOfSurfaces[count]->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;

    if (!writeMesh<MeshType>(vectorOfSurfaces[count], fileName)) {
      EXIT_FAILURE;
    }

    // compute metrics
    typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
    PointSetType::Pointer pointSet = PointSetType::New();
    pointSet->SetPoints(vectorOfSurfaces[count]->GetPoints());

    typedef PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
    PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
    metrics->SetFixedPointSet(pointSet);
    metrics->SetMovingImage(reference);
    metrics->Compute();
    metrics->PrintReport(std::cout);

    // write metrics to report file
    if ( rfile.is_open() ) {
      std::string dlm = ";";

      if (count == 0) {
        std::string header = dlm;
        header += "Mean" + dlm;
        header += "RMSE" + dlm;
        header += "Quantile " + std::to_string(metrics->GetLevelOfQuantile()) + dlm;
        header += "Maximal" + dlm;

        rfile << header << std::endl;
      }

      std::string scores = getFileNameFromPath(fileName) + dlm;
      scores += std::to_string(metrics->GetMeanValue()) + dlm;
      scores += std::to_string(metrics->GetRMSEValue()) + dlm;
      scores += std::to_string(metrics->GetQuantileValue()) + dlm;
      scores += std::to_string(metrics->GetMaximalValue()) + dlm;

      rfile << scores << std::endl;
    }
  }

  rfile.close();

  return EXIT_SUCCESS;
}


