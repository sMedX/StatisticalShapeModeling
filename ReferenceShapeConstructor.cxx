#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageToVTKImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "SurfaceToLevelSetImageFilter.h"
#include "SurfaceToImageRegistrationFilter.h"
#include "PointSetToImageMetrics.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string refFile;
  parser->GetCommandLineArgument("-reference", refFile);

  int numberOfStages = 3;
  parser->GetCommandLineArgument("-stage", numberOfStages);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << "reference shape constructor" << std::endl;
  std::cout << "    list of files " << listFile << std::endl;
  std::cout << " output reference " << refFile << std::endl;
  std::cout << " number of stages " << numberOfStages << std::endl;
  std::cout << "       iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  // read surface in to vector
  StringList listOfFiles = getFileList(listFile);
  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (StringList::const_iterator it = listOfFiles.begin(); it != listOfFiles.end(); ++it) {
    std::string surfaceFile = it->c_str();
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, surfaceFile)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(surfaceFile);

    std::cout << "input surface polydata info " << surfaceFile << std::endl;
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
  EnumTransformType typeOfTransform = EnumTransformType::Rotation;

  for (int stage = 0; stage < numberOfStages; ++stage) {
    if (stage > 0) {
      typeOfTransform = EnumTransformType::Affine;
    }

    std::cout << "stage " << stage << std::endl;
    std::cout << "type of transform " << typeOfTransform << std::endl;
    std::cout << std::endl;

    // allocate image to update reference image
    FloatImageType::Pointer updateReference = FloatImageType::New();
    updateReference->SetRegions(reference->GetLargestPossibleRegion());
    updateReference->Allocate();
    updateReference->FillBuffer(0);
    updateReference->CopyInformation(reference);

    for (int count= 0; count < vectorOfSurfaces.size(); ++count) {
      std::cout << "stage " << stage << ", " << count << ", file " << vectorOfFiles[count] << std::endl;

      // perform surface to image registration
      SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
      surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
      surfaceToImageRegistration->SetInput(vectorOfSurfaces[count]);
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
      vectorOfSurfaces[count] = surfaceToImageRegistration->GetOutput();

      // compute level set image
      SurfaceToLevelSetImageFilterType::Pointer levelset = SurfaceToLevelSetImageFilterType::New();
      levelset->SetOrigin(reference->GetOrigin());
      levelset->SetSpacing(reference->GetSpacing());
      levelset->SetSize(reference->GetLargestPossibleRegion().GetSize());
      levelset->SetInput(vectorOfSurfaces[count]);

      try {
        levelset->Update();
      }
      catch (itk::ExceptionObject& excep) {
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
      }

      // add current levelset to update reference image
      double const1 = 1. - 1. / (count + 1);
      double const2 = 1. / (count + 1);

      typedef itk::ImageRegionIterator<FloatImageType> IteratorType;
      IteratorType it1(updateReference, updateReference->GetLargestPossibleRegion());
      IteratorType it2(levelset->GetOutput(), levelset->GetOutput()->GetLargestPossibleRegion());

      for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
        double pixelValue = const1 * it1.Get() + const2*it2.Get();
        it1.Set(pixelValue);
      }
    }

    //update reference image
    reference = updateReference;
  }

  //----------------------------------------------------------------------------
  // write reference level set image
  std::cout << "output reference image " << refFile << std::endl;
  std::cout << "   size " << reference->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << reference->GetSpacing() << std::endl;
  std::cout << " origin " << reference->GetOrigin() << std::endl;
  std::cout << std::endl;

  if (!writeImage<FloatImageType>(reference, refFile)) {
    return EXIT_FAILURE;
  }

  // write alignment surfaces
  for (int count = 0; count < vectorOfSurfaces.size(); ++count) {
    std::string fileName = getDirectoryFromPath(refFile) + "/" + getFileNameFromPath(vectorOfFiles[count]);

    std::cout << "output surface polydata info" << std::endl;
    std::cout << fileName << std::endl;
    std::cout << " number of cells " << vectorOfSurfaces[count] ->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << vectorOfSurfaces[count]->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;

    if (!writeMesh<MeshType>(vectorOfSurfaces[count], fileName)) {
      EXIT_FAILURE;
    }
  }

  //----------------------------------------------------------------------------
  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    std::string reportFile;
    parser->GetCommandLineArgument("-report", reportFile);

    std::cout << "write report to the file: " << reportFile << std::endl;

    //open file and write header
    std::ofstream ofile;
    ofile.open(reportFile, std::ofstream::out);

    std::string dlm = ";";
    std::string header = dlm;
    header += "Mean" + dlm;
    header += "RMSE" + dlm;
    header += "Maximal" + dlm;

    ofile << header << std::endl;

    for (int count = 0; count < vectorOfSurfaces.size(); ++count) {

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

      // write metrics to the file
      std::string scores = getFileNameFromPath(vectorOfFiles[count]) + dlm;
      scores += std::to_string(metrics->GetMeanValue()) + dlm;
      scores += std::to_string(metrics->GetRMSEValue()) + dlm;
      scores += std::to_string(metrics->GetMaximalValue()) + dlm;

      ofile << scores << std::endl;
    }

    ofile.close();
  }

  return EXIT_SUCCESS;
}


