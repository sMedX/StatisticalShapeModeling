#include <boost/filesystem.hpp>

#include "SurfaceToImageRegistrationFilter.h"
#include "PointSetToImageMetrics.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv)
{
  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string levelsetFile;
  parser->GetCommandLineArgument("-levelset", levelsetFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << "  input surface file " << surfaceFile << std::endl;
  std::cout << " input levelset file " << levelsetFile << std::endl;
  std::cout << " output surface file " << outputFile << std::endl;
  std::cout << "          iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read surface
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "   input surface polydata " << outputFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  // read image
  FloatImageType::Pointer levelset = FloatImageType::New();
  if (!readImage<FloatImageType>(levelset, levelsetFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input level set image " << levelsetFile << std::endl;
  std::cout << "       size " << levelset->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << levelset->GetSpacing() << std::endl;
  std::cout << "     origin " << levelset->GetOrigin() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  //perform surface to image registration
  typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
  typedef SurfaceToImageRegistrationFilterType::EnumTransformType EnumTransformType;
  EnumTransformType typeOfTransform = EnumTransformType::Affine;

  SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
  surfaceToImageRegistration->SetTypeOfTransform(typeOfTransform);
  surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
  surfaceToImageRegistration->SetInput(surface);
  surfaceToImageRegistration->SetLevelsetImage(levelset);

  try {
    surfaceToImageRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  surfaceToImageRegistration->PrintReport(std::cout);

  // write output mesh
  if (!writeMesh<MeshType>(surfaceToImageRegistration->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // write transform
  if (parser->ArgumentExists("-transform")) {
    std::string fileName;
    parser->GetCommandLineArgument("-transform", fileName);

    typedef SurfaceToImageRegistrationFilterType::TransformType TransformType;
    if (!writeTransform<TransformType>(surfaceToImageRegistration->GetTransform(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  //Compute metrics
  typedef SurfaceToImageRegistrationFilterType::LevelsetImageType LevelsetImageType;
  typedef SurfaceToImageRegistrationFilterType::PointSetType PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->SetPoints(surfaceToImageRegistration->GetOutput()->GetPoints());

  typedef PointSetToImageMetrics<PointSetType, LevelsetImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetMovingImage(surfaceToImageRegistration->GetLevelsetImage());
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    float value = surfaceToImageRegistration->GetOptimizer()->GetValue();

    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);

    std::cout << "write report to the file: " << fileName << std::endl;

    std::string dlm = ";";
    std::string header = dlm;

    int idx1 = surfaceFile.find_last_of("\\/");
    int idx2 = surfaceFile.find_last_of(".");
    std::string scores = surfaceFile.substr(idx1 + 1, idx2 - idx1 - 1) + dlm;

    header += "Cost function" + dlm;
    scores += std::to_string(value) + dlm;

    header += "Mean" + dlm;
    scores += std::to_string(metrics->GetMeanValue()) + dlm;

    header += "RMSE" + dlm;
    scores += std::to_string(metrics->GetRMSEValue()) + dlm;

    header += "Maximal" + dlm;
    scores += std::to_string(metrics->GetMaximalValue()) + dlm;

    bool exist = boost::filesystem::exists(fileName);
    std::ofstream ofile;
    ofile.open(fileName, std::ofstream::out | std::ofstream::app);

    if (!exist) {
      ofile << header << std::endl;
    }

    ofile << scores << std::endl;
    ofile.close();
  }

  return EXIT_SUCCESS;
}


