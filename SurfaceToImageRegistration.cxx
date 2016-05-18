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

  std::string potentialFile;
  parser->GetCommandLineArgument("-potential", potentialFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << "  input surface file " << surfaceFile << std::endl;
  std::cout << "input potential file " << potentialFile << std::endl;
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
  FloatImageType::Pointer potential = FloatImageType::New();
  if (!readImage<FloatImageType>(potential, potentialFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input potential image " << potentialFile << std::endl;
  std::cout << "       size " << potential->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << potential->GetSpacing() << std::endl;
  std::cout << "     origin " << potential->GetOrigin() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  //perform GP model to image registration
  typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
  SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
  surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
  surfaceToImageRegistration->SetInput(surface);
  surfaceToImageRegistration->SetPotentialImage(potential);

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
  typedef SurfaceToImageRegistrationFilterType::PotentialImageType PotentialImageType;
  typedef SurfaceToImageRegistrationFilterType::PointSetType PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->SetPoints(surfaceToImageRegistration->GetOutput()->GetPoints());

  typedef PointSetToImageMetrics<PointSetType, PotentialImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetMovingImage(surfaceToImageRegistration->GetPotentialImage());
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


