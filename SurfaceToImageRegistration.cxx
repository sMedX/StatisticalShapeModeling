#include <boost/filesystem.hpp>
#include "SurfaceToImageRegistrationFilter.h"
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

  std::string potentalFile;
  parser->GetCommandLineArgument("-potential", potentalFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << "  input surface file " << surfaceFile << std::endl;
  std::cout << "input potential file " << potentalFile << std::endl;
  std::cout << " output surface file " << outputFile << std::endl;
  std::cout << "          iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read mesh
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
  if (!readImage<FloatImageType>(potential, potentalFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input potential image " << potentalFile << std::endl;
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
    std::string outputFile;
    parser->GetCommandLineArgument("-transform", outputFile);

    typedef SurfaceToImageRegistrationFilterType::TransformType TransformType;
    if (!writeTransform<TransformType>(surfaceToImageRegistration->GetTransform(), outputFile)) {
      return EXIT_FAILURE;
    }
  }

  //define interpolator
  typedef itk::LinearInterpolateImageFunction<SurfaceToImageRegistrationFilterType::PotentialImageType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  SurfaceToImageRegistrationFilterType::PotentialImageType::ConstPointer potentialImage = surfaceToImageRegistration->GetPotentialImage();
  interpolator->SetInputImage(potentialImage);

  //define point set
  typedef itk::PointSet<float, Dimension> PointSetType;
  PointSetType::PointsContainer::Pointer points = surfaceToImageRegistration->GetOutput()->GetPoints();

  size_t numberOfPoints = 0;
  double mean = 0;
  double rmse = 0;
  double maxd = 0;

  for( auto it = points->Begin(); it != points->End(); ++it) {
    itk::ContinuousIndex<double, Dimension> index;
    bool isInside = potentialImage->TransformPhysicalPointToContinuousIndex(it.Value(), index);

    if ( isInside ) {
      double value = std::abs(interpolator->EvaluateAtContinuousIndex(index));

      numberOfPoints++;
      mean += value;
      rmse += value*value;
      maxd = std::max(maxd, value);
    }
  }

  mean = mean / numberOfPoints;
  rmse = std::sqrt(rmse / numberOfPoints);

  std::cout << "distance metrics" << std::endl;
  std::cout << "   mean " << mean << std::endl;
  std::cout << "    rms " << rmse << std::endl;
  std::cout << "maximal " << maxd << std::endl;

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    float value = surfaceToImageRegistration->GetOptimizer()->GetValue();

    std::string report;
    parser->GetCommandLineArgument("-report", report);

    std::cout << "write report to the file: " << report << std::endl;

    std::string dlm = ";";
    std::string header = dlm;

    int idx1 = potentalFile.find_last_of("\\/");
    int idx2 = potentalFile.find_last_of(".");
    std::string scores = potentalFile.substr(idx1 + 1, idx2 - idx1 - 1) + dlm;

    header += "cost function" + dlm;
    scores += std::to_string(value) + dlm;

    header += "mean" + dlm;
    scores += std::to_string(mean) + dlm;

    header += "rmse" + dlm;
    scores += std::to_string(rmse) + dlm;

    header += "maxd" + dlm;
    scores += std::to_string(maxd) + dlm;

    bool exist = boost::filesystem::exists(report);
    std::ofstream ofile;
    ofile.open(report, std::ofstream::out | std::ofstream::app);

    if (!exist) {
      ofile << header << std::endl;
    }

    ofile << scores << std::endl;
    ofile.close();
  }

  return EXIT_SUCCESS;
}


