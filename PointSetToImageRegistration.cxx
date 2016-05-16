#include <boost/filesystem.hpp>
#include "PointSetToImageRegistration.h"
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

  std::string labelFile;
  parser->GetCommandLineArgument("-image", labelFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  std::string transformFile;
  parser->GetCommandLineArgument("-transform", transformFile);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << " input surface file " << surfaceFile << std::endl;
  std::cout << "   input image file " << labelFile << std::endl;
  std::cout << "output surface file " << outputFile << std::endl;
  std::cout << "     transform file " << transformFile << std::endl;
  std::cout << "         iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read mesh
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "   input surface " << outputFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  // read image
  BinaryImageType::Pointer label = BinaryImageType::New();
  if (!readImage<BinaryImageType>(label, labelFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input image " << labelFile << std::endl;
  std::cout << "       size " << label->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << label->GetSpacing() << std::endl;
  std::cout << "     origin " << label->GetOrigin() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  //perform GP model to image registration
  typedef PointSetToImageRegistration<MeshType> ShapeModelToImageRegistrationType;
  ShapeModelToImageRegistrationType::Pointer shapeModelToImageRegistration = ShapeModelToImageRegistrationType::New();
  shapeModelToImageRegistration->SetNumberOfIterations(numberOfIterations);
  shapeModelToImageRegistration->SetInput(surface);
  shapeModelToImageRegistration->SetLabel(label);

  try {
    shapeModelToImageRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  shapeModelToImageRegistration->PrintReport(std::cout);

  // write output mesh
  if (!writeMesh<MeshType>(shapeModelToImageRegistration->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // write transform
  typedef ShapeModelToImageRegistrationType::TransformType TransformType;
  if (!writeTransform<TransformType>(shapeModelToImageRegistration->GetTransform(), transformFile)) {
    return EXIT_FAILURE;
  }

  //define interpolator
  typedef itk::LinearInterpolateImageFunction<ShapeModelToImageRegistrationType::PotentialImageType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  ShapeModelToImageRegistrationType::PotentialImageType::ConstPointer potentialImage = shapeModelToImageRegistration->GetPotentialImage();
  interpolator->SetInputImage(potentialImage);

  //define point set
  typedef itk::PointSet<float, Dimension> PointSetType;
  PointSetType::PointsContainer::Pointer points = shapeModelToImageRegistration->GetOutput()->GetPoints();

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
    float value = shapeModelToImageRegistration->GetOptimizer()->GetValue();

    std::string report;
    parser->GetCommandLineArgument("-report", report);

    std::cout << "write report to the file: " << report << std::endl;

    std::string dlm = ";";
    std::string header = dlm;

    int idx1 = labelFile.find_last_of("\\/");
    int idx2 = labelFile.find_last_of(".");
    std::string scores = labelFile.substr(idx1 + 1, idx2 - idx1 - 1) + dlm;

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


