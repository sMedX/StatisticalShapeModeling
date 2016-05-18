#include <boost/filesystem.hpp>
#include <itkStandardMeshRepresenter.h>

#include "ShapeModelToSurfaceRegistrationFilter.h"
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

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string modelFile;
  parser->GetCommandLineArgument("-model", modelFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  double mscale = 1;
  parser->GetCommandLineArgument("-mscale", mscale);

  double regularization = 0.1;
  parser->GetCommandLineArgument("-regularization", regularization);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  std::cout << std::endl;
  std::cout << "shape model to image registration" << std::endl;
  std::cout << "         model file " << modelFile << std::endl;
  std::cout << " input surface file " << surfaceFile << std::endl;
  std::cout << "output surface file " << outputFile << std::endl;
  std::cout << "        model scale " << mscale << std::endl;
  std::cout << "     regularization " << regularization << std::endl;
  std::cout << "         iterations " << numberOfIterations << std::endl;
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

  // read statistical shape model
  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  RepresenterType::Pointer representer = RepresenterType::New();

  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model = StatisticalModelType::New();

  try {
    model->Load(representer, modelFile.c_str());
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "input model " << modelFile << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "    number of points " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  //shape model to image registration
  typedef ShapeModelToSurfaceRegistrationFilter<MeshType> ShapeModelToSurfaceRegistrationFilterType;
  ShapeModelToSurfaceRegistrationFilterType::Pointer shapeModelToSurfaceRegistration = ShapeModelToSurfaceRegistrationFilterType::New();
  shapeModelToSurfaceRegistration->SetNumberOfIterations(numberOfIterations);
  shapeModelToSurfaceRegistration->SetModelScale(mscale);
  shapeModelToSurfaceRegistration->SetRegularizationParameter(regularization);
  shapeModelToSurfaceRegistration->SetShapeModel(model);
  shapeModelToSurfaceRegistration->SetInput(surface);

  try {
    shapeModelToSurfaceRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  shapeModelToSurfaceRegistration->PrintReport(std::cout);

  // write surface
  if (!writeMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // write moved surface
  if (parser->ArgumentExists("-moved")) {
    std::string outputFile;
    parser->GetCommandLineArgument("-moved", outputFile);

    if (!writeMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput(), outputFile)) {
      return EXIT_FAILURE;
    }
  }

  //define interpolator
  typedef itk::LinearInterpolateImageFunction<ShapeModelToSurfaceRegistrationFilterType::PotentialImageType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  ShapeModelToSurfaceRegistrationFilterType::PotentialImageType::ConstPointer potentialImage = shapeModelToSurfaceRegistration->GetPotentialImage();
  interpolator->SetInputImage(potentialImage);

  //define point set
  typedef itk::PointSet<float, Dimension> PointSetType;
  PointSetType::PointsContainer::Pointer points = shapeModelToSurfaceRegistration->GetOutput()->GetPoints();

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
    float value = shapeModelToSurfaceRegistration->GetOptimizer()->GetValue();

    std::string report;
    parser->GetCommandLineArgument("-report", report);

    std::cout << "write report to the file: " << report << std::endl;

    std::string dlm = ";";
    std::string header = dlm;

    int idx1 = surfaceFile.find_last_of("\\/");
    int idx2 = surfaceFile.find_last_of(".");
    std::string scores = surfaceFile.substr(idx1 + 1, idx2 - idx1 - 1) + dlm;

    header += "cost function" + dlm;
    scores += std::to_string(value) + dlm;

    header += "mean" + dlm;
    scores += std::to_string(mean) + dlm;

    header += "rmse" + dlm;
    scores += std::to_string(rmse) + dlm;

    header += "maxd" + dlm;
    scores += std::to_string(maxd) + dlm;

    header += "regularization" + dlm;
    scores += std::to_string(shapeModelToSurfaceRegistration->GetRegularizationParameter()) + dlm;

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


