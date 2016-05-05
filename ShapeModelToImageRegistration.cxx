#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <boost/filesystem.hpp>
#include <itkStandardMeshRepresenter.h>

#include "ShapeModelToImageRegistration.h"
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

  std::string imageFile;
  parser->GetCommandLineArgument("-image", imageFile);

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
  std::cout << "     image file " << imageFile << std::endl;
  std::cout << "     model file " << modelFile << std::endl;
  std::cout << "    output file " << outputFile << std::endl;
  std::cout << "    model scale " << mscale << std::endl;
  std::cout << " regularization " << regularization << std::endl;
  std::cout << "     iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image

  BinaryImageType::Pointer image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, imageFile)) {
    return EXIT_FAILURE;
  }

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

  //----------------------------------------------------------------------------
  //perform GP model to image registration
  typedef ShapeModelToImageRegistration<BinaryImageType, MeshType> ShapeModelToImageRegistrationType;
  ShapeModelToImageRegistrationType::Pointer shapeModelToImageRegistration = ShapeModelToImageRegistrationType::New();
  shapeModelToImageRegistration->SetShapeModel(model);
  shapeModelToImageRegistration->SetModelScale(mscale);
  shapeModelToImageRegistration->SetRegularizationParameter(regularization);
  shapeModelToImageRegistration->SetInput(image);
  shapeModelToImageRegistration->SetNumberOfIterations(numberOfIterations);

  try {
    shapeModelToImageRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  shapeModelToImageRegistration->PrintReport(std::cout);

  // write surface
  if (!writeMesh<MeshType>(shapeModelToImageRegistration->GetDeformedOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // write posterior surface
  if (parser->ArgumentExists("-deformed")) {
    std::string outputFile;
    parser->GetCommandLineArgument("-deformed", outputFile);

    if (!writeMesh<MeshType>(shapeModelToImageRegistration->GetOutput(), outputFile)) {
      return EXIT_FAILURE;
    }
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
      double value = interpolator->EvaluateAtContinuousIndex(index);

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

    int idx1 = imageFile.find_last_of("\\/");
    int idx2 = imageFile.find_last_of(".");
    std::string scores = imageFile.substr(idx1 + 1, idx2 - idx1 - 1) + dlm;

    header += "cost function" + dlm;
    scores += std::to_string(value) + dlm;

    header += "mean" + dlm;
    scores += std::to_string(mean) + dlm;

    header += "rmse" + dlm;
    scores += std::to_string(rmse) + dlm;

    header += "maxd" + dlm;
    scores += std::to_string(maxd) + dlm;

    header += "regularization" + dlm;
    scores += std::to_string(shapeModelToImageRegistration->GetRegularizationParameter()) + dlm;

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


