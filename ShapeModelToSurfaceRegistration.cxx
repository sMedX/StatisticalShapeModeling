#include <boost/filesystem.hpp>
#include <itkStandardMeshRepresenter.h>

#include "ShapeModelToSurfaceRegistrationFilter.h"
#include "utils/PointSetToImageMetrics.h"
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
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << " shape model to image registration" << std::endl;
  std::cout << "          model file " << modelFile << std::endl;
  std::cout << "  input surface file " << surfaceFile << std::endl;
  std::cout << " output surface file " << outputFile << std::endl;
  std::cout << "         model scale " << mscale << std::endl;
  std::cout << "      regularization " << regularization << std::endl;
  std::cout << "number of iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read surface
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "   input surface polydata " << surfaceFile << std::endl;
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

  //Compute metrics
  typedef ShapeModelToSurfaceRegistrationFilterType::LevelsetImageType LevelsetImageType;
  typedef ShapeModelToSurfaceRegistrationFilterType::PointSetType PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->SetPoints(shapeModelToSurfaceRegistration->GetOutput()->GetPoints());

  typedef PointSetToImageMetrics<PointSetType, LevelsetImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetMovingImage(shapeModelToSurfaceRegistration->GetLevelsetImage());
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    float value = shapeModelToSurfaceRegistration->GetOptimizer()->GetValue();

    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);

    std::cout << "write report to the file: " << fileName << std::endl;

    std::string dlm = ";";
    std::string header = dlm;

    std::string scores = getFileNameFromPath(surfaceFile);

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


