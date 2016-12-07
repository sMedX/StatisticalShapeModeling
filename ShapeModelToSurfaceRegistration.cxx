#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>

#include "itkShapeModelRegistrationMethod.h"
#include "ShapeModelToSurfaceRegistrationFilter.h"
#include "utils/PointSetToImageMetrics.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer BuildGPShapeModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions);

int main(int argc, char** argv)
{
  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string modelFile;
  parser->GetCommandLineArgument("-model", modelFile);

  std::string referenceSurfaceFile;
  parser->GetCommandLineArgument("-reference", referenceSurfaceFile);

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

  // read reference surface
  MeshType::Pointer referenceSurface = MeshType::New();
  if (!readMesh<MeshType>(referenceSurface, referenceSurfaceFile)) {
    return EXIT_FAILURE;
  };


  //----------------------------------------------------------------------------
  //shape model to image registration
  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
  SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
  levelset->SetMargin(0.3);
  levelset->SetSpacing(1);
  levelset->SetInput(surface);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::ShapeModelRegistrationMethod<StatisticalModelType, MeshType> ShapeModelRegistrationMethod;
  ShapeModelRegistrationMethod::Pointer shapeModelToSurfaceRegistration;
  
  std::vector<double> parameters;
  parameters.push_back(30);
  parameters.push_back(20);
  parameters.push_back(10);
  parameters.push_back(5);
  double scale = 100;
  int numberOfBasisFunctions = 200;

  for (int n = 0; n < parameters.size(); ++n) {
    StatisticalModelType::Pointer model = BuildGPShapeModel(referenceSurface, parameters[n], scale, numberOfBasisFunctions);

    shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetNumberOfIterations(numberOfIterations);
    shapeModelToSurfaceRegistration->SetModelScale(mscale);
    shapeModelToSurfaceRegistration->SetRegularizationParameter(regularization);
    try {
      shapeModelToSurfaceRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    shapeModelToSurfaceRegistration->PrintReport();

    referenceSurface = const_cast<MeshType*> (shapeModelToSurfaceRegistration->GetOutput());
  }

  // write surface
  if (!writeMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // compute metrics
  typedef ShapeModelRegistrationMethod::LevelSetImageType LevelsetImageType;
  typedef ShapeModelRegistrationMethod::PointSetType PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->SetPoints(const_cast<PointSetType::PointsContainer*> (shapeModelToSurfaceRegistration->GetOutput()->GetPoints()));

  typedef PointSetToImageMetrics<PointSetType, LevelsetImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetMovingImage(shapeModelToSurfaceRegistration->GetLevelSetImage());
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write levelset image
  if (parser->ArgumentExists("-levelset")) {
    std::string fileName;
    parser->GetCommandLineArgument("-levelset", fileName);
    if (!writeImage(shapeModelToSurfaceRegistration->GetLevelSetImage(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    float value = shapeModelToSurfaceRegistration->GetOptimizer()->GetValue();

    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);

    std::cout << "write report to the file: " << fileName << std::endl;

    std::string dlm = ";";

    std::string header = dlm;
    std::string scores = getBaseNameFromPath(surfaceFile) + dlm;

    header += "Cost function" + dlm;
    scores += std::to_string(value) + dlm;

    header += "Mean" + dlm;
    scores += std::to_string(metrics->GetMeanValue()) + dlm;

    header += "RMSE" + dlm;
    scores += std::to_string(metrics->GetRMSEValue()) + dlm;

    header += "Quantile " + std::to_string(metrics->GetLevelOfQuantile()) + dlm;
    scores += std::to_string(metrics->GetQuantileValue()) + dlm;

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

StatisticalModelType::Pointer BuildGPShapeModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions)
{
  // create kernel
  typedef DataTypeShape::PointType PointType;
  auto gaussianKernel = new GaussianKernel<PointType>(parameters);

  typedef std::shared_ptr<const statismo::ScalarValuedKernel<PointType>> MatrixPointerType;
  MatrixPointerType kernel(gaussianKernel);

  typedef std::shared_ptr<statismo::MatrixValuedKernel<PointType>> KernelPointerType;
  KernelPointerType unscaledKernel(new statismo::UncorrelatedMatrixValuedKernel<PointType>(kernel.get(), 3));
  KernelPointerType modelBuildingKernel(new statismo::ScaledKernel<PointType>(unscaledKernel.get(), scale));

  typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
  typedef RepresenterType::DatasetPointerType DatasetPointerType;
  RepresenterType::Pointer representer = RepresenterType::New();
  representer->SetReference(surface);

  // build model
  typedef itk::LowRankGPModelBuilder<MeshType> ModelBuilderType;
  ModelBuilderType::Pointer modelBuilder = ModelBuilderType::New();
  modelBuilder->SetRepresenter(representer);

  // build and save model to file
  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model;
  try {
    model = modelBuilder->BuildNewModel(representer->IdentitySample(), *modelBuildingKernel.get(), numberOfBasisFunctions);
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  return model;
}


