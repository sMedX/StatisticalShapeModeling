#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

#include "utils/PointSetToImageMetrics.h"
#include "itkShapeModelRegistrationMethod.h"
#include "SurfaceToLevelSetImageFilter.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions);

int main(int argc, char** argv)
{
  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string modelFile;
  parser->GetCommandLineArgument("-model", modelFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  unsigned int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  double mscale = 1;
  parser->GetCommandLineArgument("-mscale", mscale);

  std::vector<double> regularization(0.1);
  parser->GetCommandLineArgument("-regularization", regularization);

  std::vector<double> parameters;
  parser->GetCommandLineArgument("-parameters", parameters);

  double scale = 100;
  parser->GetCommandLineArgument("-scale", scale);

  unsigned int degree = 2;
  parser->GetCommandLineArgument("-degree", degree);

  std::cout << std::endl;
  std::cout << " shape model to image registration" << std::endl;
  std::cout << "            model file " << modelFile << std::endl;
  std::cout << "    input surface file " << surfaceFile << std::endl;
  std::cout << "   output surface file " << outputFile << std::endl;
  std::cout << "           model scale " << mscale << std::endl;
  std::cout << "  number of iterations " << numberOfIterations << std::endl;
  std::cout << "                degree " << degree << std::endl;

  unsigned int numberOfStages = parameters.size() + 1;
  for (int n = regularization.size(); n < numberOfStages; ++n) {
    regularization.push_back(*regularization.end());
  }

  std::cout << "        regularization ";
  for (int n = 0; n < regularization.size(); ++n) {
    std::cout << regularization[n] << " ";
  }
  std::cout << std::endl;

  std::cout << "            parameters ";
  for (int n = 0; n < parameters.size(); ++n) {
    std::cout << parameters[n] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read surface
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input surface polydata " << surfaceFile << std::endl;
  std::cout << "number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
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
  // compute level set image
  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
  SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
  levelset->SetMargin(0.10);
  levelset->SetSpacing(0.5);
  levelset->SetInput(surface);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // shape model to image registration
  typedef itk::ShapeModelRegistrationMethod<StatisticalModelType, MeshType> ShapeModelRegistrationMethod;
  ShapeModelRegistrationMethod::Pointer shapeModelToSurfaceRegistration;

  for (int n = 0; n < numberOfStages; ++n) {
    std::cout << "---------- stage (" << n + 1 << " / " << numberOfStages << ") ----------" << std::endl;

    // run registration
    shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetNumberOfIterations(numberOfIterations);
    shapeModelToSurfaceRegistration->SetModelScale(mscale);
    shapeModelToSurfaceRegistration->SetRegularizationParameter(regularization[n]);
    shapeModelToSurfaceRegistration->SetDegree(degree);
    try {
      shapeModelToSurfaceRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    shapeModelToSurfaceRegistration->PrintReport(std::cout);

    // build new model
    if (n + 1 < numberOfStages) {
      MeshType::Pointer reference = const_cast<MeshType*> (shapeModelToSurfaceRegistration->GetOutput());
      model = BuildGPModel(reference, parameters[n], scale, model->GetNumberOfPrincipalComponents());
    }
  }

  //----------------------------------------------------------------------------
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

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);
    std::cout << "print report to the file: " << fileName << std::endl;
    metrics->PrintReportToFile(fileName, getBaseNameFromPath(surfaceFile));
  }

  // write surface
  if (!writeMesh<MeshType>(shapeModelToSurfaceRegistration->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }

  // write levelset image
  if (parser->ArgumentExists("-levelset")) {
    std::string fileName;
    parser->GetCommandLineArgument("-levelset", fileName);
    if (!writeImage(shapeModelToSurfaceRegistration->GetLevelSetImage(), fileName)) {
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions)
{
  // create kernel
  typedef DataTypeShape::PointType PointType;
  auto gaussianKernel = new GaussianKernel<PointType>(parameters);

  typedef std::shared_ptr<const statismo::ScalarValuedKernel<PointType>> MatrixPointerType;
  MatrixPointerType kernel(gaussianKernel);

  typedef std::shared_ptr<statismo::MatrixValuedKernel<PointType>> KernelPointerType;
  KernelPointerType unscaledKernel(new statismo::UncorrelatedMatrixValuedKernel<PointType>(kernel.get(), Dimension));
  KernelPointerType modelBuildingKernel(new statismo::ScaledKernel<PointType>(unscaledKernel.get(), scale));

  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
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


