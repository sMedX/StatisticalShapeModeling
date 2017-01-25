#include <boost/program_options.hpp>
#include <itkImageMomentsCalculator.h>
#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmShapeModelToImageRegistrationMethod.h"
#include "ssm/ssmBinaryImageToLevelSetImageFilter.h"

int main(int argc, char** argv)
{
  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string modelFile;
  parser->GetCommandLineArgument("-model", modelFile);

  std::string imageFile;
  parser->GetCommandLineArgument("-image", imageFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  size_t typeOfransform = 1;
  parser->GetCommandLineArgument("-transform", typeOfransform);

  size_t numberOfIterations = 500;
  parser->GetCommandLineArgument("-iterations", numberOfIterations);

  double regularization = 0.1;
  parser->GetCommandLineArgument("-regularization", regularization);

  size_t degree = 2;
  parser->GetCommandLineArgument("-degree", degree);

  std::cout << std::endl;
  std::cout << " shape model to image registration" << std::endl;
  std::cout << "            model file " << modelFile << std::endl;
  std::cout << "      input image file " << imageFile << std::endl;
  std::cout << "   output surface file " << outputFile << std::endl;
  std::cout << "  number of iterations " << numberOfIterations << std::endl;
  std::cout << "     type of transform " << typeOfransform << std::endl;
  std::cout << "                degree " << degree << std::endl;
  std::cout << "        regularization " << regularization << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read surface
  BinaryImageType::Pointer image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, imageFile)) {
    return EXIT_FAILURE;
  }
  std::cout << "image " << imageFile << std::endl;
  std::cout << "size    " << image->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << image->GetSpacing() << std::endl;
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

  std::cout << "model " << modelFile << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // compute level set image
  typedef ssm::BinaryImageToLevelSetImageFilter<BinaryImageType, FloatImageType> BinaryImageToLevelSetImageType;
  BinaryImageToLevelSetImageType::Pointer levelset = BinaryImageToLevelSetImageType::New();
  levelset->SetInput(image);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  // initialize spatial transform
  MeshType::BoundingBoxType::ConstPointer boundingBox = model->DrawMean()->GetBoundingBox();
  BinaryImageType::SpacingType spacing(1);
  BinaryImageType::PointType origin = boundingBox->GetMinimum();
  BinaryImageType::SizeType size;
  for (size_t n = 0; n < Dimension; ++n) {
    size[n] = (boundingBox->GetMaximum()[n] - boundingBox->GetMinimum()[n]) / spacing[n];
  }

  typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> ShapeToBinaryImageFilterType;
  ShapeToBinaryImageFilterType::Pointer shapeToImage = ShapeToBinaryImageFilterType::New();
  shapeToImage->SetInput(model->DrawMean());
  shapeToImage->SetSize(size);
  shapeToImage->SetOrigin(origin);
  shapeToImage->SetSpacing(spacing);
  shapeToImage->SetOutsideValue(0);
  shapeToImage->SetInsideValue(1);
  try {
    shapeToImage->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }

  // moment calculators
  typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
  ImageCalculatorType::Pointer movingCalculator = ImageCalculatorType::New();
  movingCalculator->SetImage(shapeToImage->GetOutput());
  movingCalculator->Compute();

  ImageCalculatorType::Pointer fixedCalculator = ImageCalculatorType::New();
  fixedCalculator->SetImage(image);
  fixedCalculator->Compute();

  typedef ImageCalculatorType::VectorType VectorType;
  VectorType center = movingCalculator->GetCenterOfGravity();
  VectorType translation = fixedCalculator->GetCenterOfGravity() - movingCalculator->GetCenterOfGravity();

  typedef ssm::TransformInitializer<double> TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTypeOfTransform(typeOfransform);
  initializer->SetCenter(center);
  initializer->SetTranslation(translation);
  try {
    initializer->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    return EXIT_FAILURE;
  }
  initializer->PrintReport(std::cout);

  //----------------------------------------------------------------------------
  // shape model to image registration
  typedef ssm::ShapeModelToImageRegistrationMethod<StatisticalModelType, MeshType> ShapeModelRegistrationMethod;
  ShapeModelRegistrationMethod::Pointer shapeModelToSurfaceRegistration;

  // perform registration
  shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
  shapeModelToSurfaceRegistration->SetShapeModel(model);
  shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
  shapeModelToSurfaceRegistration->SetNumberOfIterations(numberOfIterations);
  shapeModelToSurfaceRegistration->SetRegularizationParameter(regularization);
  shapeModelToSurfaceRegistration->SetDegree(degree);
  shapeModelToSurfaceRegistration->SetTransformInitializer(initializer);
  try {
    shapeModelToSurfaceRegistration->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  shapeModelToSurfaceRegistration->PrintReport(std::cout);

  typedef std::pair<std::string, std::string> PairType;
  std::vector<PairType> info;
  info.push_back(PairType("Metric", std::to_string(shapeModelToSurfaceRegistration->GetOptimizer()->GetValue())));

  //----------------------------------------------------------------------------
  // compute metrics
  typedef ShapeModelRegistrationMethod::LevelSetImageType LevelsetImageType;
  typedef ShapeModelRegistrationMethod::PointSetType PointSetType;
  PointSetType::Pointer points = PointSetType::New();
  points->SetPoints(const_cast<PointSetType::PointsContainer*> (shapeModelToSurfaceRegistration->GetOutput()->GetPoints()));

  typedef ssm::PointSetToImageMetrics<PointSetType, LevelsetImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(points);
  metrics->SetMovingImage(shapeModelToSurfaceRegistration->GetLevelSetImage());
  metrics->SetInfo(info);
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);
    std::cout << "print report to the file: " << fileName << std::endl;
    metrics->PrintReportToFile(fileName, getBaseNameFromPath(imageFile));
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
