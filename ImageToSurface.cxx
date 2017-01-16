#include <itkMesh.h>
#include <vtkPolyData.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <itkImageToVTKImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;


double AverageLengthOfEdges(vtkPolyData*poly);
double AverageAreaOfCells(vtkPolyData*poly);
FloatImageType::Pointer ComputeDistanceMapImage(FloatImageType::Pointer image, float background, float foreground);

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string imageFile;
  parser->GetCommandLineArgument("-image", imageFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-output", surfaceFile);

  double sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  double relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  unsigned int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  unsigned int numberOfPoints = 0;
  parser->GetCommandLineArgument("-points", numberOfPoints);

  bool isbinary = true;
  parser->GetCommandLineArgument("-binary", isbinary);

  std::cout << std::endl;
  std::cout << "      parameters " << std::endl;
  std::cout << "           sigma " << sigma << std::endl;
  std::cout << "      relaxation " << relaxation << std::endl;
  std::cout << "      iterations " << numberOfIterations << std::endl;
  std::cout << "number of points " << numberOfPoints << std::endl;
  std::cout << " image is binary " << isbinary << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image
  FloatImageType::Pointer image = FloatImageType::New();
  if (!readImage<FloatImageType>(image, imageFile)) {
    return EXIT_FAILURE;
  }

  // compute the minimum and the maximum intensity values of label
  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
  labelValues->SetImage(image);
  labelValues->Compute();

  double levelValue;
  if (parser->ArgumentExists("-level")) {
    parser->GetCommandLineArgument("-level", levelValue);
  }
  else {
    levelValue = 0.5*(labelValues->GetMinimum() + labelValues->GetMaximum());
  }

  std::cout << "input image " << imageFile << std::endl;
  std::cout << "       size " << image->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << image->GetSpacing() << std::endl;
  std::cout << "     origin " << image->GetOrigin() << std::endl;
  std::cout << "level value " << levelValue << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // image processing

  // smoothing
  typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(image);
  gaussian->SetSigma(sigma);

  // fill holes after smoothing
  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(gaussian->GetOutput());
  fillholes->SetFullyConnected(true);
  try {
    fillholes->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(fillholes->GetOutput());
  convertor->Update();

  typedef vtkSmartPointer<vtkMarchingCubes> MarchingCubes;
  MarchingCubes mcubes = MarchingCubes::New();
  mcubes->SetInputData(convertor->GetOutput());
  mcubes->SetValue(0, levelValue);
  try {
    mcubes->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  vtkSmartPointer<vtkPolyData> surface = mcubes->GetOutput();

  // decimate surface
  if (numberOfPoints > 0) {
    double reduction = 1 - numberOfPoints / (double)surface->GetNumberOfPoints();
    std::cout << "reduction to decimate surface " << reduction << std::endl;
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(surface);
    decimate->SetTargetReduction(reduction);
    decimate->SetPreserveTopology(true);
    decimate->SetSplitting(false);
    try {
      decimate->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    surface = decimate->GetOutput();
  }

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(surface);
  smoother->SetNumberOfIterations(numberOfIterations);
  smoother->SetRelaxationFactor(relaxation);
  try {
    smoother->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  typedef vtkSmartPointer<vtkPolyDataNormals> PolyDataNormals;
  PolyDataNormals normals = PolyDataNormals::New();
  normals->SetInputData(smoother->GetOutput());
  normals->AutoOrientNormalsOn();
  normals->FlipNormalsOff();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->SplittingOff();
  try {
    normals->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  surface = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  typedef std::pair<std::string, std::string> PairType;
  std::vector<PairType> surfaceInfo;
  surfaceInfo.push_back(PairType("number of points", std::to_string(surface->GetNumberOfPoints())));
  surfaceInfo.push_back(PairType("number of cells", std::to_string(surface->GetNumberOfCells())));
  surfaceInfo.push_back(PairType("length of edges", std::to_string(AverageLengthOfEdges(surface))));
  surfaceInfo.push_back(PairType("area of cells", std::to_string(AverageAreaOfCells(surface))));

  std::cout << "output surface polydata " << surfaceFile << std::endl;
  std::cout << "       number of cells  " << surface->GetNumberOfCells() << std::endl;
  std::cout << "       number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << "average length of edges " << AverageLengthOfEdges(surface) << endl;
  std::cout << "  average area of cells " << AverageAreaOfCells(surface) << endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // setup point set
  typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
    PointSetType::PointType point;
    point.CastFrom<double>(surface->GetPoints()->GetPoint(n));
    pointSet->SetPoint(n, point);
  }

  // compute metrics
  typedef PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetInfo(surfaceInfo);
  if (isbinary) {
    metrics->SetMovingImage(ComputeDistanceMapImage(image, labelValues->GetMinimum(), labelValues->GetMaximum()));
  }
  else {
    metrics->SetMovingImage(image);
  }
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (parser->ArgumentExists("-report")) {
    std::string fileName;
    parser->GetCommandLineArgument("-report", fileName);
    std::cout << "print report to the file: " << fileName << std::endl;
    metrics->PrintReportToFile(fileName, getBaseNameFromPath(surfaceFile));
  }

  return EXIT_SUCCESS;
}


double AverageLengthOfEdges(vtkPolyData*poly)
{
  const unsigned int numberOfCells = poly->GetNumberOfCells();
  double sum = 0;

  for (int n = 0; n < poly->GetNumberOfCells(); ++n) {
    double p1[3], p2[3], p3[3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(n);
    vtkSmartPointer<vtkPoints> points = cell->GetPoints();

    points->GetPoint(0, p1);
    points->GetPoint(1, p2);
    points->GetPoint(2, p3);

    sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + 
           sqrt(vtkMath::Distance2BetweenPoints(p1, p3)) +
           sqrt(vtkMath::Distance2BetweenPoints(p2, p3));
  }

  return sum / (3*numberOfCells);
}


double AverageAreaOfCells(vtkPolyData*poly)
{
  const unsigned int numberOfCells = poly->GetNumberOfCells();
  double sum = 0;

  for (int n = 0; n < numberOfCells; ++n) {
    double a[3], b[3], c[3];
    double p1[3], p2[3], p3[3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(n);
    vtkSmartPointer<vtkPoints> points = cell->GetPoints();

    points->GetPoint(0, p1);
    points->GetPoint(1, p2);
    points->GetPoint(2, p3);
    
    vtkMath::Subtract(p1, p2, a);
    vtkMath::Subtract(p1, p3, b);
    vtkMath::Cross(a, b, c);

    sum += sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / 2;
  }

  return sum / numberOfCells;
}


FloatImageType::Pointer ComputeDistanceMapImage(FloatImageType::Pointer image, float background, float foreground)
{
  typedef itk::SignedMaurerDistanceMapImageFilter<FloatImageType, FloatImageType> DistanceFilterType;
  DistanceFilterType::Pointer distanceToForeground = DistanceFilterType::New();
  distanceToForeground->SetUseImageSpacing(true);
  distanceToForeground->SetInput(image);
  distanceToForeground->SetBackgroundValue(background);
  distanceToForeground->SetInsideIsPositive(false);

  DistanceFilterType::Pointer distanceToBackground = DistanceFilterType::New();
  distanceToBackground->SetUseImageSpacing(true);
  distanceToBackground->SetInput(image);
  distanceToBackground->SetBackgroundValue(foreground);
  distanceToBackground->SetInsideIsPositive(true);

  typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
  AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
  addfilter->SetInput1(distanceToForeground->GetOutput());
  addfilter->SetInput2(distanceToBackground->GetOutput());

  typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
  FilterType::Pointer multiply = FilterType::New();
  multiply->SetInput(addfilter->GetOutput());
  multiply->SetConstant(0.5);
  try {
    multiply->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  return multiply->GetOutput();
}

