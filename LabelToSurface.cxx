#include <boost/filesystem.hpp>

#include <itkMesh.h>
#include <vtkPolyData.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <itkImageToVTKImageFilter.h>
#include <itkVotingBinaryHoleFillingImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkBinaryMorphologicalClosingImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;


double EdgeAverageLength(vtkPolyData*poly);
double SquareAverage(vtkPolyData*poly);
FloatImageType::Pointer ComputeDistanceMapImage(FloatImageType::Pointer image, float background, float foreground);

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string imageFile;
  parser->GetCommandLineArgument("-label", imageFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  float spacing = 1;
  parser->GetCommandLineArgument("-spacing", spacing);

  int radius = 0;
  parser->GetCommandLineArgument("-radius", radius);

  int repair = 0;
  parser->GetCommandLineArgument("-repair", repair);

  float sigma = 0;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  int numberOfPoints = 0;
  parser->GetCommandLineArgument("-points", numberOfPoints);

  std::cout << std::endl;
  std::cout << "      parameters " << std::endl;
  std::cout << "         spacing " << spacing << std::endl;
  std::cout << "          radius " << radius << std::endl;
  std::cout << "          repair " << repair << std::endl;
  std::cout << "           sigma " << sigma << std::endl;
  std::cout << "      relaxation " << relaxation << std::endl;
  std::cout << "      iterations " << numberOfIterations << std::endl;
  std::cout << "number of points " << numberOfPoints << std::endl;
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

  //--------------------------------------------------------
  //Padding to remove possible holes

  FloatImageType::SizeType paddingSize;
  paddingSize.Fill(7.0);


  typedef itk::ConstantPadImageFilter<FloatImageType, FloatImageType> PadImageFilter;
  PadImageFilter::Pointer padding = PadImageFilter::New();
  padding->SetInput(image);
  padding->SetPadBound(paddingSize);
  padding->SetConstant(labelValues->GetMinimum());
  padding->Update();

  image = padding->GetOutput();

  //----------------------------------------------------------------------------
  // resampling of input image
  if (sigma > 0) {
    typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
    RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
    gaussian->SetInput(image);
    gaussian->SetSigma(sigma);
    gaussian->Update();
    image = gaussian->GetOutput();
  }

  BinaryImageType::SizeType outSize;
  BinaryImageType::SpacingType outSpacing;
  outSpacing.Fill(spacing);

  for (int i = 0; i < Dimension; ++i) {
    outSize[i] = image->GetLargestPossibleRegion().GetSize()[i] * image->GetSpacing()[i] / outSpacing[i];
  }


  typedef itk::NearestNeighborInterpolateImageFunction< FloatImageType, double >  InterpolatorType;

  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetInput(image);
  resampler->SetDefaultPixelValue(labelValues->GetMinimum());
  resampler->SetSize(outSize);
  resampler->SetOutputSpacing(outSpacing);
  resampler->SetOutputOrigin(image->GetOrigin());
  resampler->SetOutputDirection(image->GetDirection());
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resampler->SetInterpolator(interpolator);
  resampler->Update();

  //----------------------------------------------------------------------------
  // image processing
  typedef itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer threshold = ThresholdImageFilterType::New();
  threshold->SetInput(resampler->GetOutput());
  threshold->SetLowerThreshold(levelValue);
  threshold->SetUpperThreshold(std::numeric_limits<FloatImageType::PixelType>::max());
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);

  threshold->Update();

  image = threshold->GetOutput();
  if (repair > 0) {

    typedef itk::BinaryBallStructuringElement < FloatImageType::PixelType, FloatImageType::ImageDimension >
      StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(repair);
    structuringElement.CreateStructuringElement();

    typedef itk::BinaryMorphologicalClosingImageFilter < FloatImageType, FloatImageType, StructuringElementType >
      BinaryMorphologicalClosingImageFilterType;

    BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
    closingFilter->SetInput(threshold->GetOutput());
    closingFilter->SetKernel(structuringElement);
    closingFilter->SetForegroundValue(1.0);

    closingFilter->Update();
    image = closingFilter->GetOutput();
  }

  BinaryImageType::SizeType nhoodRadius;
  nhoodRadius.Fill(radius);

  typedef itk::VotingBinaryHoleFillingImageFilter<FloatImageType, FloatImageType> VotingBinaryIterativeHoleFillingImageFilterType;
  VotingBinaryIterativeHoleFillingImageFilterType::Pointer holefilling = VotingBinaryIterativeHoleFillingImageFilterType::New();
  holefilling->SetInput(image);
  holefilling->SetBackgroundValue(0);
  holefilling->SetForegroundValue(1);
  holefilling->SetRadius(nhoodRadius);
  holefilling->Update();
  image = holefilling->GetOutput();
  //image->Update();

  if (sigma > 0) {
    typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
    RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
    gaussian->SetInput(image);
    gaussian->SetSigma(sigma);
    image = gaussian->GetOutput();
  }

  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(image);
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
  mcubes->SetValue(0, 0.5);
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
  normals->ComputeCellNormalsOn();
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

  std::cout << "output surface polydata " << surfaceFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  cout << "edge average length " << EdgeAverageLength(surface) << endl;
  cout << "square average  " << SquareAverage(surface) << endl;
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

  // compute distance map
  FloatImageType::Pointer distancemap = image;

  bool isBinary = true;
  parser->GetCommandLineArgument("-binary", isBinary);

  if (isBinary) {
    distancemap = ComputeDistanceMapImage(image, labelValues->GetMinimum(), labelValues->GetMaximum());

  }

  // compute metrics
  typedef PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetMovingImage(distancemap);
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


double EdgeAverageLength(vtkPolyData*poly)
{



  double length = 0;
  const int N = poly->GetNumberOfCells();

  for (int c = 0; c < N; c++)
  {

    double pc[3][3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(c);

    vtkSmartPointer<vtkPoints> p = cell->GetPoints();
    p->GetPoint(0, pc[0]);
    p->GetPoint(1, pc[1]);
    p->GetPoint(2, pc[2]);
    length += sqrt(vtkMath::Distance2BetweenPoints(pc[0], pc[1]));
    length += sqrt(vtkMath::Distance2BetweenPoints(pc[0], pc[2]));
    length += sqrt(vtkMath::Distance2BetweenPoints(pc[1], pc[2]));




  }
  return length / (N * 3);

}
double SquareAverage(vtkPolyData*poly)
{


  int edgesN = 0;
  double length = 0;
  const int N = poly->GetNumberOfCells();


  for (int c = 0; c < N; c++)
  {

    double pc[3][3];
    double a[3], b[3], cv[3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(c);

    vtkSmartPointer<vtkPoints> p = cell->GetPoints();
    p->GetPoint(0, pc[0]);
    p->GetPoint(1, pc[1]);
    p->GetPoint(2, pc[2]);


    vtkMath::Subtract(pc[0], pc[1], a);
    vtkMath::Subtract(pc[0], pc[2], b);
    vtkMath::Cross(a, b, cv);


    length += sqrt(cv[0] * cv[0] + cv[1] * cv[1] + cv[2] * cv[2]) / 2;


  }
  return length / N;

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

