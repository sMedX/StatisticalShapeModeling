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
#include "vtkMath.h"
#include "vtkCell.h"

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"



const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;


double EdgeAvarageLength(vtkPolyData*poly)
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


int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string labelFile;
  parser->GetCommandLineArgument("-label", labelFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  float spacing = 1;
  parser->GetCommandLineArgument("-spacing", spacing);

  int radius = 5;
  parser->GetCommandLineArgument("-radius", radius);

  float sigma = 5;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << "parameters " << std::endl;
  std::cout << "   spacing " << spacing << std::endl;
  std::cout << "    radius " << radius << std::endl;
  std::cout << "     sigma " << sigma << std::endl;
  std::cout << "relaxation " << relaxation << std::endl;
  std::cout << "iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image

  FloatImageType::Pointer label = FloatImageType::New();
  if (!readImage<FloatImageType>(label, labelFile)) {
    return EXIT_FAILURE;
  }

  // compute the minimum and the maximum intensity values of label
  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
  labelValues->SetImage(label);
  labelValues->Compute();
  float levelValue = 0.5*(labelValues->GetMinimum() + labelValues->GetMaximum());

  std::cout << "input label " << labelFile << std::endl;
  std::cout << "       size " << label->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << label->GetSpacing() << std::endl;
  std::cout << "     origin " << label->GetOrigin() << std::endl;
  std::cout << "level value " << levelValue << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // resampling of input image
  typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian1 = RecursiveGaussianImageFilterType::New();
  gaussian1->SetInput(label);
  gaussian1->SetSigma(sigma);

  BinaryImageType::SizeType outSize;
  BinaryImageType::SpacingType outSpacing;
  outSpacing.Fill(spacing);

  for (int i = 0; i < Dimension; ++i) {
    outSize[i] = label->GetLargestPossibleRegion().GetSize()[i] * label->GetSpacing()[i] / outSpacing[i];
  }

  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetInput(gaussian1->GetOutput());
  resampler->SetDefaultPixelValue(labelValues->GetMinimum());
  resampler->SetSize(outSize);
  resampler->SetOutputSpacing(outSpacing);
  resampler->SetOutputOrigin(label->GetOrigin());
  resampler->SetOutputDirection(label->GetDirection());

  typedef itk::BinaryThresholdImageFilter<FloatImageType, FloatImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer threshold = ThresholdImageFilterType::New();
  threshold->SetInput(resampler->GetOutput());
  threshold->SetLowerThreshold(levelValue);
  threshold->SetUpperThreshold(std::numeric_limits<FloatImageType::PixelType>::max());
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);

  BinaryImageType::SizeType nhoodRadius;
  nhoodRadius.Fill(radius);

  typedef itk::VotingBinaryHoleFillingImageFilter<FloatImageType, FloatImageType> VotingBinaryIterativeHoleFillingImageFilterType;
  VotingBinaryIterativeHoleFillingImageFilterType::Pointer holefilling = VotingBinaryIterativeHoleFillingImageFilterType::New();
  holefilling->SetInput(threshold->GetOutput());
  holefilling->SetBackgroundValue(0);
  holefilling->SetForegroundValue(1);
  holefilling->SetRadius(nhoodRadius);

  RecursiveGaussianImageFilterType::Pointer gaussian2 = RecursiveGaussianImageFilterType::New();
  gaussian2->SetInput(holefilling->GetOutput());
  gaussian2->SetSigma(sigma);

  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(gaussian2->GetOutput());
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

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(mcubes->GetOutput());
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

  vtkSmartPointer<vtkPolyData> surface = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, surfaceFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "output surface polydata " << surfaceFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  cout << "edge avarage length " << EdgeAvarageLength(surface) << endl;
  std::cout << std::endl;

  // compute metrics
  typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();

  for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
    pointSet->SetPoint(n, surface->GetPoints()->GetPoint(n));
  }

  // compute distance map
  typedef itk::SignedMaurerDistanceMapImageFilter<FloatImageType, FloatImageType> DistanceFilterType;
  DistanceFilterType::Pointer distanceToForeground = DistanceFilterType::New();
  distanceToForeground->SetUseImageSpacing(true);
  distanceToForeground->SetInput(label);
  distanceToForeground->SetBackgroundValue(labelValues->GetMinimum());
  distanceToForeground->SetInsideIsPositive(false);
  distanceToForeground->Update();

  DistanceFilterType::Pointer distanceToBackground = DistanceFilterType::New();
  distanceToBackground->SetUseImageSpacing(true);
  distanceToBackground->SetInput(label);
  distanceToBackground->SetBackgroundValue(labelValues->GetMaximum());
  distanceToBackground->SetInsideIsPositive(true);
  distanceToBackground->Update();

  typedef itk::AddImageFilter <FloatImageType> AddImageFilterType;
  AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
  addfilter->SetInput1(distanceToForeground->GetOutput());
  addfilter->SetInput2(distanceToBackground->GetOutput());
  addfilter->Update();

  typedef itk::MultiplyImageFilter <FloatImageType> FilterType;
  FilterType::Pointer multiply = FilterType::New();
  multiply->SetInput(addfilter->GetOutput());
  multiply->SetConstant(0.5);
  multiply->Update();

  FloatImageType::Pointer distancemap = multiply->GetOutput();

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

    std::cout << "write report to the file: " << fileName << std::endl;

    std::string dlm = ";";

    std::string header = dlm;
    std::string scores = getBaseNameFromPath(surfaceFile) + dlm;

    header += "Mean" + dlm;
    scores += std::to_string(metrics->GetMeanValue()) + dlm;

    header += "RMSE" + dlm;
    scores += std::to_string(metrics->GetRMSEValue()) + dlm;

    header += "Quantile " + std::to_string(metrics->GetLevelOfQuantile()) + dlm;
    scores += std::to_string(metrics->GetQuantileValue()) + dlm;

    header += "Maximal" + dlm;
    scores += std::to_string(metrics->GetMaximalValue()) + dlm;

    header += dlm;
    scores += dlm;

    header += "Number of points" + dlm;
    scores += std::to_string(surface->GetNumberOfPoints()) + dlm;

    header += "Number of cells" + dlm;
    scores += std::to_string(surface->GetNumberOfCells()) + dlm;

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


