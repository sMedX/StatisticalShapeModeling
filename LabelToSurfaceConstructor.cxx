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

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string labelFile;
  parser->GetCommandLineArgument("-label", labelFile);

  std::string surfaceFile;
  parser->GetCommandLineArgument("-surface", surfaceFile);

  size_t numberOfPoints = 1e+05;
  parser->GetCommandLineArgument("-point", numberOfPoints);

  int radius = 5;
  parser->GetCommandLineArgument("-radius", radius);

  float spacing = 1;
  parser->GetCommandLineArgument("-spacing", spacing);

  float sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << "parameters " << std::endl;
  std::cout << "    points " << numberOfPoints << std::endl;
  std::cout << "    radius " << radius << std::endl;
  std::cout << "   spacing " << spacing << std::endl;
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
  MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
  calculator->SetImage(label);
  calculator->Compute();
  float levelValue = 0.5*(calculator->GetMinimum() + calculator->GetMaximum());

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
  resampler->SetDefaultPixelValue(calculator->GetMinimum());
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

  // decimate surface
  double reduction = 1 - numberOfPoints / (double) mcubes->GetOutput()->GetNumberOfPoints();
  std::cout << "reduction to decimate surface " << reduction << std::endl;

  vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
  decimate->SetInputData(mcubes->GetOutput());
  decimate->SetTargetReduction(reduction);
  decimate->SetPreserveTopology(true);
  decimate->SetSplitting(false);
  decimate->Update();

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(decimate->GetOutput());
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
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


