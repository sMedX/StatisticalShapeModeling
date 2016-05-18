#include <itkMesh.h>
#include <itkResampleImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <vtkPolyData.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <itkImageToVTKImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
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

  std::string outputSurfaceFile;
  parser->GetCommandLineArgument("-surface", outputSurfaceFile);

  std::string outputPotentialFile;
  parser->GetCommandLineArgument("-potential", outputPotentialFile);

  float spacing = 1;
  parser->GetCommandLineArgument("-spacing", spacing);

  float sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int iterations = 10;
  parser->GetCommandLineArgument("-iteration", iterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << "    spacing " << spacing << std::endl;
  std::cout << "      sigma " << sigma << std::endl;
  std::cout << " relaxation " << relaxation << std::endl;
  std::cout << " iterations " << iterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image

  BinaryImageType::Pointer label = BinaryImageType::New();
  if (!readImage<BinaryImageType>(label, labelFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input label info" << std::endl;
  std::cout << labelFile << std::endl;
  std::cout << "   size " << label->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << label->GetSpacing() << std::endl;
  std::cout << " origin " << label->GetOrigin() << std::endl;
  std::cout << std::endl;

  typedef itk::RecursiveGaussianImageFilter<BinaryImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(label);
  gaussian->SetSigma(sigma);
  gaussian->Update();

  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(gaussian->GetOutput());
  fillholes->SetFullyConnected(true);
  fillholes->Update();

  // resampling of input image
  BinaryImageType::SpacingType outputSpacing;
  outputSpacing.Fill(spacing);

  BinaryImageType::SizeType outputSize;

  for (int i = 0; i < Dimension; ++i) {
    outputSize[i] = label->GetLargestPossibleRegion().GetSize()[i] * label->GetSpacing()[i] / spacing;
  }

  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetDefaultPixelValue(0);
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(fillholes->GetOutput()->GetOrigin());
  resampler->SetOutputDirection(fillholes->GetOutput()->GetDirection());
  resampler->SetInput(fillholes->GetOutput());
  resampler->Update();

  FloatImageType::Pointer processedImage = resampler->GetOutput();

  //----------------------------------------------------------------------------
  //compute surface
  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
  calculator->SetImage(processedImage);
  calculator->Compute();
  float levelValue = 0.5*(calculator->GetMinimum() + calculator->GetMaximum());

  //convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(processedImage);
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

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(mcubes->GetOutput());
  smoother->SetNumberOfIterations(iterations);
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
  if (!writeVTKPolydata(surface, outputSurfaceFile)) {
    return EXIT_FAILURE;
  }
  
  std::cout << "output surface polydata info" << std::endl;
  std::cout << outputSurfaceFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  //compute potential map image
  if (!parser->ArgumentExists("-potential")) {
    return EXIT_SUCCESS;
  }
  typedef itk::BinaryThresholdImageFilter <FloatImageType, BinaryImageType> BinaryThresholdImageFilterType;
  BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
  threshold->SetInput(processedImage);
  threshold->SetLowerThreshold(levelValue);
  threshold->SetUpperThreshold(std::numeric_limits<float>::max());
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);

  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, FloatImageType> SignedMaurerDistanceMapImageFilterType;
  SignedMaurerDistanceMapImageFilterType::Pointer distancemap = SignedMaurerDistanceMapImageFilterType::New();
  distancemap->SetUseImageSpacing(true);
  distancemap->SetInput(threshold->GetOutput());
  distancemap->Update();

  FloatImageType::Pointer potentialImage = distancemap->GetOutput();

  // write potential image to the file
  if (!writeImage<FloatImageType>(potentialImage, outputPotentialFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "output potential image info" << std::endl;
  std::cout << outputPotentialFile << std::endl;
  std::cout << "   size " << potentialImage->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << potentialImage->GetSpacing() << std::endl;
  std::cout << " origin " << potentialImage->GetOrigin() << std::endl;

  return EXIT_SUCCESS;
}


