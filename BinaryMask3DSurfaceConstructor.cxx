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

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string imageFile;
  parser->GetCommandLineArgument("-image", imageFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  float spacing = 2;
  parser->GetCommandLineArgument("-spacing", spacing);

  float sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxationFactor = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxationFactor);

  int numberOfIterations = 10;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << "input parameters" << std::endl;
  std::cout << "    spacing " << spacing << std::endl;
  std::cout << "      sigma " << sigma << std::endl;
  std::cout << " relaxation " << relaxationFactor << std::endl;
  std::cout << " iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image

  BinaryImageType::Pointer image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, imageFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input image " << imageFile << std::endl;
  std::cout << "   size " << image->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << image->GetSpacing() << std::endl;
  std::cout << " origin " << image->GetOrigin() << std::endl;
  std::cout << std::endl;

  // resampling of input image
  BinaryImageType::SpacingType outputSpacing;
  outputSpacing.Fill(spacing);

  BinaryImageType::SizeType outputSize;

  for (int i = 0; i < Dimension; ++i) {
    outputSize[i] = image->GetLargestPossibleRegion().GetSize()[i] * image->GetSpacing()[i] / spacing;
  }
  
  typedef itk::RecursiveGaussianImageFilter<BinaryImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetSigma(sigma);
  gaussian->SetInput(image);
  gaussian->Update();

  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetDefaultPixelValue(0);
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(gaussian->GetOutput()->GetOrigin());
  resampler->SetOutputDirection(gaussian->GetOutput()->GetDirection());
  resampler->SetInput(gaussian->GetOutput());
  resampler->Update();

  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(resampler->GetOutput());
  fillholes->SetFullyConnected(true);
  fillholes->Update();

  FloatImageType::Pointer processedImage = fillholes->GetOutput();

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
  smoother->SetNumberOfIterations(numberOfIterations);
  smoother->SetRelaxationFactor(relaxationFactor);

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

  vtkSmartPointer<vtkPolyData> polyData = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(polyData, outputFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "output polydata " << outputFile << std::endl;
  std::cout << " number of cells " << polyData->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << polyData->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


