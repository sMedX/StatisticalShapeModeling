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

  int radius = 1;
  parser->GetCommandLineArgument("-radius", radius);

  float sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  int smoothingIterations = 10;
  parser->GetCommandLineArgument("-iteration", smoothingIterations);

  std::cout << " image file " << imageFile << std::endl;
  std::cout << "output file " << outputFile << std::endl;
  std::cout << "    spacing " << spacing << std::endl;
  std::cout << "     radius " << radius << std::endl;
  std::cout << "      sigma " << sigma << std::endl;
  std::cout << " iterations " << smoothingIterations << std::endl;

  //----------------------------------------------------------------------------
  // read image

  BinaryImageType::Pointer itkImage = BinaryImageType::New();
  if (!readImage<BinaryImageType>(itkImage, imageFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "information about the input image " << std::endl;
  std::cout << itkImage->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << itkImage->GetSpacing() << std::endl;

  // resampling of input image
  BinaryImageType::SpacingType outputSpacing;
  outputSpacing.Fill(spacing);

  BinaryImageType::SizeType outputSize;

  for (int i = 0; i < Dimension; ++i) {
    outputSize[i] = itkImage->GetLargestPossibleRegion().GetSize()[i] * itkImage->GetSpacing()[i] / spacing;
  }

  typedef itk::ResampleImageFilter<BinaryImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetDefaultPixelValue(0);
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(itkImage->GetOrigin());
  resampler->SetOutputDirection(itkImage->GetDirection());
  resampler->SetInput(itkImage);
  resampler->Update();

  typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(resampler->GetOutput());
  gaussian->SetSigma(sigma);
  gaussian->Update();

  // Create structuring element
  typedef itk::BinaryBallStructuringElement<FloatImageType::PixelType, FloatImageType::ImageDimension> StructuringElementType;  
  StructuringElementType  strel;
  strel.SetRadius(radius);
  strel.CreateStructuringElement();

  typedef itk::GrayscaleDilateImageFilter<FloatImageType, FloatImageType, StructuringElementType> GrayscaleDilateImageFilterType;
  GrayscaleDilateImageFilterType::Pointer dilate = GrayscaleDilateImageFilterType::New();
  dilate->SetInput(gaussian->GetOutput());
  dilate->SetKernel(strel);
  dilate->Update();

  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(dilate->GetOutput());
  fillholes->SetFullyConnected(true);
  fillholes->Update();

  typedef itk::GrayscaleErodeImageFilter<FloatImageType, FloatImageType, StructuringElementType> GrayscaleErodeImageFilterType;
  GrayscaleErodeImageFilterType::Pointer erode = GrayscaleErodeImageFilterType::New();
  erode->SetInput(fillholes->GetOutput());
  erode->SetKernel(strel);
  erode->Update();
  FloatImageType::Pointer processedImage = erode->GetOutput();

  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
  calculator->SetImage(erode->GetOutput());
  calculator->Compute();
  float levelSetValue = 0.5*(calculator->GetMinimum() + calculator->GetMaximum());

  //convert ITK image to VTK image
  //VtkImageDataPointer vtkImage = ItkVtkFloatImage3DConverter::itkToVtk(erode->GetOutput());
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(processedImage);
  convertor->Update();

  typedef vtkSmartPointer<vtkMarchingCubes> MarchingCubes;
  MarchingCubes mcubes = MarchingCubes::New();
  mcubes->SetInputData(convertor->GetOutput());
  mcubes->SetValue(0, levelSetValue);

  try {
    mcubes->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  float smoothingRelaxation = 0.2;
  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(mcubes->GetOutput());
  smoother->SetNumberOfIterations(smoothingIterations);
  smoother->SetRelaxationFactor(smoothingRelaxation);

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

  // write surface to file
  if (!writeVTKPolydata(polyData, outputFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "information about the output mesh " << std::endl;
  std::cout << " number of cells " << polyData->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << polyData->GetNumberOfPoints() << std::endl;

  return EXIT_SUCCESS;
}


