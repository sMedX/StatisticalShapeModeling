#include <limits>

#include <itkResampleImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkBinaryThresholdImageFilter.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string imageFile;
  parser->GetCommandLineArgument("-image", imageFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  float spacing = 1;
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
  std::cout << "      sigma " << sigma << std::endl;
  std::cout << "     radius " << radius << std::endl;

  //----------------------------------------------------------------------------
  // read image

  FloatImageType::Pointer image = FloatImageType::New();
  if (!readImage<FloatImageType>(image, imageFile)) {
    return EXIT_FAILURE;
  }

  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer calculator = MinimumMaximumImageCalculatorType::New();
  calculator->SetImage(image);
  calculator->Compute();
  double levelValue = 0.5*(calculator->GetMinimum() + calculator->GetMaximum());

  std::cout << "information about the input image " << std::endl;
  std::cout << "    image size " << image->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << " image spacing " << image->GetSpacing() << std::endl;
  std::cout << " minimum value " << calculator->GetMinimum() << std::endl;
  std::cout << " maximum value " << calculator->GetMaximum() << std::endl;

  // resampling of input image
  FloatImageType::SpacingType outputSpacing;
  outputSpacing.Fill(spacing);

  FloatImageType::SizeType outputSize;

  for (int i = 0; i < Dimension; ++i) {
    outputSize[i] = image->GetLargestPossibleRegion().GetSize()[i] * image->GetSpacing()[i] / spacing;
  }

  typedef itk::ResampleImageFilter<FloatImageType, FloatImageType> ResampleImageFilterType;
  ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
  resampler->SetDefaultPixelValue(0);
  resampler->SetSize(outputSize);
  resampler->SetOutputSpacing(outputSpacing);
  resampler->SetOutputOrigin(image->GetOrigin());
  resampler->SetOutputDirection(image->GetDirection());
  resampler->SetInput(image);
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

  typedef itk::BinaryThresholdImageFilter<FloatImageType, BinaryImageType> ThresholdImageFilterType;
  ThresholdImageFilterType::Pointer threshold = ThresholdImageFilterType::New();
  threshold->SetInput(erode->GetOutput());
  threshold->SetLowerThreshold(levelValue);
  threshold->SetUpperThreshold(std::numeric_limits<float>::max());
  threshold->SetOutsideValue(0);
  threshold->SetInsideValue(1);
  threshold->Update();

  if (!writeImage<BinaryImageType>(threshold->GetOutput(), outputFile)) {
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}


