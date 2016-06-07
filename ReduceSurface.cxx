#include <itkMesh.h>
#include <vtkPolyData.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <itkImageToVTKImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "SurfaceToLevelSetImageFilter.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string inputFile;
  parser->GetCommandLineArgument("-input", inputFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  float sp = 1;
  parser->GetCommandLineArgument("-spacing", sp);

  float sigma = 1;
  parser->GetCommandLineArgument("-sigma", sigma);

  float relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int iterations = 10;
  parser->GetCommandLineArgument("-iteration", iterations);

  std::cout << std::endl;
  std::cout << "input parameters" << std::endl;
  std::cout << "    spacing " << sp << std::endl;
  std::cout << "      sigma " << sigma << std::endl;
  std::cout << " relaxation " << relaxation << std::endl;
  std::cout << " iterations " << iterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read image

  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, inputFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "output surface polydata info" << std::endl;
  std::cout << inputFile << std::endl;
  std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  // compute level set image
  FloatImageType::SpacingType spacing;
  spacing.Fill(sp);

  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToPotentialImageFilterType;
  SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
  levelset->SetMargin(0.5);
  levelset->SetSpacing(spacing);
  levelset->SetInput(surface);

  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  //compute surface

  //convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(levelset->GetOutput());
  convertor->Update();

  float levelValue = sp/2;
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

  vtkSmartPointer<vtkPolyData> output = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(output, outputFile)) {
    return EXIT_FAILURE;
  }
  
  std::cout << "output surface polydata info" << std::endl;
  std::cout << inputFile << std::endl;
  std::cout << " number of cells " << output->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << output->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


