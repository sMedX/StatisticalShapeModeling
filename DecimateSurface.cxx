#include <itkMesh.h>
#include <vtkPolyData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

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

  size_t numberOfPoints = 1e+05;
  parser->GetCommandLineArgument("-point", numberOfPoints);

  double relaxation = 0.2;
  parser->GetCommandLineArgument("-relaxation", relaxation);

  int iterations = 10;
  parser->GetCommandLineArgument("-iteration", iterations);

  std::cout << std::endl;
  std::cout << "parameters" << std::endl;
  std::cout << "     points " << numberOfPoints << std::endl;
  std::cout << " relaxation " << relaxation << std::endl;
  std::cout << " iterations " << iterations << std::endl;
  std::cout << std::endl;

  // read input polydata
  vtkSmartPointer<vtkPolyData> inputSurface = vtkSmartPointer<vtkPolyData>::New();
  readVTKPolydata(inputSurface, inputFile);

  std::cout << "output surface polydata info" << std::endl;
  std::cout << inputFile << std::endl;
  std::cout << " number of cells " << inputSurface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << inputSurface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  // decimate surface
  double reduction = 1 - numberOfPoints / (double)inputSurface->GetNumberOfPoints();
  std::cout << "reduction to decimate surface " << reduction << std::endl;

  vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
  decimate->SetInputData(inputSurface);
  decimate->SetTargetReduction(reduction);
  decimate->SetPreserveTopology(true);
  decimate->SetSplitting(false);
  decimate->Update();

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(decimate->GetOutput());
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

  vtkSmartPointer<vtkPolyData> outputSurface = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(outputSurface, outputFile)) {
    return EXIT_FAILURE;
  }
  
  std::cout << "output surface polydata info" << std::endl;
  std::cout << outputFile << std::endl;
  std::cout << " number of cells " << outputSurface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << outputSurface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


