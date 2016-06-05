#include <itkMesh.h>
#include <itkImage.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <itkImageToVTKImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "SurfaceToLevelSetImageFilter.h"
#include "SurfaceToImageRegistrationFilter.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string outputLevelset;
  parser->GetCommandLineArgument("-levelset", outputLevelset);

  std::string outputSurface;
  parser->GetCommandLineArgument("-surface", outputSurface);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  std::cout << std::endl;
  std::cout << "reference shape constructor" << std::endl;
  std::cout << "   list of files " << listFile << std::endl;
  std::cout << " output levelset " << outputLevelset << std::endl;
  std::cout << "  output surface " << outputSurface << std::endl;
  std::cout << "      iterations " << numberOfIterations << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read list of files
  StringList fileNames = getFileList(listFile);

  std::string fileName = fileNames.begin()->c_str();
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, fileName)) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // compute level set image
  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToPotentialImageFilterType;
  SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
  levelset->SetMargin(0.5);
  levelset->SetInput(surface);

  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  FloatImageType::Pointer reference = levelset->GetOutput();

  // reference image data
  FloatImageType::PointType origin = reference->GetOrigin();
  FloatImageType::SpacingType spacing = reference->GetSpacing();
  FloatImageType::SizeType size = reference->GetLargestPossibleRegion().GetSize();

  // allocate image
  FloatImageType::Pointer updateReference = FloatImageType::New();
  updateReference->SetRegions(reference->GetLargestPossibleRegion());
  updateReference->Allocate();
  updateReference->FillBuffer(0);
  updateReference->CopyInformation(reference);

  //----------------------------------------------------------------------------
  // compute reference image
  double count = 0;

  for (StringList::const_iterator it = fileNames.begin(); it != fileNames.end(); ++it) {
    std::string fileName = it->c_str();
    MeshType::Pointer surface = MeshType::New();
    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }

    std::cout << std::endl;
    std::cout << "surface polydata info " << fileName << std::endl;
    std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;

    // perform registration
    typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
    SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
    surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
    surfaceToImageRegistration->SetInput(surface);
    surfaceToImageRegistration->SetPotentialImage(reference);

    try {
      surfaceToImageRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    surfaceToImageRegistration->PrintReport(std::cout);

    // compute level set image
    typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToPotentialImageFilterType;
    SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
    levelset->SetOrigin(origin);
    levelset->SetSpacing(spacing);
    levelset->SetSize(size);
    levelset->SetInput(surfaceToImageRegistration->GetOutput());

    try {
      levelset->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    // add current levelset to the reference image
    count++;
    double const1 = 1 - 1 / count;
    double const2 = 1 / count;

    typedef itk::ImageRegionIterator<FloatImageType> IteratorType;
    IteratorType it1(updateReference, updateReference->GetLargestPossibleRegion());
    IteratorType it2(levelset->GetOutput(), levelset->GetOutput()->GetLargestPossibleRegion());

    for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2) {
      double pixelValue = const1 * it1.Get() + const2*it2.Get();
      it1.Set(pixelValue);
    }
  }

  std::cout << "output reference image " << outputLevelset << std::endl;
  std::cout << "   size " << updateReference->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << updateReference->GetSpacing() << std::endl;
  std::cout << " origin " << updateReference->GetOrigin() << std::endl;
  std::cout << std::endl;

  if (!writeImage<FloatImageType>(updateReference, outputLevelset)) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // compute reference surface

  //convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(updateReference);
  convertor->Update();

  double levelValue = 0;
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
  float relaxation = 0.2;
  int iterations = 10;

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

  // write surface to the file
  std::cout << "output reference surface " << outputSurface << std::endl;
  std::cout << " number of cells " << normals->GetOutput()->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << normals->GetOutput()->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  if (!writeVTKPolydata(normals->GetOutput(), outputSurface)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


