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

#include "SurfaceToImageRegistrationFilter.h"
#include "SurfaceToLevelSetImageFilter.h"

#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  int numberOfIterations = 100;
  parser->GetCommandLineArgument("-iteration", numberOfIterations);

  //----------------------------------------------------------------------------
  // read list of files
  StringList fileNames = getFileList(listFile);

  std::string fileName = fileNames.begin()->c_str();
  MeshType::Pointer surface = MeshType::New();
  if (!readMesh<MeshType>(surface, fileName)) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // compute distance map image
  typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToPotentialImageFilterType;
  SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
  levelset->SetInput(surface);
  levelset->Update();

  //----------------------------------------------------------------------------
  //
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

    typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
    SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
    surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
    surfaceToImageRegistration->SetInput(surface);
    surfaceToImageRegistration->SetPotentialImage(levelset->GetOutput());
    surfaceToImageRegistration->Update();
    surfaceToImageRegistration->PrintReport(std::cout);
  }


    /*
  BinaryImageType::Pointer mask = BinaryImageType::New();
  if (!readImage<BinaryImageType>(mask, maskFile)) {
    return EXIT_FAILURE;
  }

  std::cout << "input mask info" << std::endl;
  std::cout << maskFile << std::endl;
  std::cout << "   size " << mask->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "spacing " << mask->GetSpacing() << std::endl;
  std::cout << " origin " << mask->GetOrigin() << std::endl;
  std::cout << std::endl;
  */

  return EXIT_SUCCESS;
}


