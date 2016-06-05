#include <itkMesh.h>
#include <itkImage.h>
#include <itkBoundingBox.h>
#include <itkAddImageFilter.h>
#include <itkImageDuplicator.h>
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

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

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

  //allocate image
  FloatImageType::Pointer updateReference = FloatImageType::New();
  updateReference->SetRegions(reference->GetLargestPossibleRegion());
  updateReference->Allocate();
  updateReference->FillBuffer(0);
  updateReference->CopyInformation(reference);

  //----------------------------------------------------------------------------
  // compute reference image
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

    //perform registration
    typedef SurfaceToImageRegistrationFilter<MeshType> SurfaceToImageRegistrationFilterType;
    SurfaceToImageRegistrationFilterType::Pointer surfaceToImageRegistration = SurfaceToImageRegistrationFilterType::New();
    surfaceToImageRegistration->SetNumberOfIterations(numberOfIterations);
    surfaceToImageRegistration->SetInput(surface);
    surfaceToImageRegistration->SetPotentialImage(reference);
    surfaceToImageRegistration->Update();
    surfaceToImageRegistration->PrintReport(std::cout);

    //compute level set image
    typedef SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToPotentialImageFilterType;
    SurfaceToPotentialImageFilterType::Pointer levelset = SurfaceToPotentialImageFilterType::New();
    levelset->SetInput(surfaceToImageRegistration->GetOutput());
    levelset->SetOrigin(reference->GetOrigin());
    levelset->SetSpacing(reference->GetSpacing());
    levelset->SetSize(reference->GetLargestPossibleRegion().GetSize());
    levelset->Update();

    //add images
    typedef itk::AddImageFilter <FloatImageType, FloatImageType> AddImageFilterType;
    AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
    addFilter->SetInput1(updateReference);
    addFilter->SetInput2(levelset->GetOutput());

    try {
      addFilter->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    updateReference = addFilter->GetOutput();
  }

  /*
  typedef itk::BinaryThresholdImageFilter <FloatImageType, BinaryImageType> BinaryThresholdImageFilterType;
  BinaryThresholdImageFilterType::Pointer threshold = BinaryThresholdImageFilterType::New();
  threshold->SetInput(baseImageUpdate);
  threshold->SetLowerThreshold(std::numeric_limits<FloatImageType::PixelType>::lowest());
  threshold->SetUpperThreshold(0);
  threshold->SetInsideValue(1);
  threshold->SetOutsideValue(0);
  threshold->Update();*/

  if (writeImage<FloatImageType>(updateReference, outputFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


