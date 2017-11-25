#include <itkDataManager.h>
#include <itkDirectory.h>
#include <itkImage.h>
#include <itkLandmarkBasedTransformInitializer.h>
#include <itkMesh.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include <itkPCAModelBuilder.h>
#include <itkRigid3DTransform.h>
#include <itkStandardMeshRepresenter.h>
#include <itkStatisticalModel.h>
#include <itkTransformMeshFilter.h>
#include "utils/statismo-build-models-utils.h"

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmModelBuildingOptions.h"

void buildAndSaveShapeModel(const StringVector & listOfFiles, const ssm::ModelBuildingOptions & options);

int main(int argc, char** argv)
{
  // read options from config file
  ssm::ModelBuildingOptions options;
  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  if (!options.ParseConfigFile()) {
    return EXIT_FAILURE;
  }
  options.PrintConfig();

  //----------------------------------------------------------------------------
  // read list of files
  StringVector listOfFiles;
  try {
    listOfFiles = readListFromFile(options.GetInputFileName());
  }
  catch (std::ifstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // build shape model
  try {
    buildAndSaveShapeModel(listOfFiles, options);
  }
  catch (itk::ExceptionObject & e) {
    std::cerr << "Could not build the model:" << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

void buildAndSaveShapeModel(const StringVector & listOfFiles, const ssm::ModelBuildingOptions & options)
{
  const unsigned Dimensions = 3;

  typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
  RepresenterType::Pointer representer = RepresenterType::New();

  typedef itk::Mesh<float, Dimensions> MeshType;
  typedef itk::DataManager<MeshType> DataManagerType;
  DataManagerType::Pointer dataManager = DataManagerType::New();

  typedef itk::MeshFileReader<MeshType> MeshReaderType;
  typedef std::vector<MeshReaderType::Pointer> MeshReaderVector;
  MeshReaderVector meshes;

  for (const auto & fileName : listOfFiles) {
    MeshReaderType::Pointer reader = MeshReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();
    meshes.push_back(reader);
  }

  if (meshes.size() == 0) {
    itkGenericExceptionMacro(<< "The specified list of surfaces is empty.");
  }

  if (options.GetAlignmentMode() == "reference") {
    typedef itk::MeshFileReader<MeshType> MeshReaderType;
    MeshReaderType::Pointer reader = MeshReaderType::New();
    reader->SetFileName(options.GetReferenceFileName());
    reader->Update();
    representer->SetReference(reader->GetOutput());
  }
  else {
    std::vector<MeshType::Pointer> originalMeshes;

    for (MeshReaderVector::iterator it = meshes.begin(); it != meshes.end(); ++it) {
      MeshReaderType::Pointer reader = *it;
      originalMeshes.push_back(reader->GetOutput());
    }

    const size_t numberOfGPAIterations = 100;
    const size_t numberOfPoints = 1000;
    const double breakIfChangeBelow = 0.001;

    typedef itk::VersorRigid3DTransform< float > Rigid3DTransformType;
    typedef itk::Image<float, Dimensions> ImageType;
    typedef itk::LandmarkBasedTransformInitializer<Rigid3DTransformType, ImageType, ImageType> LandmarkBasedTransformInitializerType;
    typedef itk::TransformMeshFilter< MeshType, MeshType, Rigid3DTransformType > FilterType;
    MeshType::Pointer reference = calculateProcrustesMeanMesh<MeshType, LandmarkBasedTransformInitializerType, Rigid3DTransformType, FilterType>(originalMeshes, numberOfGPAIterations, numberOfPoints, breakIfChangeBelow);
    representer->SetReference(reference);
  }

  dataManager->SetRepresenter(representer);

  for (MeshReaderVector::const_iterator it = meshes.begin(); it != meshes.end(); ++it) {
    MeshReaderType::Pointer reader = *it;
    dataManager->AddDataset(reader->GetOutput(), reader->GetFileName());
  }

  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model;
  typedef itk::PCAModelBuilder<MeshType> PCAModelBuilder;
  PCAModelBuilder::Pointer pcaModelBuilder = PCAModelBuilder::New();
  model = pcaModelBuilder->BuildNewModel(dataManager->GetData(), options.GetNoise());

  std::cout << "shape model saved to the file " << options.GetOutputFileName() << std::endl;
  std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
  std::cout << "number of cells      " << model->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
  std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;

  model->Save(options.GetOutputFileName().c_str());
}
