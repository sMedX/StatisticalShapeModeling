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

const unsigned Dimensions = 3;
typedef itk::Mesh<float, Dimensions> MeshType;
typedef itk::StatisticalModel<MeshType> ShapeModelType;

ShapeModelType::Pointer shapeModelBuilder(const StringVector & list, const ssm::ModelBuildingOptions & options);

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

  //----------------------------------------------------------------------------
  // read list of files
  StringVector list;
  try {
    list = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // build shape model
  try {
    auto model = shapeModelBuilder(list, options);

    std::cout << "shape model saved to the file " << options.GetOutputFileName() << std::endl;
    std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
    std::cout << "number of cells      " << model->GetRepresenter()->GetReference()->GetNumberOfCells() << std::endl;
    std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
    model->Save(options.GetOutputFileName().c_str());
  }
  catch (itk::ExceptionObject & e) {
    std::cerr << "Could not build or save the shape model." << std::endl;
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

ShapeModelType::Pointer shapeModelBuilder(const StringVector & list, const ssm::ModelBuildingOptions & options)
{
  typedef itk::StandardMeshRepresenter<float, Dimensions> RepresenterType;
  auto representer = RepresenterType::New();

  typedef itk::DataManager<MeshType> DataManagerType;
  auto dataManager = DataManagerType::New();

  typedef itk::MeshFileReader<MeshType> MeshReaderType;
  typedef std::vector<MeshReaderType::Pointer> MeshReaderVector;
  MeshReaderVector meshes;

  for (const auto & fileName : list) {
    auto reader = MeshReaderType::New();
    reader->SetFileName(fileName);
    reader->Update();
    meshes.push_back(reader);
  }

  if (meshes.size() == 0) {
    itkGenericExceptionMacro(<< "The specified list of surfaces is empty.");
  }

  if (options.GetMode() == "reference") {
    typedef itk::MeshFileReader<MeshType> MeshReaderType;
    auto reader = MeshReaderType::New();
    reader->SetFileName(options.GetReferenceFileName());
    reader->Update();
    representer->SetReference(reader->GetOutput());
  }
  else {
    std::vector<MeshType::Pointer> originalMeshes;
    for (const auto &reader : meshes) {
      originalMeshes.push_back(reader->GetOutput());
    }

    const size_t numberOfIterations = 100;
    const size_t numberOfPoints = 1000;
    const double breakIfChangeBelow = 0.001;

    typedef itk::VersorRigid3DTransform< float > Rigid3DTransformType;
    typedef itk::Image<float, Dimensions> ImageType;
    typedef itk::LandmarkBasedTransformInitializer<Rigid3DTransformType, ImageType, ImageType> LandmarkBasedTransformInitializerType;
    typedef itk::TransformMeshFilter< MeshType, MeshType, Rigid3DTransformType > FilterType;
    auto reference = calculateProcrustesMeanMesh<MeshType, LandmarkBasedTransformInitializerType, Rigid3DTransformType, FilterType>(originalMeshes, numberOfIterations, numberOfPoints, breakIfChangeBelow);
    representer->SetReference(reference);
  }

  dataManager->SetRepresenter(representer);

  for (const auto &reader : meshes) {
    dataManager->AddDataset(reader->GetOutput(), reader->GetFileName());
  }

  typedef itk::PCAModelBuilder<MeshType> ModelBuilderType;
  auto modelBuilder = ModelBuilderType::New();
  auto model = modelBuilder->BuildNewModel(dataManager->GetData(), options.GetNoise());

  return model;
}
