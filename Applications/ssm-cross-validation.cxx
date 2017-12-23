#include <boost/scoped_ptr.hpp>

#include <DataManager.h>
#include <itkPCAModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <itkReducedVarianceModelBuilder.h>
#include <itkStatisticalModel.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToPointSetMetrics.h"
#include "ssmModelCrossValidationOptions.h"

typedef statismo::DataManager<MeshType> DataManagerType;
typedef itk::PCAModelBuilder<MeshType> ModelBuilderType;
typedef itk::ReducedVarianceModelBuilder<MeshType> ReducedVarianceModelBuilderType;
typedef DataManagerType::CrossValidationFoldListType CVFoldListType;
typedef DataManagerType::DataItemListType DataItemListType;

int main(int argc, char** argv)
{
  // parse options
  ssm::ModelCrossValidationOptions options;

  if (!options.ParseOptions(argc, argv)) {
    return EXIT_FAILURE;
  }

  if (!checkFileName(options.GetReportFileName())) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // read list of files
  StringVector listOfFiles;
  try {
    listOfFiles = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  try {
    auto fileName = listOfFiles.begin()->c_str();
    auto surface = MeshType::New();
    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }

    // create a data manager and add a number of datasets for model building
    auto representer = RepresenterType::New();
    representer->SetReference(surface);

    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer));

    for (const auto &fileName : listOfFiles) {
      auto surface = MeshType::New();
      if (!readMesh<MeshType>(surface, fileName)) {
        return EXIT_FAILURE;
      }
      dataManager->AddDataset(surface, fileName);
    }

    CVFoldListType cvFoldList = dataManager->GetCrossValidationFolds(dataManager->GetNumberOfSamples(), true);
    std::cout << "number of surfaces " << dataManager->GetNumberOfSamples() << std::endl;
    std::cout << "number of folds    " << cvFoldList.size() << std::endl;
    std::cout << std::endl;

    size_t components = options.GetNumberOfComponents();
    double generalizationAbility = 0;

    for (const auto & it : cvFoldList) {
      // create the model
      auto modelBuilder = ModelBuilderType::New();
      auto model = modelBuilder->BuildNewModel(it.GetTrainingData(), 0);

      // reduce the number of components
      auto reducedVarianceModelBuilder = ReducedVarianceModelBuilderType::New();
      if (components >= 0 && components < model->GetNumberOfPrincipalComponents()) {
        model = reducedVarianceModelBuilder->BuildNewModelWithLeadingComponents(model, components);
      }

      std::cout << "number of samples    " << it.GetTrainingData().size() << std::endl;
      std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
      std::cout << std::endl;

      // Now we can iterate over the test data and do whatever validation we would like to do.
      const DataItemListType testSamplesList = it.GetTestingData();

      for (const auto & it : testSamplesList) {
        auto surfaceFileName = it->GetDatasetURI();
        auto testSample = it->GetSample();
        auto outputSample = model->DrawSample(model->ComputeCoefficientsForDataset(testSample));

        typedef std::pair<std::string, std::string> PairType;
        std::vector<PairType> info;
        info.push_back(PairType("components", std::to_string(model->GetNumberOfPrincipalComponents())));

        // compute metrics
        typedef itk::PointSet<MeshType::PointType, MeshType::PointDimension> PointSetType;
        auto pointSet1 = PointSetType::New();
        pointSet1->SetPoints(testSample->GetPoints());

        auto pointSet2 = PointSetType::New();
        pointSet2->SetPoints(outputSample->GetPoints());

        typedef ssm::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
        auto metrics = PointSetToPointSetMetricsType::New();
        metrics->SetFixedPointSet(pointSet1);
        metrics->SetMovingPointSet(pointSet2);
        metrics->SetInfo(info);
        metrics->Compute();
        metrics->PrintReport(std::cout);

        generalizationAbility += metrics->GetRMSEValue();

        // print report to *.csv file
        metrics->PrintReportToFile(options.GetReportFileName(), getBaseNameFromPath(surfaceFileName));

        // write samples
        if (options.GetWrite()) {
          auto suffix = "-cvtest";
          auto fileName = addFileNameSuffix(surfaceFileName, suffix);
          if (!writeMesh<MeshType>(outputSample, fileName)) {
            return EXIT_FAILURE;
          }
        }
      }
    }

    generalizationAbility = generalizationAbility / dataManager->GetNumberOfSamples();
    std::cout << "generalization ability " << generalizationAbility << std::endl;

    if (options.GetReportFileName().size() > 0) {
      std::ofstream file(options.GetReportFileName(), std::ofstream::out | std::ofstream::app);
      file << "generalization ability; " << generalizationAbility << std::endl;
      file.close();
    }
  }
  catch (statismo::StatisticalModelException& e) {
    std::cerr << "Exception occurred while building the shape model" << std::endl;
    std::cerr << e.what() << std::endl;
  }
}
