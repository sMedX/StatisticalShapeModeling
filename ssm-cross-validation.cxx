﻿#include <boost/program_options.hpp>
#include <boost/scoped_ptr.hpp>

#include <utils/statismo-build-models-utils.h>
#include <DataManager.h>
#include <itkPCAModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <itkReducedVarianceModelBuilder.h>
#include <itkStatisticalModel.h>
#include <itkPointsLocator.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "ssm/ssmPointSetToPointSetMetrics.h"

typedef MeshType::PointType PointType;
typedef statismo::DataManager<MeshType> DataManagerType;
typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;

typedef itk::PCAModelBuilder<MeshType> ModelBuilderType;
typedef itk::ReducedVarianceModelBuilder<MeshType> ReducedVarianceModelBuilderType;
typedef statismo::VectorType VectorType;
typedef DataManagerType::CrossValidationFoldListType CVFoldListType;
typedef DataManagerType::DataItemListType DataItemListType;

struct ProgramOptions
{
  bool help;
  bool write;
  std::string listFile;
  std::string reportFile;
  std::vector<size_t> components;
  size_t numberOfSamples = 5;
};

typedef boost::filesystem::path fp;
namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);

double computeSpecificity(const DataItemListType& testMeshes, unsigned numberOfSamples);

int main(int argc, char** argv)
{
  ProgramOptions options;
  po::options_description description = initializeProgramOptions(options);
  po::variables_map vm;
  try {
    po::parsed_options parsedOptions = po::command_line_parser(argc, argv).options(description).run();
    po::store(parsedOptions, vm);
    po::notify(vm);
  }
  catch (po::error& e) {
    cerr << "An exception occurred while parsing the command line:" << endl;
    cerr << e.what() << endl << endl;
    cout << description << endl;
    return EXIT_FAILURE;
  }
  if (options.help == true) {
    cout << description << endl;
    return EXIT_SUCCESS;
  }

  StringList fileNames;
  try {
    fileNames = getFileList(options.listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "Could not read the data-list:" << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  try {
    std::string fileName = fileNames.begin()->c_str();
    MeshType::Pointer surface = MeshType::New();
    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }

    // create a data manager and add a number of datasets for model building
    RepresenterType::Pointer representer = RepresenterType::New();
    representer->SetReference(surface);

    boost::scoped_ptr<DataManagerType> dataManager(DataManagerType::Create(representer));

    for (StringList::const_iterator it = fileNames.begin(); it != fileNames.end(); ++it) {
      std::string fileName = it->c_str();
      MeshType::Pointer surface = MeshType::New();
      if (!readMesh<MeshType>(surface, fileName)) {
        return EXIT_FAILURE;
      }
      dataManager->AddDataset(surface, fileName);
    }

    CVFoldListType cvFoldList = dataManager->GetCrossValidationFolds(dataManager->GetNumberOfSamples(), true);
    std::cout << "number of surfaces " << dataManager->GetNumberOfSamples() << std::endl;
    std::cout << "number of folds    " << cvFoldList.size() << std::endl;
    std::cout << std::endl;

    if (options.components.size() == 0) {
      options.components.push_back(dataManager->GetNumberOfSamples());
    }

    double specificity = itk::NumericTraits<double>::Zero;
    double generalization = itk::NumericTraits<double>::Zero;

    // iterate over cvFoldList to get all the folds
    for (const auto & components : options.components) {
      for (CVFoldListType::const_iterator it = cvFoldList.begin(); it != cvFoldList.end(); ++it) {
        // create the model
        ModelBuilderType::Pointer pcaModelBuilder = ModelBuilderType::New();
        StatisticalModelType::Pointer model = pcaModelBuilder->BuildNewModel(it->GetTrainingData(), 0);

        // reduce the number of components
        ReducedVarianceModelBuilderType::Pointer reducedVarModelBuilder = ReducedVarianceModelBuilderType::New();
        if (components < model->GetNumberOfPrincipalComponents()) {
          model = reducedVarModelBuilder->BuildNewModelWithLeadingComponents(model, components);
        }

        std::cout << "built model from     " << it->GetTrainingData().size() << " samples" << std::endl;
        std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
        std::cout << std::endl;

        // Now we can iterate over the test data and do whatever validation we would like to do.
        const DataItemListType testSamplesList = it->GetTestingData();

        for (DataItemListType::const_iterator it = testSamplesList.begin(); it != testSamplesList.end(); ++it) {
          std::string surfaceName = (*it)->GetDatasetURI();

          MeshType::Pointer testSample = (*it)->GetSample();
          MeshType::Pointer outputSample = model->DrawSample(model->ComputeCoefficientsForDataset(testSample));

          typedef std::pair<std::string, std::string> PairType;
          std::vector<PairType> info;
          info.push_back(PairType("Components", std::to_string(model->GetNumberOfPrincipalComponents())));
          info.push_back(PairType("Probability", std::to_string(model->ComputeProbabilityOfDataset(testSample))));

          // compute metrics
          typedef itk::PointSet<MeshType::PointType, MeshType::PointDimension> PointSetType;
          PointSetType::Pointer pointSet1 = PointSetType::New();
          pointSet1->SetPoints(testSample->GetPoints());

          PointSetType::Pointer pointSet2 = PointSetType::New();
          pointSet2->SetPoints(outputSample->GetPoints());

          typedef ssm::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
          PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
          metrics->SetFixedPointSet(pointSet1);
          metrics->SetMovingPointSet(pointSet2);
          metrics->SetInfo(info);
          metrics->Compute();
          metrics->PrintReport(std::cout);

          generalization += metrics->GetRMSEValue();

          // print report to *.csv file
          if (options.reportFile != "") {
            std::cout << "print report to the file: " << options.reportFile << std::endl;
            metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(surfaceName));
          }

          // write samples
          if (options.write) {
            std::string suffix = "-cv-" + std::to_string(model->GetNumberOfPrincipalComponents());
            std::string fileName = addFileNameSuffix(surfaceName, suffix);
            if (!writeMesh<MeshType>(outputSample, fileName)) {
              return EXIT_FAILURE;
            }
          }
        }
      }
    }

    DataItemListType allMeshes;
    DataItemListType trainData = cvFoldList.begin()->GetTrainingData();
    DataItemListType testData = cvFoldList.begin()->GetTestingData();
    allMeshes.insert(allMeshes.end(), trainData.begin(), trainData.end());
    allMeshes.insert(allMeshes.end(), testData.begin(), testData.end());

    specificity = computeSpecificity(allMeshes, options.numberOfSamples);
    generalization = generalization / dataManager->GetNumberOfSamples();

    std::cout << "specificity:      " << specificity << std::endl;
    std::cout << "generalization:   " << generalization << std::endl;

    if (options.reportFile != "") {
      std::ofstream file(options.reportFile, std::ofstream::out | std::ofstream::app);

      file << "Specificity: " << specificity << std::endl;
      file << "Generalization: " << generalization << std::endl;
      file.close();
    }
  }
  catch (statismo::StatisticalModelException& e) {
    std::cout << "Exception occurred while building the shape model" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

double computeEuclideanPointDist(MeshType::PointType pt1, MeshType::PointType pt2) {
  double sumSquares = 0;
  for (unsigned i = 0; i < pt1.Dimension; i++)  {
    double dist = pt1[i] - pt2[i];
    sumSquares += dist * dist;
  }
  return std::sqrt(sumSquares);
}

std::vector<double> computeDistance(MeshType::Pointer testMesh, MeshType::Pointer sample)
{
  std::vector<double> distanceValues;
  typedef itk::PointsLocator< MeshType::PointsContainer > PointsLocatorType;

  PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
  ptLocator->SetPoints(sample->GetPoints());
  ptLocator->Initialize();

  for (unsigned i = 0; i < testMesh->GetNumberOfPoints(); i++) {
    MeshType::PointType sourcePt = testMesh->GetPoint(i);
    int closestPointId = ptLocator->FindClosestPoint(sourcePt);
    MeshType::PointType targetPt = sample->GetPoint(closestPointId);
    distanceValues.push_back(computeEuclideanPointDist(sourcePt, targetPt));
  }
  return distanceValues;
}

double computeSpecificity(const DataItemListType& testMeshes, unsigned numberOfSamples)
{
  double accumulatedDistToClosestTrainingShape = 0;

  ModelBuilderType::Pointer pcaModelBuilder = ModelBuilderType::New();
  StatisticalModelType::Pointer model = pcaModelBuilder->BuildNewModel(testMeshes, 0);

  for (unsigned i = 0; i < numberOfSamples; i++) {
    MeshType::Pointer sample = model->DrawSample();
    double minDist = std::numeric_limits<double>::max();

    for (DataItemListType::const_iterator it = testMeshes.begin(); it != testMeshes.end(); ++it) {
      MeshType::Pointer testMesh = (*it)->GetSample();
      std::vector<double> distanceValues = computeDistance(testMesh, sample);
      double dist = std::accumulate(distanceValues.begin(), distanceValues.end(), 0.0) / testMesh->GetNumberOfPoints();

      if (dist < minDist) {
        minDist = dist;
      }
    }
    accumulatedDistToClosestTrainingShape += minDist;
  }
  double specificity = accumulatedDistToClosestTrainingShape / numberOfSamples;

  return specificity;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("list", po::value<std::string>(&options.listFile), "The path to the file with list of surfaces.")
    ;

  po::options_description input("Optional input options");
  input.add_options()
    ("write", po::value<bool>(&options.write), "Write surfaces.")
    ("components", po::value<std::vector<size_t>>(&options.components)->multitoken(), "The number of components for GP shape model.")
    ("samples", po::value<size_t>(&options.numberOfSamples)->default_value(options.numberOfSamples), "The number of samples for specificity calculation.")
    ;

  po::options_description report("Optional report options");
  report.add_options()
    ("report,r", po::value<std::string>(&options.reportFile), "The path to the file to print report.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(input).add(report).add(help);

  return description;
}
