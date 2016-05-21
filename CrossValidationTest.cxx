#include <boost/scoped_ptr.hpp>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <DataManager.h>
#include <PCAModelBuilder.h>
#include <StatisticalModel.h>
#include <itkStandardMeshRepresenter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

// All the statismo classes have to be parameterized with the RepresenterType.
const unsigned int Dimension = 3;
typedef itk::Mesh<float, Dimension> MeshType;
typedef MeshType::PointType PointType;
typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
typedef statismo::DataManager<MeshType> DataManagerType;
typedef statismo::StatisticalModel<MeshType> StatisticalModelType;
typedef statismo::PCAModelBuilder<MeshType> ModelBuilderType;
typedef statismo::VectorType VectorType;
typedef DataManagerType::CrossValidationFoldListType CVFoldListType;
typedef DataManagerType::DataItemListType DataItemListType;

int main(int argc, char** argv)
{
  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();
  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  unsigned int folds = 0;
  parser->GetCommandLineArgument("-folds", folds);

  try {
    StringList fileNames = getFileList(listFile);

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

      // We provide the filename as a second argument. It will be written as metadata, and allows us to more easily figure out what we did later.
      dataManager->AddDataset(surface, fileName);
    }

    // create the model builder
    boost::scoped_ptr<ModelBuilderType> pcaModelBuilder(ModelBuilderType::Create());

    if (folds == 0) {
      folds = dataManager->GetNumberOfSamples();
    }

    CVFoldListType cvFoldList = dataManager->GetCrossValidationFolds(folds, true);

    std::cout << "successfully loaded " << dataManager->GetNumberOfSamples() << " samples " << std::endl;
    std::cout << "    number of folds " << cvFoldList.size() << std::endl;

    // iterate over cvFoldList to get all the folds
    for (CVFoldListType::const_iterator it = cvFoldList.begin(); it != cvFoldList.end(); ++it) {

      // build the model as usual
      boost::scoped_ptr<StatisticalModelType> model(pcaModelBuilder->BuildNewModel(it->GetTrainingData(), 0));
      std::cout << "built model with  " << model->GetNumberOfPrincipalComponents() << " principal components" << std::endl;

      // Now we can iterate over the test data and do whatever validation we would like to do.
      const DataItemListType testSamplesList = it->GetTestingData();

      for (DataItemListType::const_iterator it = testSamplesList.begin(); it != testSamplesList.end(); ++it) {
        std::string sampleName = (*it)->GetDatasetURI();

        MeshType::Pointer testSample = (*it)->GetSample();
        VectorType coefficients = model->ComputeCoefficientsForDataset(testSample);
        MeshType::Pointer outputSample = model->DrawSample(coefficients);

        double probability = model->ComputeProbabilityOfDataset(testSample);
        double mean = 0;
        double rmse = 0;
        double maximal = 0;
        PointType p1, p2;

        for (size_t n = 0; n < testSample->GetNumberOfPoints(); ++n) {
          p1.CastFrom<float>(testSample->GetPoint(n));
          p2.CastFrom<float>(outputSample->GetPoint(n));
          double dist = p1.EuclideanDistanceTo(p2);

          mean += dist;
          rmse += dist * dist;
          maximal = std::max(maximal, dist);
        }

        mean = mean / testSample->GetNumberOfPoints();
        rmse = std::sqrt(rmse / testSample->GetNumberOfPoints());

        std::cout << "Metric values for the sample " << sampleName << std::endl;
        std::cout << "Probability " << probability << std::endl;
        std::cout << "       Mean " << mean << " mm" << std::endl;
        std::cout << "       RMSE " << rmse << " mm" << std::endl;
        std::cout << "    Maximal " << maximal << " mm" << std::endl;

        if (parser->ArgumentExists("-report")) {
          std::string reportFile;
          parser->GetCommandLineArgument("-report", reportFile);

          std::string dlm = ";";
          std::string header = dlm;
          std::string scores = sampleName + dlm;

          header += "Probability" + dlm;
          scores += std::to_string(probability) + dlm;

          header += "Mean, mm" + dlm;
          scores += std::to_string(mean) + dlm;

          header += "RMSE, mm" + dlm;
          scores += std::to_string(rmse) + dlm;

          header += "Maximal, mm" + dlm;
          scores += std::to_string(maximal) + dlm;

          bool fileExist = boost::filesystem::exists(reportFile);

          std::ofstream ofile;
          ofile.open(reportFile, std::ofstream::out | std::ofstream::app);

          if (!fileExist) {
            ofile << header << std::endl;
          }
          ofile << scores << std::endl;
          ofile.close();
        }
      }
    }
  }
  catch (statismo::StatisticalModelException& e) {
    std::cout << "Exception occurred while building the shape model" << std::endl;
    std::cout << e.what() << std::endl;
  }
}

