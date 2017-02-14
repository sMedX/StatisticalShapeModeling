#include <boost/program_options.hpp>
#include <boost/scoped_ptr.hpp>

#include <utils/statismo-build-models-utils.h>
#include <DataManager.h>
#include <itkStandardMeshRepresenter.h>
#include <itkStatisticalModel.h>
#include <itkPointsLocator.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "ssm/ssmPointSetToPointSetMetrics.h"

typedef MeshType::PointType PointType;
typedef statismo::DataManager<MeshType> DataManagerType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef DataManagerType::DataItemListType DataItemListType;
typedef std::vector<MeshType::Pointer> MeshVectorType;

struct ProgramOptions
{
  bool help;
  std::string modelFile;
  std::string listFile;
  std::string reportFile;
  size_t numberOfSamples = 10;
};

typedef boost::filesystem::path fp;
namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);
double computeSpecificity(StatisticalModelType::Pointer model, const MeshVectorType testMeshes, unsigned numberOfSamples);

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

  // read statistical shape model
  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  RepresenterType::Pointer representer = RepresenterType::New();

  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model = StatisticalModelType::New();
  try {
    model->Load(representer, options.modelFile.c_str());
  }
  catch (itk::ExceptionObject & excp) {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
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

  MeshVectorType vectorOfSurfaces;

  for (const auto & fileName : fileNames) {
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
  }

  std::cout << "model  " << options.modelFile << std::endl;
  std::cout << "number of surfaces  " << vectorOfSurfaces.size() << std::endl;
  std::cout << "number of samples  " << options.numberOfSamples << std::endl;
  std::cout << std::endl;

  //compute specificity
  double specificity = computeSpecificity(model, vectorOfSurfaces, options.numberOfSamples);

  std::cout << "Specificity of model:   " << specificity << std::endl;
  std::cout << std::endl;

  // print report to *.csv file
  if (options.reportFile != "") {
    std::cout << "print report to the file: " << options.reportFile << std::endl;
    std::ofstream file(options.reportFile, std::ofstream::out | std::ofstream::app);
    file << getBaseNameFromPath(options.modelFile) << ";" << specificity << ";" << std::endl;
    file.close();
  }
}

double computeSpecificity(StatisticalModelType::Pointer model, const MeshVectorType testMeshes, unsigned numberOfSamples)
{
  double accumulatedDistToClosestTrainingShape = 0;

  for (unsigned i = 0; i < numberOfSamples; i++) {
    MeshType::Pointer sample = model->DrawSample();

    std::cout << "  Sample: " << i + 1 << std::endl;

    double minDist = std::numeric_limits<double>::max();

    for (MeshVectorType::const_iterator it = testMeshes.begin(); it != testMeshes.end(); ++it) {
      MeshType::Pointer testMesh = it->GetPointer();

      typedef std::pair<std::string, std::string> PairType;
      std::vector<PairType> info;
      info.push_back(PairType("Components", std::to_string(model->GetNumberOfPrincipalComponents())));
      info.push_back(PairType("Probability", std::to_string(model->ComputeProbabilityOfDataset(testMesh))));

      // compute metrics
      typedef itk::PointSet<MeshType::PointType, MeshType::PointDimension> PointSetType;
      PointSetType::Pointer pointSet1 = PointSetType::New();
      pointSet1->SetPoints(testMesh->GetPoints());

      PointSetType::Pointer pointSet2 = PointSetType::New();
      pointSet2->SetPoints(sample->GetPoints());

      typedef ssm::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
      PointSetToPointSetMetricsType::Pointer metrics = PointSetToPointSetMetricsType::New();
      metrics->SetFixedPointSet(pointSet1);
      metrics->SetMovingPointSet(pointSet2);
      metrics->SetInfo(info);
      metrics->Compute();

      double dist = metrics->GetRMSEValue();

      if (dist < minDist) {
        minDist = dist;
      }
    }
    accumulatedDistToClosestTrainingShape += minDist;
    std::cout << "    Minimal average distance: " << minDist << std::endl;
    std::cout << std::endl;
  }

  double specificity = accumulatedDistToClosestTrainingShape / numberOfSamples;
  return specificity;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("list,l", po::value<std::string>(&options.listFile), "The path to the file with list of surfaces.")
    ("model,m", po::value<std::string>(&options.modelFile), "The path to the input shape model file.")
    ;

  po::options_description input("Optional input options");
  input.add_options()
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