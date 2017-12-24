#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToPointSetMetrics.h"
#include "ssmModelSpecificityOptions.h"

double computeSpecificity(ShapeModelType::Pointer model, const MeshVectorType & vectorOfSurfaces, const size_t & numberOfSamples);

int main(int argc, char** argv)
{
  // parse options
  ssm::ModelSpecificityOptions options;

  if (!options.ParseOptions(argc, argv)) {
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // read statistical shape model
  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  auto model = StatisticalModelType::New();
  try {
    model->Load(RepresenterType::New(), options.GetModelFileName().c_str());
  }
  catch (itk::ExceptionObject & excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  // read list of files
  StringVector listOfFiles;
  try {
    listOfFiles = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  MeshVectorType vectorOfShapes;
  for (const auto & fileName : listOfFiles) {
    auto surface = MeshType::New();
    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfShapes.push_back(surface);
  }

  std::cout << "model " << options.GetModelFileName() << std::endl;
  std::cout << "number of random shapes " << options.GetNumberOfSamples() << std::endl;
  std::cout << "number of train shapes  " << vectorOfShapes.size() << std::endl;
  std::cout << std::endl;

  // compute specificity
  double specificity = computeSpecificity(model, vectorOfShapes, options.GetNumberOfSamples());

  std::cout << "specificity of the model " << specificity << std::endl;
  std::cout << std::endl;

  // print report to *.csv file
  if (options.GetReportFileName().size() > 0) {
    std::cout << "print report to the file " << options.GetReportFileName() << std::endl;
    std::ofstream file(options.GetReportFileName(), std::ofstream::out);
    file << "model;" << options.GetModelFileName() << std::endl;
    file << "number of train shapes;" << vectorOfShapes.size() << std::endl;
    file << "number of random shapes;" << options.GetNumberOfSamples() << std::endl;
    file << "specificity;" << specificity << std::endl;
    file.close();
  }
}

double computeSpecificity(ShapeModelType::Pointer model, const MeshVectorType & vectorOfShapes, const size_t & numberOfSamples)
{
  double specificity = 0;

  for (size_t count = 0; count < numberOfSamples; ++count) {

    double value = std::numeric_limits<double>::max();

    // generate random point set
    auto randomPoints = PointSetType::New();
    randomPoints->SetPoints(model->DrawSample()->GetPoints());

    for (const auto & shape : vectorOfShapes) {

      // get point set from train dataset
      auto trainPoints = PointSetType::New();
      trainPoints->SetPoints(shape->GetPoints());

      // compute metrics
      typedef ssm::PointSetToPointSetMetrics<PointSetType> PointSetToPointSetMetricsType;
      auto metrics = PointSetToPointSetMetricsType::New();
      metrics->SetFixedPointSet(randomPoints);
      metrics->SetMovingPointSet(trainPoints);
      metrics->Compute();

      value = std::min(value, metrics->GetRMSEValue());
    }

    specificity += value;

    std::cout << count + 1 << "/" << numberOfSamples << " distance: " << value << std::endl;
  }

  return specificity / numberOfSamples;
}
