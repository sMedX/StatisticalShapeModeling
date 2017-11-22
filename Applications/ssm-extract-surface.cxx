#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmBinaryImageToLevelSetImageFilter.h"
#include "ssmBinaryMask3DMeshSource.h"
#include "ssmOptions.h"

double averageLengthOfEdges(vtkPolyData* poly);
double averageAreaOfCells(vtkPolyData* poly);
int extractSurface(const ProgramOptions & options);

namespace po = boost::program_options;
namespace pt = boost::property_tree;

po::options_description initializeProgramOptions(ProgramOptions& poParameters);

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
  catch (const po::error& e) {
    cerr << "An exception occurred while parsing the command line:" << endl;
    cerr << e.what() << endl;
    cout << description << endl;
    return EXIT_FAILURE;
  }
  if (options.help == true) {
    cout << description << endl;
    return EXIT_SUCCESS;
  }

  bool configIsDesabled = vm["config"].empty();

  if ( configIsDesabled ) {
    extractSurface(options);
    return EXIT_SUCCESS;
  }

  // read options from config file
  if (!options.ReadOptions()) {
    return EXIT_FAILURE;
  }
  options.PrintOptions();

  // read list of files
  StringList listOfInputFiles;

  try {
    listOfInputFiles = readListOfFiles(options.inp_list);
  }
  catch (std::ifstream::failure & e) {
    std::cerr << "Could not read list of files to the file: " << options.inp_list << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // extract surfaces
  StringList listOfOutputFiles;

  for (const auto & inputFile : listOfInputFiles) {
    options.inputFile = inputFile;
    options.outputFile = options.FormatOutput(inputFile);

    if (extractSurface(options)) {
      listOfOutputFiles.push_back(options.outputFile);
    }
  }

  // write list of files to file
  try {
    std::ofstream file(options.out_list, std::ofstream::out);
    for (const auto & outputFile : listOfOutputFiles) {
      std::cout << outputFile << std::endl;
      file << outputFile << std::endl;
    }
    file.close();
  }
  catch (std::ofstream::failure & e) {
    std::cerr << "Could not write the list of files to the file: " << options.out_list << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int extractSurface(const ProgramOptions & options )
{
  // read image
  auto image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, options.inputFile)) {
    return 0;
  }
  printImageInfo<BinaryImageType>(image, options.inputFile);

  typedef ssm::BinaryMask3DMeshSource<BinaryImageType, vtkPolyData> BinaryMask3DMeshSourceType;
  auto binaryMaskToSurface = BinaryMask3DMeshSourceType::New();
  binaryMaskToSurface->SetInput(image);
  try {
    binaryMaskToSurface->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return 0;
  }
  auto surface = binaryMaskToSurface->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, options.outputFile)) {
    return 0;
  }

  //----------------------------------------------------------------------------
  // compute report
  typedef std::pair<std::string, std::string> PairType;
  std::vector<PairType> surfaceInfo;
  surfaceInfo.push_back(PairType("number of points", std::to_string(surface->GetNumberOfPoints())));
  surfaceInfo.push_back(PairType(" number of cells", std::to_string(surface->GetNumberOfCells())));
  surfaceInfo.push_back(PairType(" length of edges", std::to_string(averageLengthOfEdges(surface))));
  surfaceInfo.push_back(PairType("   area of cells", std::to_string(averageAreaOfCells(surface))));

  //----------------------------------------------------------------------------
  // compute level set image
  typedef ssm::BinaryImageToLevelSetImageFilter<BinaryImageType, FloatImageType> BinaryImageToLevelSetImageType;
  auto binaryImageToLevelset = BinaryImageToLevelSetImageType::New();
  binaryImageToLevelset->SetInput(image);

  // setup point set
  typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
  auto pointSet = PointSetType::New();
  for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
    PointSetType::PointType point;
    point.CastFrom<double>(surface->GetPoints()->GetPoint(n));
    pointSet->SetPoint(n, point);
  }

  // compute metrics
  typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  auto metrics = PointSetToImageMetricsType::New();
  metrics->SetMovingImage(binaryImageToLevelset->GetOutput());
  metrics->SetFixedPointSet(pointSet);
  metrics->SetInfo(surfaceInfo);
  try {
    metrics->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return 0;
  }
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (options.reportFile!="") {
    std::cout << "print report to the file: " << options.reportFile << std::endl; 
    std::cout << std::endl;
    metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(options.outputFile));
  }

  return 1;
}


double averageLengthOfEdges(vtkPolyData*poly)
{
  const unsigned int numberOfCells = poly->GetNumberOfCells();
  double sum = 0;

  for (int n = 0; n < poly->GetNumberOfCells(); ++n) {
    double p1[3], p2[3], p3[3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(n);
    vtkSmartPointer<vtkPoints> points = cell->GetPoints();

    points->GetPoint(0, p1);
    points->GetPoint(1, p2);
    points->GetPoint(2, p3);

    sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + 
           sqrt(vtkMath::Distance2BetweenPoints(p1, p3)) +
           sqrt(vtkMath::Distance2BetweenPoints(p2, p3));
  }

  return sum / (3 * numberOfCells);
}

double averageAreaOfCells(vtkPolyData*poly)
{
  const size_t numberOfCells = poly->GetNumberOfCells();
  double sum = 0;

  for (size_t n = 0; n < numberOfCells; ++n) {
    double a[3], b[3], c[3];
    double p1[3], p2[3], p3[3];

    vtkSmartPointer<vtkCell> cell = poly->GetCell(n);
    vtkSmartPointer<vtkPoints> points = cell->GetPoints();

    points->GetPoint(0, p1);
    points->GetPoint(1, p2);
    points->GetPoint(2, p3);
    
    vtkMath::Subtract(p1, p2, a);
    vtkMath::Subtract(p1, p3, b);
    vtkMath::Cross(a, b, c);

    sum += sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) / 2;
  }

  return sum / numberOfCells;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("config,c", po::value<std::string>(&options.configFile), "The path to the config ini file.")
    ("image,i", po::value<std::string>(&options.inputFile), "The path to the input image file.")
    ("surface,s", po::value<std::string>(&options.outputFile), "The path for the output surface file.")
    ;

  po::options_description output("Optional input options");
  output.add_options()
    ("sigma", po::value<double>(&options.sigma)->default_value(options.sigma), "The sigma of the Gaussian kernel measured in world coordinates.")
    ("factor", po::value<double>(&options.relaxation)->default_value(options.relaxation), "The relaxation factor for Laplacian smoothing.")
    ("iterations", po::value<size_t>(&options.iterations)->default_value(options.iterations), "The number of iterations.")
    ("points", po::value<size_t>(&options.points)->default_value(options.points), "The number of points in the output surface.")
    ;

  po::options_description report("Optional report options");
  report.add_options()
    ("report,r", po::value<std::string>(&options.reportFile), "The path for the file to print report.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(output).add(report).add(help);

  return description;
}
