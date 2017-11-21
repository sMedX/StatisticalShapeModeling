#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/format.hpp>

#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmBinaryImageToLevelSetImageFilter.h"
#include "ssmBinaryMask3DMeshSource.h"

struct ProgramOptions
{
  bool help;
  std::string configFile;
  std::string imageFile;
  std::string outputFile;
  std::string reportFile;
  double sigma = 0;
  double relaxation = 0.2;
  size_t iterations = 100;
  size_t points = 0;
};

double averageLengthOfEdges(vtkPolyData* poly);
double averageAreaOfCells(vtkPolyData* poly);
void  extractSurface(const ProgramOptions & options);

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
  StringList listOfFiles;

  if ( configIsDesabled ) {
    extractSurface(options);
    return EXIT_SUCCESS;
  }
  else {
    const std::string group = "EXTRACTION.";
    pt::ptree ptree;

    try {
      pt::ini_parser::read_ini(options.configFile, ptree);
    }
    catch (const pt::ptree_error &e) {
      cerr << "An exception occurred while parsing the ini file:" << options.configFile << endl;
      cout << e.what() << endl;
    }

    ProgramOptions op;
    options = op;

    // read list of files
    try {
      listOfFiles = readListOfFiles(ptree.get<std::string>(group + "list"));
    }
    catch (ifstream::failure & e) {
      std::cerr << "Could not read the list of files: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }

    try {
      options.sigma = ptree.get<double>(group + "sigma");
      options.relaxation = ptree.get<double>(group + "relaxation");
      options.iterations = ptree.get<size_t>(group + "iterations");
      options.points = ptree.get<size_t>(group + "points");
      options.reportFile = ptree.get<std::string>(group + "report");
      options.outputFile = ptree.get<std::string>(group + "output");
    }
    catch (...) {
    }
  }

  // extract surfaces
  const std::string format = options.outputFile;

  for (const auto & imageFile : listOfFiles) {
    options.imageFile = imageFile;
    options.outputFile = (boost::format(format) % getBaseNameFromPath(options.imageFile)).str();

    std::cout << options.imageFile << std::endl;
    std::cout << options.outputFile << std::endl;
    extractSurface(options);
  }

  return EXIT_SUCCESS;
}

void extractSurface(const ProgramOptions & options )
{
  // read image
  auto image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, options.imageFile)) {
    throw;
  }

  typedef ssm::BinaryMask3DMeshSource<BinaryImageType, vtkPolyData> BinaryMask3DMeshSourceType;
  auto binaryMaskToSurface = BinaryMask3DMeshSourceType::New();
  binaryMaskToSurface->SetInput(image);
  try {
    binaryMaskToSurface->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }
  auto surface = binaryMaskToSurface->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, options.outputFile)) {
    throw;
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
    throw;
  }
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (options.reportFile!="") {
    std::cout << "print report to the file: " << options.reportFile << std::endl;
    metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(options.outputFile));
  }

  return;
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
    ("image,i", po::value<std::string>(&options.imageFile), "The path to the input image file.")
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
