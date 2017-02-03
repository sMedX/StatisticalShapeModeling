#include <boost/program_options.hpp>

#include <vtkPolyData.h>
#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <itkImageToVTKImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmBinaryImageToLevelSetImageFilter.h"

double averageLengthOfEdges(vtkPolyData* poly);
double averageAreaOfCells(vtkPolyData* poly);

struct ProgramOptions
{
  bool help;
  std::string imageFile;
  std::string surfaceFile;
  std::string reportFile;
  double sigma = 1.0;
  double relaxation = 0.2;
  size_t numberOfIterations = 100;
  size_t numberOfPoints = 0;
  bool isbinary = true;
  double levelValue = std::numeric_limits<double>::lowest();
};

namespace po = boost::program_options;
po::options_description initializeProgramOptions(ProgramOptions& poParameters);

int main(int argc, char** argv) {
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

  //----------------------------------------------------------------------------
  // read image
  FloatImageType::Pointer image = FloatImageType::New();
  if (!readImage<FloatImageType>(image, options.imageFile)) {
    return EXIT_FAILURE;
  }

  // compute the minimum and the maximum intensity values of label
  typedef itk::MinimumMaximumImageCalculator <FloatImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
  labelValues->SetImage(image);
  labelValues->Compute();

  if (options.levelValue < labelValues->GetMinimum() || options.levelValue > labelValues->GetMaximum()) {
    options.levelValue = 0.5*(labelValues->GetMinimum() + labelValues->GetMaximum());
  }

  std::cout << "input image " << options.imageFile << std::endl;
  std::cout << "       size " << image->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "    spacing " << image->GetSpacing() << std::endl;
  std::cout << "     origin " << image->GetOrigin() << std::endl;
  std::cout << "level value " << options.levelValue << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // image processing

  // smoothing
  typedef itk::RecursiveGaussianImageFilter<FloatImageType, FloatImageType> RecursiveGaussianImageFilterType;
  RecursiveGaussianImageFilterType::Pointer gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(image);
  gaussian->SetSigma(options.sigma);

  // fill holes after smoothing
  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(gaussian->GetOutput());
  fillholes->SetFullyConnected(true);
  try {
    fillholes->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------
  // convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(fillholes->GetOutput());
  convertor->Update();

  typedef vtkSmartPointer<vtkMarchingCubes> MarchingCubes;
  MarchingCubes mcubes = MarchingCubes::New();
  mcubes->SetInputData(convertor->GetOutput());
  mcubes->SetValue(0, options.levelValue);
  try {
    mcubes->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  vtkSmartPointer<vtkPolyData> surface = mcubes->GetOutput();

  // decimate surface
  if (options.numberOfPoints > 0) {
    double reduction = 1 - (options.numberOfPoints - 1) / (double)surface->GetNumberOfPoints();
    std::cout << "reduction to decimate surface " << reduction << std::endl;
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(surface);
    decimate->SetTargetReduction(reduction);
    decimate->SetPreserveTopology(true);
    decimate->SetSplitting(false);
    try {
      decimate->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    surface = decimate->GetOutput();
  }

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(surface);
  smoother->SetNumberOfIterations(options.numberOfIterations);
  smoother->SetRelaxationFactor(options.relaxation);
  try {
    smoother->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  typedef vtkSmartPointer<vtkPolyDataNormals> PolyDataNormals;
  PolyDataNormals normals = PolyDataNormals::New();
  normals->SetInputData(smoother->GetOutput());
  normals->AutoOrientNormalsOn();
  normals->FlipNormalsOff();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->SplittingOff();
  try {
    normals->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  surface = normals->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, options.surfaceFile)) {
    return EXIT_FAILURE;
  }

  typedef std::pair<std::string, std::string> PairType;
  std::vector<PairType> surfaceInfo;
  surfaceInfo.push_back(PairType("number of points", std::to_string(surface->GetNumberOfPoints())));
  surfaceInfo.push_back(PairType("number of cells", std::to_string(surface->GetNumberOfCells())));
  surfaceInfo.push_back(PairType("length of edges", std::to_string(averageLengthOfEdges(surface))));
  surfaceInfo.push_back(PairType("area of cells", std::to_string(averageAreaOfCells(surface))));

  std::cout << "output surface polydata " << options.surfaceFile << std::endl;
  std::cout << "       number of cells  " << surface->GetNumberOfCells() << std::endl;
  std::cout << "       number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << "average length of edges " << averageLengthOfEdges(surface) << endl;
  std::cout << "  average area of cells " << averageAreaOfCells(surface) << endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // setup point set
  typedef itk::PointSet<float, MeshType::PointDimension> PointSetType;
  PointSetType::Pointer pointSet = PointSetType::New();
  for (size_t n = 0; n < surface->GetPoints()->GetNumberOfPoints(); ++n) {
    PointSetType::PointType point;
    point.CastFrom<double>(surface->GetPoints()->GetPoint(n));
    pointSet->SetPoint(n, point);
  }

  // compute metrics
  typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
  metrics->SetFixedPointSet(pointSet);
  metrics->SetInfo(surfaceInfo);
  if (options.isbinary) {
    // compute level set image
    typedef ssm::BinaryImageToLevelSetImageFilter<FloatImageType, FloatImageType> BinaryImageToLevelSetImageType;
    BinaryImageToLevelSetImageType::Pointer levelset = BinaryImageToLevelSetImageType::New();
    levelset->SetInput(image);
    try {
      levelset->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }
    metrics->SetMovingImage(levelset->GetOutput());
  }
  else {
    metrics->SetMovingImage(image);
  }
  metrics->Compute();
  metrics->PrintReport(std::cout);

  // write report to *.csv file
  if (options.reportFile!="") {
    std::cout << "print report to the file: " << options.reportFile << std::endl;
    metrics->PrintReportToFile(options.reportFile, getBaseNameFromPath(options.surfaceFile));
  }

  return EXIT_SUCCESS;
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

  return sum / (3*numberOfCells);
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
    ("image,i", po::value<std::string>(&options.imageFile), "The path to the input image file.")
    ("output,o", po::value<std::string>(&options.surfaceFile), "The path to the output surface file.")
    ;

  po::options_description output("Optional input options");
  output.add_options()
    ("sigma,s", po::value<double>(&options.sigma), "The sigma of the Gaussian kernel measured in world coordinates.")
    ("factor,f", po::value<double>(&options.relaxation), "The relaxation factor for Laplacian smoothing.")
    ("iterations,t", po::value<size_t>(&options.numberOfIterations), "The number of iterations.")
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
  description.add(mandatory).add(output).add(report).add(help);

  return description;
}
