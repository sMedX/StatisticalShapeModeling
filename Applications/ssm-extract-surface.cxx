#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmBinaryImageToLevelSetImageFilter.h"
#include "ssmImage3DMeshSource.h"
#include "ssmExtractionOptions.h"

double averageLengthOfEdges(vtkPolyData* poly);
double averageAreaOfCells(vtkPolyData* poly);
bool extractSurface(const ssm::ExtractionOptions & options);

int main(int argc, char** argv)
{
  ssm::ExtractionOptions options;

  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  if ( !options.ConfigIsEnabled() ) {
    extractSurface(options);
    return EXIT_SUCCESS;
  }

  // read options from config file
  if (!options.ParseConfigFile()) {
    return EXIT_FAILURE;
  }
  options.PrintConfig();

  // read list of files
  StringVector listOfInputFiles;
  try {
    listOfInputFiles = readListFromFile(options.GetInputList());
  }
  catch (std::ifstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  // extract surfaces
  StringVector listOfOutputFiles;

  for (const auto & inputFile : listOfInputFiles) {
    options.SetInputFileName(inputFile);
    options.SetOutputFileName(options.FormatOutput(inputFile));

    if (extractSurface(options)) {
      listOfOutputFiles.push_back(options.GetOutputFileName());
    }
  }

  // write list of files
  try {
    writeListToFile(options.GetOutputList(), listOfOutputFiles);
  }
  catch (std::ofstream::failure & e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

bool extractSurface(const ssm::ExtractionOptions & options )
{
  // read image
  auto image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, options.GetInputFileName())) {
    return false;
  }
  printImageInfo<BinaryImageType>(image, options.GetInputFileName());

  typedef ssm::Image3DMeshSource<BinaryImageType, vtkPolyData> Image3DMeshSourceType;
  auto binaryMaskToSurface = Image3DMeshSourceType::New();
  binaryMaskToSurface->SetInput(image);
  binaryMaskToSurface->SetSigma(options.GetSigma());
  binaryMaskToSurface->SetNumberOfIterations(options.GetNumberOfIterations());
  binaryMaskToSurface->SetRelaxationFactor(options.GetFactor());
  binaryMaskToSurface->SetNumberOfPoints(options.GetNumberOfPoints());
  try {
    binaryMaskToSurface->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return false;
  }
  auto surface = binaryMaskToSurface->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, options.GetOutputFileName())) {
    return false;
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

  // compute metrics
  typedef itk::PointSet<float, 3> PointSetType;
  typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  auto metrics = PointSetToImageMetricsType::New();
  metrics->SetPointSetAsPolyData<vtkPolyData>(surface);
  metrics->SetImage(binaryImageToLevelset->GetOutput());
  metrics->SetInfo(surfaceInfo);
  try {
    metrics->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return false;
  }

  metrics->PrintReport();

  // write report to *.csv file
  std::cout << "print report to the file: " << options.GetReportFile() << std::endl;
  std::cout << std::endl;

  metrics->PrintReportToFile(options.GetReportFile(), getBaseNameFromPath(options.GetOutputFileName()));

  return true;
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
