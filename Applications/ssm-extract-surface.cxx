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
bool extractSurface(const std::string & inputFileName, const std::string & outputFileName, const ssm::ExtractionOptions & options);

int main(int argc, char** argv)
{
  // parse options
  ssm::ExtractionOptions options;

  if (!options.ParseOptions(argc, argv)) {
    return EXIT_FAILURE;
  }

  if ( !options.ConfigIsEnabled() ) {
    if (!extractSurface(options.GetInputFileName(), options.GetOutputFileName(), options)) {
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }

  //----------------------------------------------------------------------------
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

  for (const auto & inputFileName : listOfInputFiles) {
    std::string outputFileName = options.FormatOutput(inputFileName);

    if(!extractSurface(inputFileName, outputFileName, options)) {
      return EXIT_FAILURE;
    }

    listOfOutputFiles.push_back(outputFileName);
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

bool extractSurface(const std::string & inputFileName, const std::string & outputFileName, const ssm::ExtractionOptions & options )
{
  // read image
  auto image = BinaryImageType::New();
  if (!readImage<BinaryImageType>(image, inputFileName)) {
    return false;
  }
  printImageInfo<BinaryImageType>(image, inputFileName);

  typedef ssm::Image3DMeshSource<BinaryImageType, vtkPolyData> Image3DMeshSourceType;
  auto binaryMaskToSurface = Image3DMeshSourceType::New();
  binaryMaskToSurface->SetInput(image);
  try {
    binaryMaskToSurface->SetSigma(options.GetSigma());
    binaryMaskToSurface->SetNumberOfIterations(options.GetNumberOfIterations());
    binaryMaskToSurface->SetNumberOfPoints(options.GetNumberOfPoints());
  }
  catch (...) {
    return false;
  }

  try {
    binaryMaskToSurface->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return false;
  }
  auto surface = binaryMaskToSurface->GetOutput();

  // write polydata to the file
  if (!writeVTKPolydata(surface, outputFileName)) {
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
  metrics->PrintReportToFile(options.GetReportFileName(), getBaseNameFromPath(options.GetOutputFileName()));

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
