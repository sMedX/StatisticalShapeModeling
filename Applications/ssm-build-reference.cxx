#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkCell.h>

#include "ssmTypes.h"
#include "ssmUtils.h"
#include "ssmPointSetToImageMetrics.h"
#include "ssmBinaryMask3DMeshSource.h"
#include "ssmReferenceOptions.h"

double averageLengthOfEdges(vtkPolyData* poly);
double averageAreaOfCells(vtkPolyData* poly);
bool extractSurface(const ssm::ReferenceOptions & options);

int main(int argc, char** argv)
{
  ssm::ReferenceOptions options;

  if (!options.ParseCommandLine(argc, argv)) {
    return EXIT_FAILURE;
  }

  // if config is enabled read options from config file
  if ( options.ConfigIsEnabled() ) {
    if (!options.ReadConfigFile()) {
      return EXIT_FAILURE;
    }
    options.PrintConfig();
  }

  // build reference
  extractSurface(options);

  return EXIT_SUCCESS;
}

bool extractSurface(const ssm::ReferenceOptions & options )
{
  // read image
  auto image = FloatImageType::New();
  if (!readImage<FloatImageType>(image, options.GetInputFileName())) {
    return false;
  }
  printImageInfo<FloatImageType>(image, options.GetInputFileName());

  typedef ssm::BinaryMask3DMeshSource<FloatImageType, vtkPolyData> BinaryMask3DMeshSourceType;
  auto binaryMaskToSurface = BinaryMask3DMeshSourceType::New();
  binaryMaskToSurface->SetInput(image);
  binaryMaskToSurface->SetSigma(options.GetSigma());
  binaryMaskToSurface->SetLevelValue(0);
  binaryMaskToSurface->SetComputeLevelValue(false);
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
  // compute metrics
  typedef itk::PointSet<float, 3> PointSetType;
  typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
  auto metrics = PointSetToImageMetricsType::New();
  metrics->SetMovingImage(image);
  metrics->SetPointSetAsPolyData<vtkPolyData>(surface);
  metrics->SetInfo(surfaceInfo);
  try {
    metrics->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    return false;
  }

  metrics->PrintReport(std::cout);

  // write report to *.csv file
  std::cout << "print report to the file: " << options.GetReportFileName() << std::endl;
  std::cout << std::endl;

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
