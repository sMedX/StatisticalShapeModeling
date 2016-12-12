#include <itkMesh.h>
#include <itkImage.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"
#include "utils/PointSetToImageMetrics.h"
#include "SurfaceToLevelSetImageFilter.h"
#include "itkSurfaceToImageRegistrationFilter.h"

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;
using fp = boost::filesystem::path;

int main(int argc, char** argv) {

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();

  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string meanSurfaceFile;
  parser->GetCommandLineArgument("-output", meanSurfaceFile);

  std::cout << std::endl;
  std::cout << "parameters" << std::endl;
  std::cout << "    list of files " << listFile << std::endl;
  std::cout << "   output surface " << meanSurfaceFile << std::endl;
  std::cout << std::endl;

  // read surface in to vector
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "Could not read the data-list:" << endl;
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (StringList::const_iterator it = listOfFiles.begin(); it != listOfFiles.end(); ++it) {
    std::string fileName = it->c_str();
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    std::cout << "input surface polydata info " << fileName << std::endl;
    std::cout << " number of cells " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;
  }

  std::cout << "number of surfaces " << vectorOfSurfaces.size() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // build mean surface

  MeshType::Pointer meanSurface = vectorOfSurfaces[0];

  for (int n = 0; n < meanSurface->GetNumberOfPoints(); ++n) {
    MeshType::PointType point;
    point.Fill(0);

    for (int s = 0; s < vectorOfSurfaces.size(); ++s) {
      for (int d = 0; d < MeshType::PointDimension; ++d) {
        point[d] += vectorOfSurfaces[s]->GetPoint(n)[d];
      }
    }

    for (int d = 0; d < MeshType::PointDimension; ++d) {
      point[d] /= vectorOfSurfaces.size();
    }

    meanSurface->SetPoint(n, point);
  }

  std::cout << "write surface to the file " << meanSurfaceFile << std::endl;
  if (!writeMesh<MeshType>(meanSurface, meanSurfaceFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


