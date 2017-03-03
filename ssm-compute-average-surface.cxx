#include <fstream>

#include <boost/program_options.hpp>
#include <utils/statismo-build-models-utils.h>
#include "utils/ssmTypes.h"
#include "utils/io.h"

struct ProgramOptions
{
  bool help;
  std::string listFile;
  std::string surfaceFile;
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

  // read surfaces into vector
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(options.listFile);
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

    std::cout << fileName << std::endl;
    std::cout << "number of cells  " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;

    if (vectorOfSurfaces[0]->GetNumberOfPoints() != surface->GetNumberOfPoints()) {
      std::cout << "The number of points must be the same for all surfaces" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "number of surfaces " << vectorOfSurfaces.size() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // compute the average surface
  MeshType::Pointer average = vectorOfSurfaces[0];

  for (size_t n = 0; n < average->GetNumberOfPoints(); ++n) {
    MeshType::PointType point(0);

    for (size_t s = 0; s < vectorOfSurfaces.size(); ++s) {
      for (size_t d = 0; d < MeshType::PointDimension; ++d) {
        point[d] += vectorOfSurfaces[s]->GetPoint(n)[d];
      }
    }

    for (size_t d = 0; d < MeshType::PointDimension; ++d) {
      point[d] /= vectorOfSurfaces.size();
    }

    average->SetPoint(n, point);
  }

  std::cout << "write surface to the file " << options.surfaceFile << std::endl;
  if (!writeMesh<MeshType>(average, options.surfaceFile)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("list,l", po::value<std::string>(&options.listFile), "The path to the file with list of surfaces to average.")
    ("surface,s", po::value<std::string>(&options.surfaceFile), "The path for the output surface.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(help);

  return description;
}
