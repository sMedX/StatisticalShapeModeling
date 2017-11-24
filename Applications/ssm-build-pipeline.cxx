#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "ssmTypes.h"
#include "ssmUtils.h"
//#include <utils/statismo-build-models-utils.h>
#include "ssmBinaryMask3DMeshSource.h"


struct ProgramOptions
{
  bool help;
  std::string iniFile;
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

  std::cout << options.iniFile << std::endl;

  boost::property_tree::ptree ptree;
  boost::property_tree::ini_parser::read_ini(options.iniFile, ptree);

  std::cout << ptree.get<std::string>("FILES.list") << std::endl;

  // read list of files
  StringVector listOfFiles;
  try {
    listOfFiles = readListFromFile(ptree.get<std::string>("FILES.list"));
  }
  catch (ifstream::failure & e) {
    std::cerr << "Could not read the list of files: "  << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  for (const auto & fileName : listOfFiles) {
    std::cout << fileName << std::endl;

    // read image
    auto image = BinaryImageType::New();
    if (!readImage<BinaryImageType>(image, fileName)) {
      return EXIT_FAILURE;
    }

    typedef ssm::BinaryMask3DMeshSource<BinaryImageType, vtkPolyData> BinaryMask3DMeshSourceType;
    auto extract = BinaryMask3DMeshSourceType::New();
    extract->SetInput(image);
    try {
      extract->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      return EXIT_FAILURE;
    }

    vtkSmartPointer<vtkPolyData> surface = extract->GetOutput();
    std::cout << "       number of cells  " << surface->GetNumberOfCells() << std::endl;
    std::cout << "       number of points " << surface->GetNumberOfPoints() << std::endl;
  }

  return EXIT_SUCCESS;
}



po::options_description initializeProgramOptions(ProgramOptions& options)
{
  po::options_description mandatory("Mandatory options");
  mandatory.add_options()
    ("config,c", po::value<std::string>(&options.iniFile), "The path to the config ini file.")
    ;

  po::options_description help("Optional options");
  help.add_options()
    ("help,h", po::bool_switch(&options.help), "Display this help message")
    ;

  po::options_description description;
  description.add(mandatory).add(help);

  return description;
}
