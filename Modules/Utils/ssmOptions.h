#pragma once

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/format.hpp>

namespace pt = boost::property_tree;
namespace po = boost::program_options;

class SurfaceExtractionOptions
{
public:

  bool help;

  bool configIsEnabled;
  std::string configFile;
  std::string group;
  std::string format;

  std::string inp_list;
  std::string out_list;
  std::string inputFile;
  std::string outputFile;
  std::string reportFile;

  double sigma;
  double relaxation;
  size_t iterations;
  size_t points;

  bool ParseCommandLine(int argc, char** argv)
  {
    po::variables_map vm;
    try {
      po::parsed_options parsedOptions = po::command_line_parser(argc, argv).options(description).run();
      po::store(parsedOptions, vm);
      po::notify(vm);
    }
    catch (const po::error& e) {
      cerr << "An exception occurred while parsing the command line." << endl;
      cerr << e.what() << endl;
      cout << description << endl;
      return false;
    }
    if (help == true) {
      cout << description << endl;
      return false;
    }

    configIsEnabled = !vm["config"].empty();

    return true;
  }

  bool ReadOptions()
  {
    SurfaceExtractionOptions();
    pt::ptree ptree;

    try {
      pt::ini_parser::read_ini(configFile, ptree);
    }
    catch (const pt::ptree_error &e) {
      std::cerr << "An exception occurred while parsing the config file:" << configFile << std::endl;
      std::cout << e.what() << endl;
      return false;
    }

    if (ptree.find(group) == ptree.not_found()) {
      std::cerr << "The group " << group << " is not found in the config file: " << configFile << std::endl;
      return false;
    }

    try { inp_list = ptree.get<std::string>(group + ".inp_list"); }
    catch (...) {}

    try { out_list = ptree.get<std::string>(group + ".out_list"); }
    catch (...) {}

    try { sigma = ptree.get<double>(group + ".sigma"); }
    catch (...) {}

    try {relaxation = ptree.get<double>(group + ".relaxation"); }
    catch (...) {}

    try { iterations = ptree.get<size_t>(group + ".iterations"); }
    catch (...) {}

    try { points = ptree.get<size_t>(group + ".points"); }
    catch (...) {}

    try { reportFile = ptree.get<std::string>(group + ".report"); }
    catch (...) {}

    try { format = ptree.get<std::string>(group + ".output"); }
    catch (...) {}

    return true;
  }

  std::string FormatOutput(const std::string & fileName)
  {
    try {
      return (boost::format(format) % getBaseNameFromPath(fileName)).str();
    }
    catch (const boost::io::format_error &e) {
      std::cerr << "Could not format with format string " << "'" << format << "'" << std::endl;
      std::cout << e.what() << std::endl;
      throw;
    }
  };

  void PrintOptions()
  {
    std::cout << inputFile << std::endl;
    std::cout << outputFile << std::endl;
    std::cout << std::endl;
    std::cout << group << std::endl;
    std::cout << sigma << std::endl;
    std::cout << relaxation << std::endl;
    std::cout << iterations << std::endl;
    std::cout << points << std::endl;
    std::cout << std::endl;
  };

  SurfaceExtractionOptions() 
  {
    group = "EXTRACTION";

    help = false;
    sigma = 0;
    relaxation = 0.2;
    iterations = 100;
    points = 0;

    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("config,c", po::value<std::string>(&configFile), "The path to the config file.")
      ("input,i", po::value<std::string>(&inputFile), "The path to the input image file.")
      ("output,o", po::value<std::string>(&outputFile), "The path for the output surface file.")
      ;

    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("sigma", po::value<double>(&sigma)->default_value(sigma), "The sigma of the Gaussian kernel measured in world coordinates.")
      ("factor", po::value<double>(&relaxation)->default_value(relaxation), "The relaxation factor for Laplacian smoothing.")
      ("iterations", po::value<size_t>(&iterations)->default_value(iterations), "The number of iterations.")
      ("points", po::value<size_t>(&points)->default_value(points), "The number of points in the output surface.")
      ;

    po::options_description reportOptions("Optional report options");
    reportOptions.add_options()
      ("report,r", po::value<std::string>(&reportFile), "The path for the file to print report.")
      ;

    po::options_description helpOptions("Optional options");
    helpOptions.add_options()
      ("help,h", po::bool_switch(&help), "Display this help message")
      ;

    description.add(mandatoryOptions).add(inputOptions).add(reportOptions).add(helpOptions);
  };

private:
  po::options_description description;
};