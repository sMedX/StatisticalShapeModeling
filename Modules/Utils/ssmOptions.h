#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

namespace pt = boost::property_tree;

class ProgramOptions
{
public:

  bool help;

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

  int ReadOptions()
  {
    ProgramOptions();
    pt::ptree ptree;

    try {
      pt::ini_parser::read_ini(configFile, ptree);
    }
    catch (const pt::ptree_error &e) {
      std::cerr << "An exception occurred while parsing the config file:" << configFile << std::endl;
      std::cout << e.what() << endl;
      return 0;
    }

    if (ptree.find(group) == ptree.not_found()) {
      std::cerr << "The group " << group << " is not found in the config file: " << configFile << std::endl;
      return 0;
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

    return 1;
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

  ProgramOptions() 
  {
    group = "EXTRACTION";

    help = false;
    sigma = 0;
    relaxation = 0.2;
    iterations = 100;
    points = 0;
  };

private:
};