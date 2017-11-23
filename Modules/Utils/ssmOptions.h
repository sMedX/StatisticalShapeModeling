#pragma once

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/format.hpp>

namespace pt = boost::property_tree;
namespace po = boost::program_options;

//=========================================================================
// Some basic functions
//=========================================================================
void printTree(const pt::ptree & tree, std::ostream & os, unsigned int level /*=0*/)
{
  if (!tree.empty()) {
    os << std::endl;
    std::string indent(3 * level, ' ');

    for (const auto & it : tree) {
      std::string name = it.first;
      name.resize(12);
      os << indent << name << " ";

      printTree(it.second, os, level + 1);
      os << std::endl;
    }
  }
  else {
    std::cout << " " << tree.data();
  }
  return;
}

std::string AddQuotes(std::string str)
{
  return "'" + str + "'";
}

//=========================================================================
// Base options class
//=========================================================================
class OptionsBase
{
public:

  const bool & ConfigIsEnabled() const 
  { 
    return configIsEnabled; 
  }

  void SetGroup(const std::string & str) 
  { 
    group = str; 
  }

  bool ParseCommandLine(int argc, char** argv)
  {
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
    if (help) {
      cout << description << endl;
      return false;
    }

    configIsEnabled = !vm["config"].empty();
    return true;
  }

  bool ReadConfigFile()
  {
    try {
      pt::ini_parser::read_ini(config, parsedPtree);
    }
    catch (const pt::ptree_error &e) {
      std::cerr << "An exception occurred while parsing the config file:" << config << std::endl;
      std::cout << e.what() << endl;
      return false;
    }

    if (parsedPtree.find(group) == parsedPtree.not_found()) {
      std::cerr << "The group " << group << " is not found in the config file: " << config << std::endl;
      return false;
    }

    parsedPtree = parsedPtree.get_child(group);
  }

  void PrintConfig()
  {
    std::cout << std::endl;
    std::cout << "Config options for group " << group << std::endl;
    printTree(parsedPtree, std::cout, 0);
    std::cout << std::endl;
  };

  template <typename T>
  T GetDefaultValue(const std::string & str) const
  {
    return ptreeOfDefaultValues.get<T>(str);
  };

  template <typename T>
  T Get(const std::string & str) const
  {
    return parsedPtree.get<T>(str);
  };

protected:
  OptionsBase()
  {
    help = false;
    configIsEnabled = false;

    po::options_description configOptions("Optional config options");
    configOptions.add_options()("config,c", po::value<std::string>(&config), "The path to the config file.");

    po::options_description helpOptions("Optional help options");
    helpOptions.add_options()("help,h", po::bool_switch(&help)->default_value(help), "Display this help message");

    description.add(configOptions).add(helpOptions);
  }

  template <typename T>
  void Put(const std::string & str, const T & value, const bool & required = true)
  {
    ptreeOfRequireds.put(str, required);
    ptreeOfDefaultValues.put(str, value);
  };

  std::string Path(const std::string & str) const
  {
    return group + "." + str;
  }

  bool Find(const std::string & str) const
  {
    if (parsedPtree.find(str) == parsedPtree.not_found()) {
      std::cerr << "Key " << AddQuotes(Path(str)) << " is missing in the config file: " << config << std::endl;
      return false;
    }
    return true;
  };

  bool help;
  std::string config;
  std::string group;
  bool configIsEnabled;

  po::variables_map vm;
  po::options_description description;
  pt::ptree parsedPtree;
  pt::ptree ptreeOfRequireds;
  pt::ptree ptreeOfDefaultValues;
};

//=========================================================================
// Surface extraction options
//=========================================================================
class SurfaceExtractionOptions : public OptionsBase
{
public:

  void SetInputFileName(const std::string & str) { inputFileName = str; }
  const std::string & GetInputFileName() const { return inputFileName; } 

  void SetOutputFileName(const std::string & str) { outputFileName = str; }
  const std::string & GetOutputFileName() const { return outputFileName; }

  std::string GetInputList() const
  { 
    return parsedPtree.get<std::string>("inplist"); 
  }

  std::string GetOutputList() const
  {
    return parsedPtree.get<std::string>("outlist");
  }

  std::string GetReportFile() const
  {
    if (configIsEnabled)
      return parsedPtree.get<std::string>("report");
    else
      return vm["report"].as<std::string>();
  }

  double GetSigma() const
  {
    const std::string name = "sigma";
    if (configIsEnabled)
      return Get<double>(name);
    else
      return vm[name].as<double>();
  }

  double GetFactor() const
  {
    const std::string name = "factor";
    if (configIsEnabled)
      return Get<double>(name);
    else
      return vm[name].as<double>();
  }

  size_t GetNumberOfPoints() const
  {
    const std::string name = "points";
    if (configIsEnabled)
      return Get<size_t>(name);
    else
      return vm[name].as<size_t>();
  }

  size_t GetNumberOfIterations() const
  {
    const std::string name = "iterations";
    if (configIsEnabled)
      return Get<size_t>(name);
    else
      return vm[name].as<size_t>();
  }

  std::string FormatOutput(const std::string & fileName)
  {
    const auto format = Get<std::string>("output");
    try {
      return (boost::format(format) % getBaseNameFromPath(fileName)).str();
    }
    catch (const boost::io::format_error &e) {
      std::cerr << "Could not format string with format " << format << std::endl;
      std::cout << e.what() << std::endl;
      throw;
    }
  };

  SurfaceExtractionOptions() 
  {
    SetGroup("EXTRACTION");

    // initialize ptree
    Put<std::string>("inplist", "");
    Put<std::string>("outlist", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");

    Put<double>("sigma", 0, 0);
    Put<double>("factor", 0.2, 0);
    Put<size_t>("iterations", 100, 0);
    Put<size_t>("points", 0, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("input,i", po::value<std::string>(&inputFileName), "The path to the input image file.")
      ("output,o", po::value<std::string>(&outputFileName), "The path for the output surface file.")
      ;
    
    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("sigma", po::value<double>()->default_value(this->GetDefaultValue<double>("sigma")), "The sigma of the Gaussian kernel measured in world coordinates.")
      ("factor", po::value<double>()->default_value(this->GetDefaultValue<double>("factor")), "The relaxation factor for Laplacian smoothing.")
      ("iterations", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("iterations")), "The number of iterations.")
      ("points", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("points")), "The number of points in output surface.")
      ;

    po::options_description reportOptions("Optional report options");
    reportOptions.add_options()
      ("report,r", po::value<std::string>(), "The path for the file to print report.")
      ;

    description.add(mandatoryOptions).add(inputOptions).add(reportOptions);
  };

private:
  std::string inputFileName;
  std::string outputFileName;
};
