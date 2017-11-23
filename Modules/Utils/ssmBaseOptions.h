#pragma once

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/format.hpp>

namespace pt = boost::property_tree;
namespace po = boost::program_options;

namespace ssm
{
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

void checkParsedTree(const pt::ptree & ptreeOfRequireds, pt::ptree & parsedPtree, std::string & key, std::vector<std::string> & list)
{
  if (ptreeOfRequireds.empty()) {
    return;
  }

  for (const auto & it : ptreeOfRequireds) {
    const auto &name = it.first;
    const auto &tree = it.second;

    if (!tree.empty()) {
      key = key + "." + name;
      checkParsedTree(ptreeOfRequireds.get_child(name), parsedPtree.get_child(name), key, list);
      return;
    }
    else {
      if (ptreeOfRequireds.get<bool>(name)) {
        if (parsedPtree.find(name) == parsedPtree.not_found()) {
          list.push_back(key + "." + name);
          continue;
        }
      }
    }
  }
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
      std::cerr << "An exception occurred while parsing the config file:" << AddQuotes(config) << std::endl;
      std::cout << e.what() << endl;
      return false;
    }

    if (parsedPtree.find(group) == parsedPtree.not_found()) {
      std::cerr << "The group " << group << " is not found in the config file: " << AddQuotes(config) << std::endl;
      return false;
    }

    parsedPtree = parsedPtree.get_child(group);

    // check parsed ptree
    std::vector<std::string> listOfKeys;
    checkParsedTree(ptreeOfRequireds, parsedPtree, group, listOfKeys);

    if (listOfKeys.size() > 0) {
      std::cerr << "The keys are not found in the config file: " << AddQuotes(config) << std::endl;
      for (const auto & str : listOfKeys) {
        std::cerr << AddQuotes(str) << std::endl;
      }
      return false;
    }

    return true;
  }

  void PrintConfig()
  {
    std::cout << std::endl;
    std::cout << "Config options for group " << AddQuotes(group) << std::endl;
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
    if (!required) {
      parsedPtree.put(Path(str), value);
    }
  };

  std::string Path(const std::string & str) const
  {
    return group + "." + str;
  }

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
}
