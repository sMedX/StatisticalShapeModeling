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
void printTree(const pt::ptree & tree, std::ostream & os, unsigned int level = 0)
{
  if (!tree.empty()) {
    os << std::endl;
    std::string indent(3 * level, ' ');

    for (const auto & it : tree) {
      std::string path = it.first;
      path.resize(16);
      os << indent << path << " ";

      printTree(it.second, os, level + 1);
      os << std::endl;
    }
  }

  std::cout << " " << tree.data();

  return;
}

void checkParsedTree(const pt::ptree & ptreeOfDefaultValues, const pt::ptree & ptreeOfRequired, pt::ptree & parsedPtree, std::string & path, std::vector<std::string> & list)
{
  if (ptreeOfRequired.empty()) {
    return;
  }

  for (const auto & it : ptreeOfRequired) {
    const auto &name = it.first;
    const auto &tree = it.second;

    if (!tree.empty()) {
      path = path + "." + name;
      checkParsedTree(ptreeOfDefaultValues.get_child(name), ptreeOfRequired.get_child(name), parsedPtree.get_child(name), path, list);
    }

    if (parsedPtree.find(name) == parsedPtree.not_found()) {
      if (ptreeOfRequired.get<bool>(name)) {
        list.push_back(path + "." + name);
      }
      else {
        const auto & it = ptreeOfDefaultValues.find(name);
        parsedPtree.put(name, it->second.data());
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
    return m_ConfigIsEnabled;
  }

  void SetConfigFileName(const std::string & fileName)
  {
    m_Config = fileName;
  }

  bool ParseCommandLine(int argc, char** argv)
  {
    try {
      po::parsed_options parsedOptions = po::command_line_parser(argc, argv).options(m_Description).run();
      po::store(parsedOptions, m_Vm);
      po::notify(m_Vm);
    }
    catch (const po::error& e) {
      std::cerr << "An exception occurred while parsing the command line." << std::endl;
      std::cerr << e.what() << endl;
      std::cout << m_Description << std::endl;
      return false;
    }
    if (m_Help) {
      std::cout << m_Description << std::endl;
      return false;
    }

    m_ConfigIsEnabled = !m_Vm["config"].empty();
    return true;
  }

  bool ParseConfigFile()
  {
    try {
      pt::ini_parser::read_ini(m_Config, m_ParsedPtree);
    }
    catch (const pt::ptree_error &e) {
      std::cerr << "An exception occurred while parsing the config file: " << AddQuotes(m_Config) << std::endl;
      std::cout << e.what() << endl;
      return false;
    }

    if (m_ParsedPtree.find(m_NameOfGroup) == m_ParsedPtree.not_found()) {
      std::cerr << "The group " << AddQuotes(m_NameOfGroup) << " is not found in the config file: " << AddQuotes(m_Config) << std::endl;
      return false;
    }

    m_ParsedPtree = m_ParsedPtree.get_child(m_NameOfGroup);

    pt::ptree tree;
    for (const auto & it : m_ParsedPtree) {
      tree.put(it.first, (it.second).data());
    }
    m_ParsedPtree.swap(tree);

    // check parsed ptree
    std::vector<std::string> list;
    checkParsedTree(m_PtreeOfDefaultValues, m_PtreeOfRequired, m_ParsedPtree, m_NameOfGroup, list);

    // print parsed tree
    PrintConfig();

    if (list.size() > 0) {
      std::cerr << "The required keys are not found in the config file: " << AddQuotes(m_Config) << std::endl;
      for (const auto & path : list) {
        std::cout << AddQuotes(path) << std::endl;
      }
      return false;
    }

    // clear map of variables from command line
    m_Vm.clear();

    return true;
  }

  void PrintConfig()
  {
    std::cout << std::endl;
    std::cout << "Config data for group " << AddQuotes(m_NameOfGroup) << std::endl;
    printTree(m_ParsedPtree, std::cout);
    std::cout << std::endl;
  }

  template <typename T>
  T GetDefaultValue(const std::string & str) const
  {
    return m_PtreeOfDefaultValues.get<T>(str);
  }

  template <typename T>
  T Get(const std::string & path) const
  {
    T value;

    if (m_ConfigIsEnabled) {
      try {
        value = m_ParsedPtree.get<T>(path);
      }
      catch (const pt::ptree_error &e) {
        std::cerr << "An exception occurred while getting value from ptree." << std::endl;
        std::cerr << e.what() << std::endl;

        const auto & it = m_ParsedPtree.find(path);
        if (it != m_ParsedPtree.not_found()) {
          std::cerr << AddQuotes(it->first) << " " << it->second.data() << std::endl;
        }

        throw;
      }
    }
    else {
      value = m_Vm[path].as<T>();
    }

    return value;
  }

  template<typename T>
  std::vector<T> GetAsVector(const std::string & path)
  {
    std::stringstream stream(Get<std::string>(path));
    std::string item;

    std::vector<T> vector;

    while (std::getline(stream, item, m_Dlm)) {
      // skip empty spaces ''
      if (item == "")
        continue;

      // get value
      try {
        vector.push_back(std::stod(item));
      }
      catch (const std::invalid_argument &e) {
        std::cerr << "An exception occurred while getting vector from ptree." << std::endl;
        std::cerr << e.what() << " " << AddQuotes(item) << std::endl;

        const auto & it = m_ParsedPtree.find(path);
        std::cerr << AddQuotes(it->first) << " " << it->second.data() << std::endl;

        throw;
      }
    }

    return vector;
  }


private:

  bool m_Help;
  bool m_ConfigIsEnabled;
  std::string m_Config;
  std::string m_NameOfGroup;
  char m_Dlm;

  pt::ptree m_ParsedPtree;
  pt::ptree m_PtreeOfRequired;
  pt::ptree m_PtreeOfDefaultValues;
  po::variables_map m_Vm;
  po::options_description m_Description;

protected:
  OptionsBase()
  {
    m_Help = false;
    m_ConfigIsEnabled = false;
    m_Dlm = ' ';

    po::options_description configOptions("Optional config options");
    configOptions.add_options()("config,c", po::value<std::string>(&m_Config), "The path to the config file.");

    po::options_description helpOptions("Optional help options");
    helpOptions.add_options()("help,h", po::bool_switch(&m_Help)->default_value(m_Help), "Display this help message");

    m_Description.add(configOptions).add(helpOptions);
  }

  void SetNameOfGroup(const std::string & str)
  {
    m_NameOfGroup = str;
  }

  void AddToDescription(const po::options_description & options)
  {
    m_Description.add(options);
  }

  template <typename T>
  void Put(const std::string & path, const T & value, const bool & required = true)
  {
    m_PtreeOfRequired.put(path, required);
    m_PtreeOfDefaultValues.put(path, value);
  }

  std::string Path(const std::string & path) const
  {
    return m_NameOfGroup + "." + path;
  }
};
}
