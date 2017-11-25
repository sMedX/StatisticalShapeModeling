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
      name.resize(16);
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

    // check parsed ptree
    std::vector<std::string> listOfKeys;
    checkParsedTree(m_PtreeOfRequireds, m_ParsedPtree, m_NameOfGroup, listOfKeys);

    if (listOfKeys.size() > 0) {
      std::cerr << "The keys are not found in the config file: " << AddQuotes(m_Config) << std::endl;
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
    std::cout << "Config data for group " << AddQuotes(m_NameOfGroup) << std::endl;
    printTree(m_ParsedPtree, std::cout, 0);
    std::cout << std::endl;
  };

  template <typename T>
  T GetDefaultValue(const std::string & str) const
  {
    return m_PtreeOfDefaultValues.get<T>(str);
  };

  template <typename T>
  T Get(const std::string & name) const
  {
    if (m_ConfigIsEnabled)
      return m_ParsedPtree.get<T>(name);
    else
      return m_Vm[name].as<T>();
  };

  template<typename T>
  bool GetAsVector(const std::string & key, std::vector<T> & v)
  {
    std::stringstream stream(Get<std::string>(key));
    std::string item;

    while (std::getline(stream, item, m_Dlm)) {
      try {
        v.push_back(std::stod(item));
      }
      catch (...) {
        std::cerr << "Error while parsing string in the method GetAsVector." << std::endl;
        std::cout << "string: " << AddQuotes(Get<std::string>(key)) << std::endl;
        return false;
      }
    }

    return true;
  }


private:

  bool m_Help;
  bool m_ConfigIsEnabled;
  std::string m_Config;
  std::string m_NameOfGroup;
  char m_Dlm;

  pt::ptree m_ParsedPtree;
  pt::ptree m_PtreeOfRequireds;
  pt::ptree m_PtreeOfDefaultValues;
  po::variables_map m_Vm;

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

  template <typename T>
  void Put(const std::string & str, const T & value, const bool & required = true)
  {
    m_PtreeOfRequireds.put(str, required);
    m_PtreeOfDefaultValues.put(str, value);
    if (!required) {
      m_ParsedPtree.put(Path(str), value);
    }
  };

  std::string Path(const std::string & str) const
  {
    return m_NameOfGroup + "." + str;
  }

  po::options_description m_Description;
};
}
