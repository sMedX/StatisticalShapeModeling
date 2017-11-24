#pragma once

#include "ssmBaseOptions.h"

namespace ssm
{
//=========================================================================
// Surface extraction options
//=========================================================================
class CorrespondenceOptions : public OptionsBase
{
public:

  std::string GetInputList() const
  {
    return m_ParsedPtree.get<std::string>("inplist");
  }

  std::string GetOutputList() const
  {
    return m_ParsedPtree.get<std::string>("outlist");
  }

  std::string GetReportFileName() const
  {
    return m_ParsedPtree.get<std::string>("report");
  }

  std::string GetReferenceFileName() const
  {
    return m_ParsedPtree.get<std::string>("reference");
  }

  size_t GetNumberOfStages() const
  {
    return m_ParsedPtree.get<size_t>("stages");
  }

  size_t GetNumberOfIterations() const
  {
    return m_ParsedPtree.get<size_t>("iterations");
  }

  size_t GetTransform() const
  {
    return m_ParsedPtree.get<size_t>("transform");
  }

  double GetScale() const
  {
    return m_ParsedPtree.get<double>("scale");
  }

  std::vector<double> GetParameters() const
  {
    return m_Parameters;
  }

  std::vector<size_t> GetNumberOfComponents() const
  {
    return m_Components;
  }

  std::vector<double> GetRegularization() const
  {
    return m_Regularization;
  }

  std::string FormatOutput(const std::string & fileName)
  {
    const auto & format = Get<std::string>("output");
    try {
      return (boost::format(format) % getBaseNameFromPath(fileName)).str();
    }
    catch (const boost::io::format_error &e) {
      std::cerr << "Could not format string with format " << format << std::endl;
      std::cout << e.what() << std::endl;
      throw;
    }
  };

  bool ReadConfigFile()
  {
    if (!OptionsBase::ReadConfigFile()) {
      return false;
    }

    if (!this->GetAsVector<double>("regularization", m_Regularization)) {
      return false;
    }

    if (!this->GetAsVector<double>("parameters", m_Parameters) ||
        !this->GetAsVector<size_t>("components", m_Components) ) {
      return false;
    }

    if (m_Parameters.size() == 0 || m_Components.size() == 0 || m_Regularization.size() == 0) {
      std::cout << "parameters, components and regularization factors must be specified" << std::endl;
      return EXIT_FAILURE;
    }

    for (size_t n = m_Components.size(); n < m_Parameters.size(); ++n) {
      m_Components.push_back(m_Components.back());
    }

    for (size_t n = m_Regularization.size(); n < m_Parameters.size(); ++n) {
      m_Regularization.push_back(m_Regularization.back());
    }


    return true;
  }

  CorrespondenceOptions()
  {
    SetNameOfGroup("CORRESPONDENCE");

    // initialize ptree
    Put<std::string>("inplist", "");
    Put<std::string>("outlist", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");
    Put<std::string>("reference", "");

    Put<size_t>("transform", 3);
    Put<size_t>("stages", 1);
    Put<size_t>("iterations", 100);
    Put<double>("regularization", 0.10);

    Put<double>("scale", 50);
    Put<double>("parameters", 50);
    Put<size_t>("components", 100);
  };

private:
  std::string m_InputFileName;
  std::string m_OutputFileName;

  std::vector<double> m_Regularization;

  std::vector<double> m_Parameters;
  std::vector<size_t> m_Components;
};
}
