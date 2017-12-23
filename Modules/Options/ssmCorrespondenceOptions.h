#pragma once

#include "ssmOptionsBase.h"

namespace ssm
{
//=========================================================================
// Establish correspondence options
//=========================================================================
class CorrespondenceOptions : public OptionsBase
{
public:

  std::string GetInputList() const
  {
    return this->Get<std::string>("inplist");
  }

  std::string GetOutputList() const
  {
    return this->Get<std::string>("outlist");
  }

  std::string GetReportFileName() const
  {
    return this->Get<std::string>("report");
  }

  std::string GetReferenceFileName() const
  {
    return this->Get<std::string>("reference");
  }

  size_t GetNumberOfStages() const
  {
    return this->Get<size_t>("stages");
  }

  size_t GetNumberOfIterations() const
  {
    return this->Get<size_t>("iterations");
  }

  size_t GetTransform() const
  {
    return this->Get<size_t>("transform");
  }

  double GetGPModelScale() const
  {
    return this->Get<double>("gpmodel.scale");
  }

  std::vector<double> GetGPModelParameters() const
  {
    return m_Parameters;
  }

  std::vector<size_t> GetGPModelNumberOfComponents() const
  {
    return m_Components;
  }

  std::vector<double> GetGPModelRegularization() const
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
  }

  bool ParseOptions(int argc, char** argv)
  {
    if (!OptionsBase::ParseOptions(argc, argv)) {
      return false;
    }

    if (!CheckOptions()) {
      return false;
    }

    return true;
  }

  bool CheckOptions()
  {
    try {
      m_Parameters = this->GetAsVector<double>("gpmodel.parameters");
      m_Components = this->GetAsVector<size_t>("gpmodel.components");
      m_Regularization = this->GetAsVector<double>("gpmodel.regularization");
    }
    catch (...) {
      return false;
    }

    size_t numberOfStages = m_Parameters.size();

    if (numberOfStages == 0 || m_Components.size() == 0 || m_Regularization.size() == 0) {
      std::cerr << "parameters, components and regularization factors must be specified" << std::endl;
      return false;
    }

    for (size_t n = m_Components.size(); n < numberOfStages; ++n) {
      m_Components.push_back(m_Components.back());
    }

    for (size_t n = m_Regularization.size(); n < numberOfStages; ++n) {
      m_Regularization.push_back(m_Regularization.back());
    }

    return checkFileName(GetReportFileName());
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
    Put<size_t>("iterations", 1000);

    Put<double>("gpmodel.scale", 50);
    Put<double>("gpmodel.parameters", 50);
    Put<size_t>("gpmodel.components", 100);
    Put<double>("gpmodel.regularization", 0.10);
  };

private:
  std::string m_InputFileName;
  std::string m_OutputFileName;

  std::vector<double> m_Regularization;
  std::vector<double> m_Parameters;
  std::vector<size_t> m_Components;
};
}
