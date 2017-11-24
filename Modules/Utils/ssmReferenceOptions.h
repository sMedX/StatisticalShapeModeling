#pragma once

#include "ssmBaseOptions.h"

namespace ssm
{
//=========================================================================
// Surface reference options
//=========================================================================
class ReferenceOptions : public OptionsBase
{
public:

  std::string GetInputFileName() const
  {
    const std::string name = "input";
    if (this->ConfigIsEnabled())
      return this->m_ParsedPtree.get<std::string>(name);
    else
      return this->m_Vm[name].as<std::string>();
  }

  std::string GetOutputFileName() const
  {
    const std::string name = "output";
    if (this->ConfigIsEnabled())
      return this->m_ParsedPtree.get<std::string>(name);
    else
      return this->m_Vm[name].as<std::string>();
  }

  std::string GetReportFileName() const
  {
    const std::string name = "report";
    if (this->ConfigIsEnabled())
      return this->m_ParsedPtree.get<std::string>(name);
    else
      return this->m_Vm[name].as<std::string>();
  }

  double GetSigma() const
  {
    const std::string name = "sigma";
    if (this->ConfigIsEnabled())
      return this->Get<double>(name);
    else
      return this->m_Vm[name].as<double>();
  }

  double GetFactor() const
  {
    const std::string name = "factor";
    if (this->ConfigIsEnabled())
      return this->Get<double>(name);
    else
      return this->m_Vm[name].as<double>();
  }

  size_t GetNumberOfPoints() const
  {
    const std::string name = "points";
    if (this->ConfigIsEnabled())
      return Get<size_t>(name);
    else
      return this->m_Vm[name].as<size_t>();
  }

  size_t GetNumberOfIterations() const
  {
    const std::string name = "iterations";
    if (this->ConfigIsEnabled())
      return this->Get<size_t>(name);
    else
      return this->m_Vm[name].as<size_t>();
  }

  ReferenceOptions()
  {
    SetNameOfGroup("REFERENCE");

    // initialize ptree
    Put<std::string>("input", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");

    Put<double>("sigma", 0, 0);
    Put<double>("factor", 0.2, 0);
    Put<size_t>("iterations", 100, 0);
    Put<size_t>("points", 0, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("input,i", po::value<std::string>(), "The path to the input image file.")
      ("output,o", po::value<std::string>(), "The path for the output surface file.")
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

    m_Description.add(mandatoryOptions).add(inputOptions).add(reportOptions);
  }
};
}
