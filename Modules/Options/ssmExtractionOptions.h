#pragma once

#include "ssmOptionsBase.h"

namespace ssm
{
//=========================================================================
// Surface extraction options
//=========================================================================
class ExtractionOptions : public OptionsBase
{
public:

  void SetInputFileName(const std::string & str)
  {
    m_InputFileName = str;
  }

  const std::string & GetInputFileName() const
  {
    return m_InputFileName;
  }

  void SetOutputFileName(const std::string & str)
  {
    m_OutputFileName = str;
  }

  const std::string & GetOutputFileName() const
  {
    return m_OutputFileName;
  }

  std::string GetInputList() const
  {
    return this->Get<std::string>("inplist");
  }

  std::string GetOutputList() const
  {
    return this->Get<std::string>("outlist");
  }

  std::string GetReportFile() const
  {
    return this->Get<std::string>("report");
  }

  double GetSigma() const
  {
    return this->Get<double>("sigma");
  }

  double GetFactor() const
  {
    return this->Get<double>("factor");
  }

  size_t GetNumberOfPoints() const
  {
    return this->Get<size_t>("points");
  }

  size_t GetNumberOfIterations() const
  {
    return this->Get<size_t>("iterations");
  }

  std::string FormatOutput(const std::string & fileName)
  {
    const auto & format = this->Get<std::string>("output");
    try {
      return (boost::format(format) % getBaseNameFromPath(fileName)).str();
    }
    catch (const boost::io::format_error &e) {
      std::cerr << "Could not format string with format " << format << std::endl;
      std::cout << e.what() << std::endl;
      throw;
    }
  };

  ExtractionOptions()
  {
    SetNameOfGroup("EXTRACTION");

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
      ("input,i", po::value<std::string>(&m_InputFileName), "The path to the input image file.")
      ("output,o", po::value<std::string>(&m_OutputFileName), "The path for the output surface file.")
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

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
    this->AddToDescription(reportOptions);
  };

private:
  std::string m_InputFileName;
  std::string m_OutputFileName;
};
}
