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

  ExtractionOptions()
  {
    SetNameOfGroup("EXTRACTION");

    // initialize ptree
    Put<std::string>("inplist", "");
    Put<std::string>("outlist", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");

    Put<double>("sigma", 0, 0);
    Put<size_t>("iterations", 100, 0);
    Put<size_t>("points", 0, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("input,i", po::value<std::string>(), "Path to input image file.")
      ("output,o", po::value<std::string>(), "Path for output surface file.")
      ;

    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("sigma", po::value<double>()->default_value(this->GetDefaultValue<double>("sigma")), "Sigma of the Gaussian kernel for RecursiveGaussianImageFilter to smooth input image")
      ("iterations", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("iterations")), "Number of iterations to adjust point positions for output surface")
      ("points", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("points")), "The number of points in output decimated surface (default value 0, i.e. no decimation)")
      ;

    po::options_description reportOptions("Optional report options");
    reportOptions.add_options()
      ("report,r", po::value<std::string>(), "Output report file")
      ;

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
    this->AddToDescription(reportOptions);
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
  }

  const std::string & GetInputFileName() const { return this->Get<std::string>("input"); }
  const std::string & GetOutputFileName() const { return this->Get<std::string>("output"); }
  std::string GetInputList() const { return this->Get<std::string>("inplist"); }
  std::string GetOutputList() const { return this->Get<std::string>("outlist"); }
  std::string GetReportFileName() const { return this->Get<std::string>("report"); }
  double GetSigma() const { return this->Get<double>("sigma"); }
  double GetFactor() const { return this->Get<double>("factor"); }
  size_t GetNumberOfPoints() const { return this->Get<size_t>("points"); }
  size_t GetNumberOfIterations() const { return this->Get<size_t>("iterations"); }

};
}
