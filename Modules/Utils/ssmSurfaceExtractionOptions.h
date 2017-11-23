#pragma once

#include "ssmBaseOptions.h"

namespace ssm
{
//=========================================================================
// Surface extraction options
//=========================================================================
class SurfaceExtractionOptions : public OptionsBase
{
public:

  void SetInputFileName(const std::string & str)
  {
    inputFileName = str;
  }

  const std::string & GetInputFileName() const
  {
    return inputFileName;
  }

  void SetOutputFileName(const std::string & str)
  {
    outputFileName = str;
  }

  const std::string & GetOutputFileName() const
  {
    return outputFileName;
  }

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
    const std::string name = "report";
    if (configIsEnabled)
      return parsedPtree.get<std::string>(name);
    else
      return vm[name].as<std::string>();
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

  SurfaceExtractionOptions()
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
}
