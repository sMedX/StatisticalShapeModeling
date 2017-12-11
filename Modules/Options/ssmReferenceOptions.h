#pragma once

#include "ssmOptionsBase.h"

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
    return this->Get<std::string>("input");
  }

  std::string GetOutputFileName() const
  {
    return this->Get<std::string>("output");
  }

  std::string GetReportFileName() const
  {
    return this->Get<std::string>("report");
  }

  double GetSigma() const
  {
    return this->Get<double>("sigma");
  }

  double GetLevelValue() const
  {
    return this->Get<double>("level");
  }

  size_t GetNumberOfPoints() const
  {
    return this->Get<size_t>("points");
  }

  size_t GetDecimation() const
  {
    return this->Get<size_t>("decimation");
  }

  size_t GetSmoothing() const
  {
    return this->Get<size_t>("smoothing");
  }

  double GetRelaxationFactor() const
  {
    return this->Get<double>("factor");
  }

  size_t GetNumberOfIterations() const
  {
    return this->Get<size_t>("iterations");
  }

  ReferenceOptions()
  {
    SetNameOfGroup("REFERENCE");

    // initialize ptree
    Put<std::string>("input", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");

    Put<double>("sigma", 0, 0);
    Put<double>("level", 0, 0);

    Put<size_t>("smoothing", 0, 0);
    Put<double>("factor", 0.2, 0);
    Put<size_t>("iterations", 100, 0);

    Put<size_t>("decimation", 0, 0);
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
      ("level", po::value<double>()->default_value(this->GetDefaultValue<double>("level")), "The level value to extract surface from input level set image.")
      ("smoothing", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("smoothing")), "The method of surface smoothing \n 0 --- None \n 1 --- vtkWindowedSincPolyDataFilter \n 2 --- vtkSmoothPolyDataFilter.")
      ("factor", po::value<double>()->default_value(this->GetDefaultValue<double>("factor")), "The relaxation factor for Laplacian smoothing.")
      ("iterations", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("iterations")), "The number of iterations.")
      ("decimation", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("decimation")), "The method of surface decimation \n 0 --- None \n 1 --- vtkQuadricDecimation \n 2 --- vtkDecimatePro.")
      ("points", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("points")), "The number of points in output surface.")
      ;

    po::options_description reportOptions("Optional report options");
    reportOptions.add_options()
      ("report,r", po::value<std::string>(), "The path for the file to print report.")
      ;

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
    this->AddToDescription(reportOptions);
  }
};
}
