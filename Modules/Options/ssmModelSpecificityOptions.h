#pragma once

#include "ssmOptionsBase.h"

namespace ssm
{
//=========================================================================
// Model specificity options
//=========================================================================
class ModelSpecificityOptions : public OptionsBase
{
public:
  std::string GetInputList() const { return this->Get<std::string>("inplist"); }
  std::string GetModelFileName() const {return this->Get<std::string>("model");}
  std::string GetReportFileName() const {return this->Get<std::string>("specificity.report");}
  size_t GetNumberOfSamples() const {return this->Get<size_t>("specificity.samples");}

  bool ParseOptions(int argc, char** argv)
  {
    if (!OptionsBase::ParseOptions(argc, argv)) {
      return false;
    }

    return checkFileName(GetReportFileName());
  }

  ModelSpecificityOptions()
  {
    this->SetNameOfGroup("MODELQUALITY");

    // initialize ptree
    this->Put<std::string>("inplist","");
    this->Put<std::string>("model", "");
    this->Put<std::string>("specificity.report", "", 0);
    this->Put<size_t>("specificity.samples", 1000, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("inplist,i", po::value<std::string>(), "Input file with a list of files of input surfaces.")
      ("model,m", po::value<std::string>(), "Input model file.")
      ;

    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("report,r", po::value<std::string>(), "Output report file")
      ("samples", po::value<size_t>()->default_value(this->GetDefaultValue<size_t>("specificity.samples")), "The number of random samples to compute specificity.")
      ;

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
  }

private:

};
}
