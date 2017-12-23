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
  std::string GetInputList() { return this->Get<std::string>("inplist"); }
  std::string GetModelFileName()  {return this->Get<std::string>("specificity.model");}
  std::string GetReportFileName()  {return this->Get<std::string>("specificity.report");}
  size_t GetNumberOfSamples()  {return this->Get<size_t>("specificity.samples");}

  ModelSpecificityOptions()
  {
    this->SetNameOfGroup("MODELQUALITY");

    // initialize ptree
    this->Put<std::string>("inplist","");
    this->Put<std::string>("specificity.model", "");
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
