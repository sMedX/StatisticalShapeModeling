#pragma once

#include "ssmOptionsBase.h"

namespace ssm
{
//=========================================================================
// Cross validation options
//=========================================================================
class ModelCrossValidationOptions : public OptionsBase
{
public:
  std::string GetInputList() const {return this->Get<std::string>("inplist");}
  std::string GetReportFileName()const {return this->Get<std::string>("cvtest.report");}
  int GetNumberOfComponents() const {return this->Get<int>("cvtest.components");}
  bool GetWrite() const {return this->Get<bool>("cvtest.write");}
  
  bool ParseOptions(int argc, char** argv)
  {
    if (!OptionsBase::ParseOptions(argc, argv)) {
      return false;
    }

    return checkFileName(GetReportFileName());
  }

  ModelCrossValidationOptions()
  {
    SetNameOfGroup("MODELQUALITY");

    // initialize ptree
    Put<std::string>("inplist", "");
    Put<std::string>("cvtest.report", "", 0);
    Put<int>("cvtest.components", -1, 0);
    Put<bool>("cvtest.write", false, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("inplist,i", po::value<std::string>(), "Input file with a list of files of input surfaces.")
      ;

    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("report,r", po::value<std::string>(), "Output report file")
      ("components", po::value<int>()->default_value(this->GetDefaultValue<int>("cvtest.components")), "Number of components to reduce shape model for cross validation test.")
      ("write", po::value<bool>()->default_value(this->GetDefaultValue<bool>("cvtest.write")), "Write surfaces.")
      ;

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
  }
};
}
