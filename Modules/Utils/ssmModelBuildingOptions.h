#pragma once

#include "ssmOptionsBase.h"

namespace ssm
{
//=========================================================================
// Surface extraction options
//=========================================================================
class ModelBuildingOptions : public OptionsBase
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

  std::string GetReferenceFileName() const
  {
    return this->Get<std::string>("reference");
  }

  double GetNoise() const
  {
    return this->Get<double>("noise");
  }

  std::string GetAlignmentMode() const
  {
    return this->Get<std::string>("alignment");
  }

  ModelBuildingOptions()
  {
    SetNameOfGroup("MODELBUILDING");

    // initialize ptree
    Put<std::string>("input", "");
    Put<std::string>("output", "");

    Put<std::string>("alignment", "GPA", 0);
    Put<std::string>("reference", "", 0);
    Put<double>("noise", 0, 0);

    // initialize description
    po::options_description mandatoryOptions("Mandatory options");
    mandatoryOptions.add_options()
      ("input,i", po::value<std::string>(), "File containing a list of meshes to build shape model from")
      ("output,o", po::value<std::string>(), "Name of the output file")
      ;

    po::options_description inputOptions("Optional input options");
    inputOptions.add_options()
      ("alignment", po::value<std::string>()->default_value("GPA"), "Specify how the data is aligned: REFERENCE aligns all datasets rigidly to the reference and GPA alignes all datasets to the population mean.")
      ("reference", po::value<std::string>(), "Specify the reference used for model building. This is needed if --alignment is REFERENCE")
      ("noise", po::value<double>()->default_value(0), "Noise variance of the PPCA model")
      ;

    this->AddToDescription(mandatoryOptions);
    this->AddToDescription(inputOptions);
  }

private:

};
}
