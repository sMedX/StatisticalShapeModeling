#pragma once

#include "ssmBaseOptions.h"

namespace ssm
{
//=========================================================================
// Surface extraction options
//=========================================================================
class AlignmentOptions : public OptionsBase
{
public:

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
    return parsedPtree.get<std::string>("report");
  }

  std::string GetReferenceFileName() const
  {
    return parsedPtree.get<std::string>("reference");
  }

  size_t GetNumberOfStages() const
  {
    return parsedPtree.get<size_t>("stages");
  }

  size_t GetNumberOfIterations() const
  {
    return parsedPtree.get<size_t>("iterations");
  }

  size_t GetTransform() const
  {
    return parsedPtree.get<size_t>("transform");
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

  AlignmentOptions()
  {
    SetNameOfGroup("ALIGNMENT");

    // initialize ptree
    Put<std::string>("inplist", "");
    Put<std::string>("outlist", "");
    Put<std::string>("output", "");
    Put<std::string>("report", "");
    Put<std::string>("reference", "");

    Put<size_t>("transform", 2);
    Put<size_t>("stages", 3);
    Put<size_t>("iterations", 500);
  };

private:
  std::string inputFileName;
  std::string outputFileName;
};
}
