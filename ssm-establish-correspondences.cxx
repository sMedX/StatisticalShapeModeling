﻿#include <boost/program_options.hpp>
#include <itkLowRankGPModelBuilder.h>
#include <itkStandardMeshRepresenter.h>
#include <statismo-build-gp-model-kernels.h>
#include <utils/statismo-build-models-utils.h>

#include "utils/ssmTypes.h"
#include "utils/io.h"
#include "utils/itkCommandLineArgumentParser.h"

#include "ssm/ssmPointSetToImageMetrics.h"
#include "ssm/ssmShapeModelRegistrationMethod.h"
#include "ssm/ssmSurfaceToLevelSetImageFilter.h"

struct options
{
  std::vector<double> parameters;
  std::vector<double> regularization;
  double scale = 100;
  unsigned int iterations = 500;
  unsigned int degree;
};

typedef itk::StatisticalModel<MeshType> StatisticalModelType;
StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions);
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer model, options & opt);


int main(int argc, char** argv)
{
  options opt;

  itk::CommandLineArgumentParser::Pointer parser = itk::CommandLineArgumentParser::New();
  parser->SetCommandLineArguments(argc, argv);

  std::string listFile;
  parser->GetCommandLineArgument("-list", listFile);

  std::string outputFile;
  parser->GetCommandLineArgument("-output", outputFile);

  std::string referenceFileName;
  parser->GetCommandLineArgument("-reference", referenceFileName);

  parser->GetCommandLineArgument("-iterations", opt.iterations);
  parser->GetCommandLineArgument("-regularization", opt.regularization);
  parser->GetCommandLineArgument("-parameters", opt.parameters);
  parser->GetCommandLineArgument("-scale", opt.scale);
  parser->GetCommandLineArgument("-degree", opt.degree);

  unsigned int numberOfStages = 1;
  parser->GetCommandLineArgument("-stages", numberOfStages);

  unsigned int numberOfComponents = 100;
  parser->GetCommandLineArgument("-components", numberOfComponents);

  std::cout << std::endl;
  std::cout << " shape model to image registration" << std::endl;
  std::cout << "   output surface file " << outputFile << std::endl;
  std::cout << "      number of stages " << numberOfStages << std::endl;
  std::cout << "  number of iterations " << opt.iterations << std::endl;
  std::cout << "                degree " << opt.degree << std::endl;

  for (int n = opt.regularization.size(); n < opt.parameters.size(); ++n) {
    opt.regularization.push_back(opt.regularization.back());
  }

  std::cout << "        regularization ";
  for (int n = 0; n < opt.regularization.size(); ++n) {
    std::cout << opt.regularization[n] << " ";
  }
  std::cout << std::endl;

  std::cout << "            parameters ";
  for (int n = 0; n < opt.parameters.size(); ++n) {
    std::cout << opt.parameters[n] << " ";
  }
  std::cout << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read the reference surface
  MeshType::Pointer reference = MeshType::New();
  if (!readMesh<MeshType>(reference, referenceFileName)) {
    return EXIT_FAILURE;
  }
  std::cout << "reference " << referenceFileName << std::endl;
  std::cout << "number of cells " << reference->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << reference->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;

  //----------------------------------------------------------------------------
  // read list of files
  StringList listOfFiles;
  try {
    listOfFiles = getFileList(listFile);
  }
  catch (ifstream::failure & e) {
    cerr << "could not read the data-list: " << e.what() << endl;
    return EXIT_FAILURE;
  }

  std::vector<MeshType::Pointer> vectorOfSurfaces;
  std::vector<std::string> vectorOfFiles;

  for (const auto & fileName : listOfFiles) {
    MeshType::Pointer surface = MeshType::New();

    if (!readMesh<MeshType>(surface, fileName)) {
      return EXIT_FAILURE;
    }
    vectorOfSurfaces.push_back(surface);
    vectorOfFiles.push_back(fileName);

    std::cout << "surface " << fileName << std::endl;
    std::cout << "number of cells   " << surface->GetNumberOfCells() << std::endl;
    std::cout << "number of points  " << surface->GetNumberOfPoints() << std::endl;
    std::cout << std::endl;
  }

  //----------------------------------------------------------------------------
  //
  for (unsigned int stage = 0; stage < numberOfStages; ++stage) {
    std::string strInfo = "---------- stage (" + std::to_string(stage + 1) + " / " + std::to_string(numberOfStages);

    // build GP model for the reference surface
    StatisticalModelType::Pointer model = BuildGPModel(reference, opt.parameters[0], opt.scale, numberOfComponents);

    // initialize reference
    for (auto & iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
      iter.Value().Fill(0);
    }

    // preform GP model to surface registration
    for (size_t count = 0; count < vectorOfSurfaces.size(); ++count) {
      MeshType::Pointer surface = vectorOfSurfaces[count];
      MeshType::Pointer output = shapeModelToSurfaceRegistration(surface, model, opt);

      // compute the new reference surface
      for (auto & iter = reference->GetPoints()->Begin(); iter != reference->GetPoints()->End(); ++iter) {
        MeshType::PointType point = output->GetPoint(iter.Index());
        for (size_t d = 0; d < MeshType::PointDimension; ++d) {
          iter.Value()[d] += point[d] / vectorOfSurfaces.size();
        }
      }

      // compute metrics and write surface to file
      if (stage == numberOfStages - 1) {
        // compute level set image
        typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
        SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
        levelset->SetMargin(0.10);
        levelset->SetSpacing(1.0);
        levelset->SetInput(surface);
        try {
          levelset->Update();
        }
        catch (itk::ExceptionObject& excep) {
          std::cerr << excep << std::endl;
          return EXIT_FAILURE;
        }

        // compute metrics
        typedef itk::PointSet<MeshType::PixelType, MeshType::PointDimension> PointSetType;
        PointSetType::Pointer points = PointSetType::New();
        points->SetPoints(surface->GetPoints());

        typedef ssm::PointSetToImageMetrics<PointSetType, FloatImageType> PointSetToImageMetricsType;
        PointSetToImageMetricsType::Pointer metrics = PointSetToImageMetricsType::New();
        metrics->SetFixedPointSet(points);
        metrics->SetMovingImage(levelset->GetOutput());
        metrics->Compute();
        metrics->PrintReport(std::cout);

        // print report to *.csv file
        if (parser->ArgumentExists("-report")) {
          std::string fileName;
          parser->GetCommandLineArgument("-report", fileName);
          std::cout << "print report to file " << fileName << std::endl;
          metrics->PrintReportToFile(fileName, getBaseNameFromPath(vectorOfFiles[count]));
        }

        // write output surface to the file
        typedef boost::filesystem::path fp;
        fp path = fp(outputFile).parent_path() / fp(fp(vectorOfFiles[count]).stem().string() + "-" + fp(outputFile).filename().string());
        std::string fileName = path.string();
        std::cout << "write surface to file " << fileName << std::endl;
        if (!writeMesh<MeshType>(output, fileName)) {
          EXIT_FAILURE;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
//==============================================================================
// Shape model to surface registration
MeshType::Pointer shapeModelToSurfaceRegistration(MeshType::Pointer surface, StatisticalModelType::Pointer initialModel, options & opt)
{
  // compute level set image
  typedef ssm::SurfaceToLevelSetImageFilter<MeshType, FloatImageType> SurfaceToLevelSetImageFilter;
  SurfaceToLevelSetImageFilter::Pointer levelset = SurfaceToLevelSetImageFilter::New();
  levelset->SetMargin(0.10);
  levelset->SetSpacing(1.0);
  levelset->SetInput(surface);
  try {
    levelset->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  // shape model to image registration
  StatisticalModelType::Pointer model = initialModel;
  size_t numberOfStages = opt.parameters.size();

  for (size_t stage = 0; stage < numberOfStages; ++stage) {
    std::cout << "---------- stage (" << stage + 1 << " / " << numberOfStages << ") ----------" << std::endl;

    // if stage > 0 update GP model
    if (stage > 0) {
      model = BuildGPModel(surface, opt.parameters[stage], opt.scale, model->GetNumberOfPrincipalComponents());
    }

    std::cout << "built Gaussian process model" << std::endl;
    std::cout << "number of components " << model->GetNumberOfPrincipalComponents() << std::endl;
    std::cout << "number of points     " << model->GetRepresenter()->GetReference()->GetNumberOfPoints() << std::endl;
    std::cout << "parameters           " << opt.parameters[stage] << ", " << opt.scale << std::endl;
    std::cout << std::endl;

    // perform registration
    typedef ssm::ShapeModelRegistrationMethod<StatisticalModelType, MeshType> ShapeModelRegistrationMethod;
    ShapeModelRegistrationMethod::Pointer shapeModelToSurfaceRegistration;
    shapeModelToSurfaceRegistration = ShapeModelRegistrationMethod::New();
    shapeModelToSurfaceRegistration->SetShapeModel(model);
    shapeModelToSurfaceRegistration->SetLevelSetImage(levelset->GetOutput());
    shapeModelToSurfaceRegistration->SetNumberOfIterations(opt.iterations);
    shapeModelToSurfaceRegistration->SetRegularizationParameter(opt.regularization[stage]);
    shapeModelToSurfaceRegistration->SetDegree(opt.degree);
    try {
      shapeModelToSurfaceRegistration->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      throw;
    }
    shapeModelToSurfaceRegistration->PrintReport(std::cout);

    surface = const_cast<MeshType*> (shapeModelToSurfaceRegistration->GetOutput());
  }

  return surface;
}

//==============================================================================
// Build Gaussian process model
StatisticalModelType::Pointer BuildGPModel(MeshType::Pointer surface, double parameters, double scale, int numberOfBasisFunctions)
{
  // create kernel
  typedef DataTypeShape::PointType PointType;
  auto gaussianKernel = new GaussianKernel<PointType>(parameters);

  typedef std::shared_ptr<const statismo::ScalarValuedKernel<PointType>> MatrixPointerType;
  MatrixPointerType kernel(gaussianKernel);

  typedef std::shared_ptr<statismo::MatrixValuedKernel<PointType>> KernelPointerType;
  KernelPointerType unscaledKernel(new statismo::UncorrelatedMatrixValuedKernel<PointType>(kernel.get(), Dimension));
  KernelPointerType modelBuildingKernel(new statismo::ScaledKernel<PointType>(unscaledKernel.get(), scale));

  typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
  typedef RepresenterType::DatasetPointerType DatasetPointerType;
  RepresenterType::Pointer representer = RepresenterType::New();
  representer->SetReference(surface);

  // build model
  typedef itk::LowRankGPModelBuilder<MeshType> ModelBuilderType;
  ModelBuilderType::Pointer modelBuilder = ModelBuilderType::New();
  modelBuilder->SetRepresenter(representer);

  // build and save model to file
  typedef itk::StatisticalModel<MeshType> StatisticalModelType;
  StatisticalModelType::Pointer model;
  try {
    model = modelBuilder->BuildNewModel(representer->IdentitySample(), *modelBuildingKernel.get(), numberOfBasisFunctions);
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    throw;
  }

  return model;
}

