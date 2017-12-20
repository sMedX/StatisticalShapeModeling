#pragma once

#include <itkRecursiveGaussianImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>

#include "ssmImage3DMeshSource.h"

namespace ssm
{
template< typename TInputImage, typename TOutputMesh >
Image3DMeshSource< TInputImage, TOutputMesh>::Image3DMeshSource()
{
  m_Sigma = 0;

  m_Smoothing = Smoothing::WindowedSinc;
  m_NumberOfIterations = 100;
  m_RelaxationFactor = 0.2;
  m_FeatureAngle = 60.0;
  m_PassBand = 0.001;

  m_Decimation = Decimation::None;
  m_NumberOfPoints = 0;

  m_ComputeLevelValue = true;
}


template< typename TInputImage, typename TOutputMesh >
void Image3DMeshSource< TInputImage, TOutputMesh >::SetInput(const TInputImage *image)
{
  this->ProcessObject::SetNthInput( 0, const_cast< InputImageType * >( image ) );
}

template< typename TInputImage, typename TOutputMesh >
typename TInputImage::ConstPointer Image3DMeshSource<TInputImage, TOutputMesh>::GetInput()
{
  return static_cast<const TInputImage*>(this->ProcessObject::GetInput(0));
}

template< typename TInputImage, typename TOutputMesh >
TOutputMesh* Image3DMeshSource<TInputImage, TOutputMesh>::GetOutput()
{
  return m_Output;
}

/** Generate the data */
template< typename TInputImage, typename TOutputMesh >
void Image3DMeshSource< TInputImage, TOutputMesh >::GenerateData()
{
// compute the minimum and the maximum intensity values of label
  if (m_ComputeLevelValue) {
    typedef itk::MinimumMaximumImageCalculator <InputImageType> MinimumMaximumImageCalculatorType;
    auto labelValues = MinimumMaximumImageCalculatorType::New();
    labelValues->SetImage(this->GetInput());
    labelValues->Compute();
    m_LevelValue = 0.5 * (labelValues->GetMinimum() + labelValues->GetMaximum());
  }

  // smoothing
  if (std::abs(m_Sigma) < itk::NumericTraits<double>::epsilon()) {
    m_Sigma = this->GetInput()->GetSpacing().GetVnlVector().max_value();
  }

  typedef itk::RecursiveGaussianImageFilter<InputImageType, FloatImageType> RecursiveGaussianImageFilterType;
  auto gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(this->GetInput());
  gaussian->SetSigma(m_Sigma);

  // fill holes after smoothing
  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  auto fillholes = GrayscaleFillholeImageFilterType::New();
  fillholes->SetInput(gaussian->GetOutput());
  fillholes->SetFullyConnected(true);
  fillholes->Update();

  // convert ITK image to VTK image
  typedef itk::ImageToVTKImageFilter<FloatImageType> ConvertorType;
  ConvertorType::Pointer convertor = ConvertorType::New();
  convertor->SetInput(fillholes->GetOutput());
  convertor->Update();

  // extract surface
  typedef vtkSmartPointer<vtkMarchingCubes> MarchingCubes;
  auto mcubes = MarchingCubes::New();
  mcubes->SetInputData(convertor->GetOutput());
  mcubes->SetValue(0, m_LevelValue);
  mcubes->Update();
  m_Output = mcubes->GetOutput();

  // decimate surface
  this->SurfaceDecimation();

  // smoothing surface
  this->SurfaceSmoothing();

  // compute normals
  typedef vtkSmartPointer<vtkPolyDataNormals> PolyDataNormals;
  auto normals = PolyDataNormals::New();
  normals->SetInputData(m_Output);
  normals->AutoOrientNormalsOn();
  normals->FlipNormalsOff();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->SplittingOff();
  normals->Update();

  m_Output = normals->GetOutput();
}

/** Decimate surface */
template< typename TInputImage, typename TOutputMesh >
void Image3DMeshSource< TInputImage, TOutputMesh >::SurfaceDecimation()
{
  // decimation
  if (m_Decimation==Decimation::None || m_NumberOfPoints == 0) {
    return;
  }

  while (m_Output->GetNumberOfPoints() > m_NumberOfPoints) {
    m_Reduction = (m_Output->GetNumberOfPoints() - m_NumberOfPoints) / (double)m_Output->GetNumberOfPoints();

    switch (m_Decimation) {
    case Decimation::QuadricDecimation: {
      typedef vtkSmartPointer<vtkQuadricDecimation> Decimation;
      auto decimate = Decimation::New();
      decimate->SetInputData(m_Output);
      decimate->SetTargetReduction(m_Reduction);
      decimate->Update();
      m_Output = decimate->GetOutput();
      break;
    }
    case Decimation::DecimatePro: {
      typedef vtkSmartPointer<vtkDecimatePro> Decimate;
      auto decimate = Decimate::New();
      decimate->SetInputData(m_Output);
      decimate->SetSplitting(false);
      decimate->SetErrorIsAbsolute(5);
      decimate->SetFeatureAngle(m_FeatureAngle);
      decimate->SetPreserveTopology(true);
      decimate->SetBoundaryVertexDeletion(false);
      decimate->SetDegree(10); // std-value is 25!
      decimate->SetTargetReduction(m_Reduction);
      decimate->SetMaximumError(0.002);
      decimate->Update();
      m_Output = decimate->GetOutput();
      break;
    }
    }
  }
}

/** Decimate surface */
template< typename TInputImage, typename TOutputMesh >
void Image3DMeshSource< TInputImage, TOutputMesh >::SurfaceSmoothing()
{
  // smoothing
  switch (m_Smoothing) {
  case Smoothing::Laplacian: {
    typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
    auto smoother = SmoothPolyData::New();
    smoother->SetInputData(m_Output);
    smoother->SetNumberOfIterations(m_NumberOfIterations);
    smoother->SetRelaxationFactor(m_RelaxationFactor);
    smoother->SetFeatureAngle(m_FeatureAngle);
    smoother->SetConvergence(0);
    smoother->SetBoundarySmoothing(false);
    smoother->SetFeatureEdgeSmoothing(false);
    smoother->Update();
    m_Output = smoother->GetOutput();
    break;
  }
  case Smoothing::WindowedSinc: {
    typedef vtkSmartPointer<vtkWindowedSincPolyDataFilter> SmoothPolyData;
    auto smoother = SmoothPolyData::New();
    smoother->SetInputData(m_Output);
    smoother->SetNumberOfIterations(m_NumberOfIterations);
    smoother->SetFeatureEdgeSmoothing(false);
    smoother->SetFeatureAngle(m_FeatureAngle);
    smoother->SetPassBand(m_PassBand);
    smoother->SetNonManifoldSmoothing(true);
    smoother->SetNormalizeCoordinates(true);
    smoother->SetBoundarySmoothing(false);
    smoother->Update();
    m_Output = smoother->GetOutput();
    break;
  }
  }

}

/** Print report */
template< typename TInputImage, typename TOutputMesh >
void Image3DMeshSource< TInputImage, TOutputMesh >::PrintReport() const
{

}
}
