#pragma once

#include <itkRecursiveGaussianImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <vtkMarchingCubes.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDecimatePro.h>

#include "ssmBinaryMask3DMeshSource.h"

namespace ssm
{
template< typename TInputImage, typename TOutputMesh >
BinaryMask3DMeshSource< TInputImage, TOutputMesh>::BinaryMask3DMeshSource()
{
  m_Sigma = 0;
  m_NumberOfIterations = 100;
  m_RelaxationFactor = 0.2;
  m_NumberOfPoints = 0;

  m_ComputeLevelValue = true;
}


template< typename TInputImage, typename TOutputMesh >
void BinaryMask3DMeshSource< TInputImage, TOutputMesh >::SetInput(const TInputImage *image)
{
  this->ProcessObject::SetNthInput( 0, const_cast< InputImageType * >( image ) );
}

template< typename TInputImage, typename TOutputMesh >
typename TInputImage::ConstPointer BinaryMask3DMeshSource<TInputImage, TOutputMesh>::GetInput()
{
  return static_cast<const TInputImage*>(this->ProcessObject::GetInput(0));
}

template< typename TInputImage, typename TOutputMesh >
TOutputMesh* BinaryMask3DMeshSource<TInputImage, TOutputMesh>::GetOutput()
{
  return m_Output;
}

/** Generate the data */
template< typename TInputImage, typename TOutputMesh >
void BinaryMask3DMeshSource< TInputImage, TOutputMesh >::GenerateData()
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
  if (m_NumberOfPoints > 0) {
    m_Reduction = (m_Output->GetNumberOfPoints() - m_NumberOfPoints) / (double) m_Output->GetNumberOfPoints();
    typedef vtkSmartPointer<vtkDecimatePro> DecimatePolyData;
    auto decimate = DecimatePolyData::New();
    decimate->SetInputData(m_Output);
    decimate->SetTargetReduction(m_Reduction);
    decimate->SetPreserveTopology(true);
    decimate->SetSplitting(false);
    decimate->Update();
    m_Output = decimate->GetOutput();
  }

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  auto smoother = SmoothPolyData::New();
  smoother->SetInputData(m_Output);
  smoother->SetNumberOfIterations(m_NumberOfIterations);
  smoother->SetRelaxationFactor(m_RelaxationFactor);
  smoother->Update();

  typedef vtkSmartPointer<vtkPolyDataNormals> PolyDataNormals;
  auto normals = PolyDataNormals::New();
  normals->SetInputData(smoother->GetOutput());
  normals->AutoOrientNormalsOn();
  normals->FlipNormalsOff();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->SplittingOff();
  normals->Update();

  m_Output = normals->GetOutput();
}

/** Print report */
template< typename TInputImage, typename TOutputMesh >
void BinaryMask3DMeshSource< TInputImage, TOutputMesh >::PrintReport() const
{

}
}
