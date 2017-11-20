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
  m_Sigma = 1;
  m_NumberOfIterations = 100;
  m_RelaxationFactor = 0.2;
  m_NumberOfPoints = 0;
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
  typedef itk::MinimumMaximumImageCalculator <InputImageType> MinimumMaximumImageCalculatorType;
  MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
  labelValues->SetImage(this->GetInput());
  labelValues->Compute();
  m_LevelValue = 0.5*(labelValues->GetMinimum() + labelValues->GetMaximum());

  // smoothing
  typedef itk::RecursiveGaussianImageFilter<InputImageType, FloatImageType> RecursiveGaussianImageFilterType;
  auto gaussian = RecursiveGaussianImageFilterType::New();
  gaussian->SetInput(this->GetInput());
  gaussian->SetSigma(m_Sigma);

  // fill holes after smoothing
  typedef itk::GrayscaleFillholeImageFilter<FloatImageType, FloatImageType> GrayscaleFillholeImageFilterType;
  GrayscaleFillholeImageFilterType::Pointer fillholes = GrayscaleFillholeImageFilterType::New();
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
  MarchingCubes mcubes = MarchingCubes::New();
  mcubes->SetInputData(convertor->GetOutput());
  mcubes->SetValue(0, m_LevelValue);
  mcubes->Update();
  m_Output = mcubes->GetOutput();

  // decimate surface
  if (m_NumberOfPoints > 0) {
    m_Reduction = 1 - (double) m_NumberOfPoints / (double) m_Output->GetNumberOfPoints();
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInputData(m_Output);
    decimate->SetTargetReduction(m_Reduction);
    decimate->SetPreserveTopology(true);
    decimate->SetSplitting(false);
    decimate->Update();
    m_Output = decimate->GetOutput();
  }

  typedef vtkSmartPointer<vtkSmoothPolyDataFilter> SmoothPolyData;
  SmoothPolyData smoother = SmoothPolyData::New();
  smoother->SetInputData(m_Output);
  smoother->SetNumberOfIterations(m_NumberOfIterations);
  smoother->SetRelaxationFactor(m_RelaxationFactor);
  smoother->Update();

  typedef vtkSmartPointer<vtkPolyDataNormals> PolyDataNormals;
  PolyDataNormals normals = PolyDataNormals::New();
  normals->SetInputData(smoother->GetOutput());
  normals->AutoOrientNormalsOn();
  normals->FlipNormalsOff();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->SplittingOff();
  normals->Update();

  m_Output = normals->GetOutput();

  this->CreateMesh();
}

template< typename TInputImage, typename TOutputMesh >
void BinaryMask3DMeshSource< TInputImage, TOutputMesh >::CreateMesh()
{
}


/** PrintSelf */
template< typename TInputImage, typename TOutputMesh >
void BinaryMask3DMeshSource< TInputImage, TOutputMesh >::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent
     << "ObjectValue: " << static_cast<itk::NumericTraits<unsigned char>::PrintType >( m_ObjectValue )
     << std::endl;
}
}
