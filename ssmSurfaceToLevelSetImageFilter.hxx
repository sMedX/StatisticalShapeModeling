#pragma once

#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include "ssmSurfaceToLevelSetImageFilter.h"

namespace ssm
{
  template <typename TInputMesh, typename TOutputImage>
  SurfaceToLevelSetImageFilter<TInputMesh, TOutputImage>::SurfaceToLevelSetImageFilter()
  {
    m_UseOrigin = false;
    m_UseSize = false;
    m_Spacing.Fill(1);
  }

  template <typename TInputMesh, typename TOutputImage>
  void SurfaceToLevelSetImageFilter<TInputMesh, TOutputImage>::SetOrigin(typename TOutputImage::PointType origin)
  {
    m_Origin = origin;
    m_UseOrigin = true;
  }

  template <typename TInputMesh, typename TOutputImage>
  void SurfaceToLevelSetImageFilter<TInputMesh, TOutputImage>::SetSize(typename TOutputImage::SizeType size)
  {
    m_Size = size;
    m_UseSize = true;
  }

  template <typename TInputMesh, typename TOutputImage>
  void SurfaceToLevelSetImageFilter<TInputMesh, TOutputImage>::GenerateData()
  {
    m_BoundingBox = m_Input->GetBoundingBox();
    MaskImageType::SpacingType diff = m_BoundingBox->GetMaximum() - m_BoundingBox->GetMinimum();

    for (unsigned i = 0; i < Dimension; ++i) {
      //compute origin
      if (!m_UseOrigin) {
        double margin = m_Margin * diff[i];
        m_Origin[i] = m_BoundingBox->GetMinimum()[i] - margin;
      }

      //compute size
      if (!m_UseSize) {
        double margin = m_BoundingBox->GetMinimum()[i] - m_Origin[i];
        m_Size[i] = (diff[i] + 2 * margin) / m_Spacing[i];
      }
    }

    typedef itk::TriangleMeshToBinaryImageFilter<TInputMesh, MaskImageType> TriangleMeshToBinaryImageFilterType;
    TriangleMeshToBinaryImageFilterType::Pointer surfaceToImage = TriangleMeshToBinaryImageFilterType::New();
    surfaceToImage->SetInput(m_Input);
    surfaceToImage->SetInsideValue(1);
    surfaceToImage->SetSize(m_Size);
    surfaceToImage->SetSpacing(m_Spacing);
    surfaceToImage->SetOrigin(m_Origin);
    try {
      surfaceToImage->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      itkExceptionMacro(<< excep);
    }
    m_Mask = surfaceToImage->GetOutput();

    // compute minimum and maximum values
    typedef itk::MinimumMaximumImageCalculator <MaskImageType> MinimumMaximumImageCalculatorType;
    MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
    labelValues->SetImage(surfaceToImage->GetOutput());
    labelValues->Compute();
    if (labelValues->GetMaximum() < surfaceToImage->GetInsideValue()) {
      itkExceptionMacro(<< "warning: there is no region of interest in image");
    }

    m_BackgroundValue = labelValues->GetMinimum();
    m_ForegroundValue = labelValues->GetMaximum();

    // compute level set image
    typedef itk::SignedMaurerDistanceMapImageFilter<MaskImageType, TOutputImage> DistanceFilterType;
    DistanceFilterType::Pointer distanceToForeground = DistanceFilterType::New();
    distanceToForeground->SetInput(m_Mask);
    distanceToForeground->SetUseImageSpacing(true);
    distanceToForeground->SetBackgroundValue(m_BackgroundValue);
    distanceToForeground->SetInsideIsPositive(false);

    DistanceFilterType::Pointer distanceToBackground = DistanceFilterType::New();
    distanceToBackground->SetInput(m_Mask);
    distanceToBackground->SetUseImageSpacing(true);
    distanceToBackground->SetBackgroundValue(m_ForegroundValue);
    distanceToBackground->SetInsideIsPositive(true);

    typedef itk::AddImageFilter <TOutputImage> AddImageFilterType;
    AddImageFilterType::Pointer addfilter = AddImageFilterType::New();
    addfilter->SetInput1(distanceToForeground->GetOutput());
    addfilter->SetInput2(distanceToBackground->GetOutput());

    typedef itk::MultiplyImageFilter <TOutputImage> FilterType;
    FilterType::Pointer multiply = FilterType::New();
    multiply->SetInput(addfilter->GetOutput());
    multiply->SetConstant(0.5);
    try {
      multiply->Update();
    }
    catch (itk::ExceptionObject& excep) {
      std::cerr << excep << std::endl;
      itkExceptionMacro(<< excep);
    }

    m_Output = multiply->GetOutput();
  }
}
