#pragma once

#include  <itkTriangleMeshToBinaryImageFilter.h>
#include  <itkSignedMaurerDistanceMapImageFilter.h>

#include "SurfaceToLevelSetImageFilter.h"

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
  BinaryImageType::SpacingType diff = m_BoundingBox->GetMaximum() - m_BoundingBox->GetMinimum();

  for (unsigned i = 0; i < Dimension; ++i) {
    //compute origin
    if ( !m_UseOrigin ) {
      double margin = m_Margin * diff[i];
      m_Origin[i] = m_BoundingBox->GetMinimum()[i] - margin;
    }

    //compute size
    if ( !m_UseSize ) {
      double margin = m_BoundingBox->GetMinimum()[i] - m_Origin[i];
      m_Size[i] = (diff[i] + 2 * margin) / m_Spacing[i];
    }
  }

  typedef itk::TriangleMeshToBinaryImageFilter<TInputMesh, BinaryImageType> TriangleMeshToBinaryImageFilterType;
  TriangleMeshToBinaryImageFilterType::Pointer surfaceToImage = TriangleMeshToBinaryImageFilterType::New();
  surfaceToImage->SetInsideValue(1);
  surfaceToImage->SetSize(m_Size);
  surfaceToImage->SetSpacing(m_Spacing);
  surfaceToImage->SetOrigin(m_Origin);
  surfaceToImage->SetInput(m_Input);

  try {
    surfaceToImage->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, TOutputImage> DistanceFilterType;
  DistanceFilterType::Pointer distanceMapFilter = DistanceFilterType::New();
  distanceMapFilter->SetUseImageSpacing(true);
  distanceMapFilter->SetInput(surfaceToImage->GetOutput());

  try {
    distanceMapFilter->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cerr << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  m_Output = distanceMapFilter->GetOutput();
}
