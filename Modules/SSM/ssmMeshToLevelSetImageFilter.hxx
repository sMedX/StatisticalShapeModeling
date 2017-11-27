#pragma once

#include <itkTriangleMeshToBinaryImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkNumericTraits.h>

#include "ssmMeshToLevelSetImageFilter.h"

namespace ssm
{
template <typename TInputMesh, typename TOutputImage>
MeshToLevelSetImageFilter<TInputMesh, TOutputImage>::MeshToLevelSetImageFilter()
{
  this->SetNumberOfRequiredInputs(1);

  m_UseOrigin = false;
  m_UseSize = false;

  m_Margin = 0.10;
  m_Spacing.Fill(1);

  m_ForegroundValue = itk::NumericTraits< BinaryPixelype >::OneValue();
  m_BackgroundValue = itk::NumericTraits< BinaryPixelype >::ZeroValue();
}

/** Set the Input Mesh */
template< typename TInputMesh, typename TOutputImage >
void MeshToLevelSetImageFilter< TInputMesh, TOutputImage >::SetInput(TInputMesh *input)
{
  this->ProcessObject::SetNthInput(0, input);
}

/** Get the input Mesh */
template< typename TInputMesh, typename TOutputImage >
typename MeshToLevelSetImageFilter< TInputMesh, TOutputImage >::InputMeshType *
MeshToLevelSetImageFilter< TInputMesh, TOutputImage >::GetInput(void)
{
  return static_cast<TInputMesh *> (this->ProcessObject::GetInput(0));
}

/** Get the input Mesh */
template< typename TInputMesh, typename TOutputImage >
typename MeshToLevelSetImageFilter< TInputMesh, TOutputImage >::InputMeshType *
MeshToLevelSetImageFilter< TInputMesh, TOutputImage >::GetInput(unsigned int idx)
{
  return static_cast< TInputMesh * >(this->ProcessObject::GetInput(idx));
}

template <typename TInputMesh, typename TOutputImage>
void MeshToLevelSetImageFilter<TInputMesh, TOutputImage>::SetOrigin(const ImagePointType & point)
{
  m_Origin = point;
  m_UseOrigin = true;
}

template <typename TInputMesh, typename TOutputImage>
void MeshToLevelSetImageFilter<TInputMesh, TOutputImage>::SetSize(const SizeType & size)
{
  m_Size = size;
  m_UseSize = true;
}

template <typename TInputMesh, typename TOutputImage>
void MeshToLevelSetImageFilter<TInputMesh, TOutputImage>::GenerateData()
{
  auto bbox = this->GetInput()->GetBoundingBox();
  auto diff = bbox->GetMaximum() - bbox->GetMinimum();

  // compute origin and size
  for (size_t i = 0; i < Dimension; ++i) {
    if (!m_UseOrigin) {
      m_Origin[i] = bbox->GetMinimum()[i] - m_Margin * diff[i];
    }

    if (!m_UseSize) {
      m_Size[i] = diff[i] * (1 + 2 * m_Margin) / m_Spacing[i];
    }
  }

  typedef itk::TriangleMeshToBinaryImageFilter<TInputMesh, BinaryImageType> TriangleMeshToBinaryImageFilterType;
  auto surfaceToImage = TriangleMeshToBinaryImageFilterType::New();
  surfaceToImage->SetInput(this->GetInput());
  surfaceToImage->SetInsideValue(m_ForegroundValue);
  surfaceToImage->SetOutsideValue(m_BackgroundValue);
  surfaceToImage->SetSize(m_Size);
  surfaceToImage->SetSpacing(m_Spacing);
  surfaceToImage->SetOrigin(m_Origin);
  surfaceToImage->Update();
  m_Mask = surfaceToImage->GetOutput();

  // compute level set image
  typedef itk::SignedMaurerDistanceMapImageFilter<BinaryImageType, TOutputImage> DistanceFilterType;
  auto distanceToForeground = DistanceFilterType::New();
  distanceToForeground->SetInput(surfaceToImage->GetOutput());
  distanceToForeground->SetUseImageSpacing(true);
  distanceToForeground->SetBackgroundValue(m_BackgroundValue);
  distanceToForeground->SetInsideIsPositive(false);

  typename DistanceFilterType::Pointer distanceToBackground = DistanceFilterType::New();
  distanceToBackground->SetInput(surfaceToImage->GetOutput());
  distanceToBackground->SetUseImageSpacing(true);
  distanceToBackground->SetBackgroundValue(m_ForegroundValue);
  distanceToBackground->SetInsideIsPositive(true);

  typedef itk::AddImageFilter <TOutputImage> AddImageFilterType;
  auto add = AddImageFilterType::New();
  add->SetInput1(distanceToForeground->GetOutput());
  add->SetInput2(distanceToBackground->GetOutput());

  typedef itk::MultiplyImageFilter <TOutputImage> FilterType;
  auto multiply = FilterType::New();
  multiply->SetInput(add->GetOutput());
  multiply->SetConstant(0.5);
  multiply->Update();

  auto output = this->GetOutput();
  output->Graft(multiply->GetOutput());
}
}
