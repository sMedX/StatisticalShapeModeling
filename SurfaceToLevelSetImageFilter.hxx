#ifndef __SurfaceToLevelSetImageFilter_hxx
#define __SurfaceToLevelSetImageFilter_hxx

#include "SurfaceToLevelSetImageFilter.h"

template <typename TInputMesh, typename TOutputImage>
void SurfaceToLevelSetImageFilter<TInputMesh, TOutputImage>::GenerateData()
{
  typename TInputMesh::BoundingBoxType::ConstPointer boundingBox = m_Input->GetBoundingBox();
  BinaryImageType::SpacingType diff = boundingBox->GetMaximum() - boundingBox->GetMinimum();
  BinaryImageType::PointType origin;
  BinaryImageType::SpacingType spacing;
  spacing.Fill(m_Spacing);

  BinaryImageType::SizeType size;

  for (unsigned i = 0; i < Dimension; ++i) {
    // margin on each side
    double margin = m_Margin * diff[i];

    //compute size and origin
    size[i] = (diff[i] + 2 * margin) / spacing[i];
    origin[i] = boundingBox->GetMaximum()[i] - margin;
  }

  typedef itk::TriangleMeshToBinaryImageFilter<TInputMesh, BinaryImageType> TriangleMeshToBinaryImageFilterType;
  TriangleMeshToBinaryImageFilterType::Pointer surfaceToImage = TriangleMeshToBinaryImageFilterType::New();
  surfaceToImage->SetInsideValue(1);
  surfaceToImage->SetSize(size);
  surfaceToImage->SetSpacing(spacing);
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

#endif
