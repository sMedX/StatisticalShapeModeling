#pragma once

#include <itkMesh.h>
#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkLightProcessObject.h>

namespace ssm
{
  template <typename TInputMesh, typename TOutputImage>
  class SurfaceToLevelSetImageFilter : public itk::LightProcessObject
  {
  public:
    typedef SurfaceToLevelSetImageFilter   Self;
    typedef itk::LightProcessObject        Superclass;
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(SurfaceToLevelSetImageFilter, itk::LightProcessObject);

    itkStaticConstMacro(Dimension, unsigned int, TInputMesh::PointDimension);
    typedef itk::Image<unsigned char, Dimension> MaskImageType;

  public:
    //inputs and output
    itkSetObjectMacro(Input, TInputMesh);
    itkGetObjectMacro(Output, TOutputImage);
    itkGetObjectMacro(Mask, MaskImageType);

    // parameters
    itkSetMacro(Margin, double);
    itkGetMacro(Margin, double);

    void SetOrigin(typename TOutputImage::PointType origin);
    itkGetMacro(Origin, typename TOutputImage::PointType);

    void SetSize(typename TOutputImage::SizeType size);
    itkGetMacro(Size, typename TOutputImage::SizeType);

    itkSetMacro(Spacing, typename TOutputImage::SpacingType);
    itkGetMacro(Spacing, typename TOutputImage::SpacingType);

  public:
    void Update() { this->UpdateOutputData(); }

  protected:
    SurfaceToLevelSetImageFilter();
    virtual ~SurfaceToLevelSetImageFilter() {};

    // it makes all work
    void GenerateData();

  protected:
    typename TInputMesh::Pointer m_Input;
    typename TOutputImage::Pointer m_Output;
    typename MaskImageType::Pointer m_Mask;
    typename TInputMesh::BoundingBoxType::ConstPointer m_BoundingBox;

    typename TOutputImage::PointType m_Origin;
    bool m_UseOrigin;

    typename TOutputImage::SpacingType m_Spacing;
    bool m_UseSpacing;

    typename TOutputImage::SizeType m_Size;
    bool m_UseSize;
    double m_Margin = 0.25;

    float m_BackgroundValue;
    float m_ForegroundValue;

  private:
    SurfaceToLevelSetImageFilter(const Self&);
    void operator=(const Self&);
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmSurfaceToLevelSetImageFilter.hxx"
#endif
