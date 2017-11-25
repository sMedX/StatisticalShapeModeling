#pragma once

#include <itkImageSource.h>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkLightProcessObject.h>

namespace ssm
{
  template <typename TInputMesh, typename TOutputImage>
  class SurfaceToLevelSetImageFilter : public itk::ImageSource<TOutputImage>
  {
  public:
    typedef SurfaceToLevelSetImageFilter    Self;
    typedef itk::ImageSource<TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>   ConstPointer;

    itkNewMacro(Self);
    itkTypeMacro(SurfaceToLevelSetImageFilter, itk::ImageSource);

    /** Constants for the image dimensions */
    itkStaticConstMacro(Dimension, unsigned int, 3U);
    static_assert(TInputMesh::PointDimension == Dimension, "Invalid dimension of the input mesh.");
    static_assert(TOutputImage::ImageDimension == Dimension, "Invalid dimension of the output image.");

    typedef TInputMesh InputMeshType;
    typedef TOutputImage OutputImageType;
    typedef unsigned char BinaryPixelype;
    typedef itk::Image<BinaryPixelype, Dimension> BinaryImageType;
    typedef typename OutputImageType::PointType ImagePointType;
    typedef typename OutputImageType::SpacingType  SpacingType;
    typedef typename TOutputImage::SizeType SizeType;

    /** Set the mesh of this process object.  */
    using Superclass::SetInput;
    void SetInput(InputMeshType *input);

    /** Get the mesh input of this process object.  */
    InputMeshType * GetInput();
    InputMeshType * GetInput(unsigned int idx);

    itkGetObjectMacro(Mask, BinaryImageType);

    // parameters
    itkSetMacro(Margin, double);
    itkGetMacro(Margin, double);

    void SetOrigin(const ImagePointType & point);
    itkGetMacro(Origin, ImagePointType);

    void SetSize(const SizeType & size);
    itkGetMacro(Size, SizeType);

    itkSetMacro(Spacing, SpacingType);
    itkGetMacro(Spacing, SpacingType);

    itkSetMacro(ForegroundValue, BinaryPixelype);
    itkGetMacro(ForegroundValue, BinaryPixelype);

    itkSetMacro(BackgroundValue, BinaryPixelype);
    itkGetMacro(BackgroundValue, BinaryPixelype);

  protected:
    SurfaceToLevelSetImageFilter();
    virtual ~SurfaceToLevelSetImageFilter() {};

    virtual void GenerateOutputInformation() ITK_OVERRIDE {} // do nothing
    virtual void GenerateData() ITK_OVERRIDE;

    typename BinaryImageType::Pointer m_Mask;
    ImagePointType m_Origin;
    bool m_UseOrigin;

    SpacingType m_Spacing;
    bool m_UseSpacing;

    SizeType m_Size;
    bool m_UseSize;

    double m_Margin;

    BinaryPixelype m_BackgroundValue;
    BinaryPixelype m_ForegroundValue;

  private:
    SurfaceToLevelSetImageFilter(const Self&);
    void operator=(const Self&);
  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmSurfaceToLevelSetImageFilter.hxx"
#endif
