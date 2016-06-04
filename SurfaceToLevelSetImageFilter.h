#ifndef __SurfaceToLevelSetImageFilter_h
#define __SurfaceToLevelSetImageFilter_h

#include <itkMesh.h>
#include <itkImage.h>
#include <itkSmartPointer.h>
#include <itkLightProcessObject.h>

template <typename TInputMesh, typename TOutputImage>
class SurfaceToLevelSetImageFilter : public itk::LightProcessObject
{
public:
  typedef SurfaceToLevelSetImageFilter Self;
  typedef itk::LightProcessObject Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(SurfaceToLevelSetImageFilter, itk::LightProcessObject);

public:

  //inputs and output
  itkSetObjectMacro(Input, TInputMesh);
  itkGetConstObjectMacro(Output, TOutputImage);

  // parameters
  itkSetMacro(Margin, double);
  itkGetMacro(Margin, double);

  itkSetMacro(Spacing, double);
  itkGetMacro(Spacing, double);

public:
  void Update() { this->UpdateOutputData(); }

protected:
  SurfaceToLevelSetImageFilter() {};
  virtual ~SurfaceToLevelSetImageFilter() {};

  // it makes all work
  void GenerateData();

protected:
  typename TInputMesh::Pointer m_Input;
  typename TOutputImage::ConstPointer m_Output;

  double m_Margin = 0.25;
  double m_Spacing = 1;

private:
  SurfaceToLevelSetImageFilter(const Self&);
  void operator=(const Self&);
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "SurfaceToLevelSetImageFilter.hxx"
#endif

#endif
