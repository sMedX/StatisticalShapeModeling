#pragma once

#include <itkImageToImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>

namespace ssm
{
  template< typename TInputImage, typename TOutputImage>
  class BinaryImageToLevelSetImageFilter : public itk::ImageToImageFilter <TInputImage, TOutputImage>
  {
  public:
    /** Standard class typedefs. */
    typedef BinaryImageToLevelSetImageFilter                           Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef itk::SmartPointer< Self >                                  Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BinaryImageToLevelSetImageFilter, itk::ImageToImageFilter);

  protected:
    BinaryImageToLevelSetImageFilter() {};
    virtual ~BinaryImageToLevelSetImageFilter() {}

    /** Does the real work. */
    void GenerateData()
    {
      // define input image
      typename TInputImage::ConstPointer input = this->GetInput();

      // compute minimum and maximum values
      typedef itk::MinimumMaximumImageCalculator <TInputImage> MinimumMaximumImageCalculatorType;
      MinimumMaximumImageCalculatorType::Pointer labelValues = MinimumMaximumImageCalculatorType::New();
      labelValues->SetImage(input);
      labelValues->Compute();

      m_BackgroundValue = labelValues->GetMinimum();
      m_ForegroundValue = labelValues->GetMaximum();

      if (m_BackgroundValue == m_ForegroundValue) {
        itkExceptionMacro(<< "warning: there is no region of interest in the input image");
      }

      // compute level set image
      typedef typename itk::SignedMaurerDistanceMapImageFilter<TInputImage, TOutputImage> DistanceFilterType;
      DistanceFilterType::Pointer distanceToForeground = DistanceFilterType::New();
      distanceToForeground->SetInput(input);
      distanceToForeground->SetUseImageSpacing(true);
      distanceToForeground->SetBackgroundValue(m_BackgroundValue);
      distanceToForeground->SetInsideIsPositive(false);

      DistanceFilterType::Pointer distanceToBackground = DistanceFilterType::New();
      distanceToBackground->SetInput(input);
      distanceToBackground->SetUseImageSpacing(true);
      distanceToBackground->SetBackgroundValue(m_ForegroundValue);
      distanceToBackground->SetInsideIsPositive(true);

      typedef typename itk::AddImageFilter <TOutputImage> AddImageFilterType;
      AddImageFilterType::Pointer add = AddImageFilterType::New();
      add->SetInput1(distanceToForeground->GetOutput());
      add->SetInput2(distanceToBackground->GetOutput());

      typedef typename itk::MultiplyImageFilter <TOutputImage> MultiplyImageFilterType;
      MultiplyImageFilterType::Pointer multiply = MultiplyImageFilterType::New();
      multiply->SetInput(add->GetOutput());
      multiply->SetConstant(0.5);
      multiply->Update();

      // set output
      this->GraftOutput(multiply->GetOutput());
    }

  private:
    BinaryImageToLevelSetImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    bool m_UseImageSpacing = true;
    typename TInputImage::PixelType m_BackgroundValue;
    typename TInputImage::PixelType m_ForegroundValue;
  };

}
