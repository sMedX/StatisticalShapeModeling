#pragma once

#include <itkImage.h>
#include <itkStatisticalModel.h>
#include <itkTransform.h>
#include <itkLBFGSOptimizer.h>

#include "ssmShapeModelRegistrationMethodBase.h"
#include "ssmStatisticalShapeModelMultiTransform.h"
#include "ssmShapeModelToLevelSetImageMetric.h"

namespace ssm
{
  template <typename TShapeModel, typename TInputImage, typename TOutputMesh>
  class ShapeModelToImageRegistrationMethod : public ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>
  {
  public:
    /** Standard class typedefs. */
    typedef ShapeModelToImageRegistrationMethod                         Self;
    typedef ShapeModelRegistrationMethodBase<TShapeModel, TOutputMesh>  Superclass;
    typedef itk::SmartPointer<Self>                                     Pointer;
    typedef itk::SmartPointer<const Self>                               ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(ShapeModelToImageRegistrationMethod, Superclass);

    /** Constants for the image dimensions */
    itkStaticConstMacro(ImageDimension, unsigned int, 3U);
    itkStaticConstMacro(PointDimension, unsigned int, itkGetStaticConstMacro(PointDimension));
    static_assert(TInputImage::ImageDimension == ImageDimension, "Invalid dimension of input image. Dimension 3 is supported.");

    /** Typedefs of the shape model. */
    typedef typename Superclass::ShapeModelType ShapeModelType;
    typedef typename Superclass::OutputMeshType OutputMeshType;
    typedef TInputImage InputImageType;

    /** Typedefs of the input image. */
    typedef itk::Image<unsigned char, ImageDimension> BinaryImageType;
    typedef itk::Image<float, ImageDimension> LevelSetImageType;

    /**  Type of the Transform Base class */
    typedef typename Superclass::ShapeModelMultiTransformType ShapeModelMultiTransformType;
    typedef ShapeModelToLevelSetImageMetric<ShapeModelType, ShapeModelMultiTransformType, LevelSetImageType> MetricType;
    typedef itk::LBFGSOptimizer OptimizerType;

    // Set/Get level set image
    itkSetConstObjectMacro(Image, InputImageType);
    itkGetConstObjectMacro(Image, InputImageType);
    itkGetConstObjectMacro(LevelSetImage, LevelSetImageType);

    itkSetObjectMacro(Optimizer, OptimizerType);
    itkGetObjectMacro(Optimizer, OptimizerType);

    itkSetObjectMacro(Metric, MetricType);
    itkGetObjectMacro(Metric, MetricType);

    itkSetMacro(ComputeLevelSetImage, bool);
    itkGetMacro(ComputeLevelSetImage, bool);

    void PrintReport();

  protected:
    ShapeModelToImageRegistrationMethod();
    ~ShapeModelToImageRegistrationMethod() {};

    void GenerateData();
    void Initialize();
    void ComputeLevelSet();
    void InitializeMetric();

    typename InputImageType::ConstPointer m_Image;
    typename LevelSetImageType::ConstPointer m_LevelSetImage;
    typename MetricType::Pointer m_Metric;
    bool m_ComputeLevelSetImage;
  };
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmShapeModelToImageRegistrationMethod.hxx"
#endif
