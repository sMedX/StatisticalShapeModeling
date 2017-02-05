#pragma once

#include <itkTranslationTransform.h>
#include <itkEuler3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>

namespace ssm
{
  enum class EnumTransform
  {
    Translation,
    Euler3D,
    Similarity,
    ScaleSkewVersor3D
  };

  template<typename TParametersValueType = double>
  class TransformInitializer : public itk::ProcessObject
  {
  public:
    /** Standard class typedefs. */
    typedef TransformInitializer                     Self;
    typedef itk::ProcessObject                      Superclass;
    typedef itk::SmartPointer<Self>                 Pointer;
    typedef itk::SmartPointer<const Self>           ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(TransformInitializer, itk::ProcessObject);

    /** typedefs */
    itkStaticConstMacro(PointDimension, size_t, 3);
    typedef typename itk::Transform<TParametersValueType, PointDimension> TransformType;
    typedef typename TransformType::InputPointType InputPointType;
    typedef typename TransformType::OutputVectorType OutputVectorType;

    // Get transform
    itkGetObjectMacro(Transform, TransformType);

    // Set/Get type of transform
    itkSetEnumMacro(TypeOfTransform, EnumTransform)
    itkGetEnumMacro(TypeOfTransform, EnumTransform);
    void SetTypeOfTransform(const size_t & type) { m_TypeOfTransform = static_cast<EnumTransform>(type); }

    // Get scales
    itkGetMacro(Scales, itk::Array<double>);

    // Set/Get center and translation
    itkGetMacro(Center, InputPointType);
    itkSetMacro(Center, InputPointType);
    itkGetMacro(Translation, OutputVectorType);
    itkSetMacro(Translation, OutputVectorType);

    void Update()
    {
      switch (m_TypeOfTransform) {
      case EnumTransform::Translation: {
        // Translation transform
        typedef typename itk::TranslationTransform<TParametersValueType, PointDimension> TranslationTransformType;
        typename TranslationTransformType::Pointer translation = TranslationTransformType::New();
        translation->Translate(m_Translation);
        m_Transform = translation;
        m_Scales.set_size(m_Transform->GetNumberOfParameters());

        // define scales
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;
        this->SetScales(count, m_NumberOfTranslationComponents, m_TranslationScale);

        break;
      }
      case EnumTransform::Euler3D:{
        // Euler3DTransform
        typedef itk::Euler3DTransform<TParametersValueType> Euler3DTransformType;
        typename Euler3DTransformType::Pointer Euler3D = Euler3DTransformType::New();
        Euler3D->SetIdentity();
        Euler3D->SetCenter(m_Center);
        Euler3D->SetTranslation(m_Translation);
        m_Transform = Euler3D;
        m_Scales.set_size(m_Transform->GetNumberOfParameters());

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;
        this->SetScales(count, m_NumberOfRotationComponents, m_RotationScale);
        this->SetScales(count, m_NumberOfTranslationComponents, m_TranslationScale);

        break;
      }

      case EnumTransform::Similarity:{
        // Similarity3DTransform
        typedef itk::Similarity3DTransform<TParametersValueType> Similarity3DTransformType;
        typename Similarity3DTransformType::Pointer similarity3D = Similarity3DTransformType::New();
        similarity3D->SetIdentity();
        similarity3D->SetCenter(m_Center);
        similarity3D->SetTranslation(m_Translation);
        m_Transform = similarity3D;
        m_Scales.set_size(m_Transform->GetNumberOfParameters());

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 1;

        size_t count = 0;
        this->SetScales(count, m_NumberOfRotationComponents, m_RotationScale);
        this->SetScales(count, m_NumberOfTranslationComponents, m_TranslationScale);
        this->SetScales(count, m_NumberOfScalingComponents, m_ScalingScale);

        break;
      }

      case EnumTransform::ScaleSkewVersor3D:{
        typedef itk::ScaleSkewVersor3DTransform<TParametersValueType> ScaleSkewVersor3DTransformType;
        typename ScaleSkewVersor3DTransformType::Pointer scaleskewversor3D = ScaleSkewVersor3DTransformType::New();
        scaleskewversor3D->SetIdentity();
        scaleskewversor3D->SetCenter(m_Center);
        scaleskewversor3D->SetTranslation(m_Translation);
        m_Transform = scaleskewversor3D;
        m_Scales.set_size(m_Transform->GetNumberOfParameters());

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 3;
        m_NumberOfSkewComponents = 6;

        size_t count = 0;
        this->SetScales(count, m_NumberOfRotationComponents, m_RotationScale);
        this->SetScales(count, m_NumberOfTranslationComponents, m_TranslationScale);
        this->SetScales(count, m_NumberOfScalingComponents, m_ScalingScale);
        this->SetScales(count, m_NumberOfSkewComponents, m_SkewScale);

        break;
      }

      default:
        itkExceptionMacro(<< "Invalid type of transform");
      }
    }

    void PrintReport(std::ostream& os) const
    {
      os << this->GetNameOfClass() << std::endl;
      os << m_Transform->GetTransformTypeAsString() << ", category " << m_Transform->GetTransformCategory() << std::endl;
      os << "center               " << m_Center << std::endl;
      os << "translation          " << m_Translation << std::endl;
      os << "scales               " << m_Scales << std::endl;
      os << "fixed parameters     " << m_Transform->GetFixedParameters() << ", " << m_Transform->GetNumberOfFixedParameters() << std::endl;
      os << "transform parameters " << m_Transform->GetParameters() << ", " << m_Transform->GetNumberOfParameters() << std::endl;
      os << std::endl;
    }

  protected:
    EnumTransform m_TypeOfTransform = EnumTransform::Euler3D;
    typename TransformType::Pointer m_Transform;
    InputPointType m_Center;
    OutputVectorType m_Translation;
    itk::Array<double> m_Scales;
    double m_TranslationScale = 1;
    double m_RotationScale = 0.1;
    double m_ScalingScale = 0.1;
    double m_SkewScale = 0.1;

    size_t m_NumberOfTranslationComponents = 0;
    size_t m_NumberOfRotationComponents = 0;
    size_t m_NumberOfScalingComponents = 0;
    size_t m_NumberOfSkewComponents = 0;

    void SetScales(size_t & count, size_t numberOfComponents, double scale)
    {
      for (size_t i = 0; i < numberOfComponents; ++i, ++count) {
        m_Scales[count] = scale;
      }
    }

    TransformInitializer()
    {
      this->SetNumberOfRequiredInputs(0);
      this->SetNumberOfRequiredOutputs(0);
      m_Transform = nullptr;
    }
    ~TransformInitializer() {}
  };
}
