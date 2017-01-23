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
        typedef itk::TranslationTransform<TParametersValueType, PointDimension> TranslationTransformType;
        TranslationTransformType::Pointer translation = TranslationTransformType::New();
        translation->Translate(m_Translation);
        m_Transform = translation;

        //define scales
        m_Scales.set_size(m_Transform->GetNumberOfParameters());

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i) {
          m_Scales[i] = m_TranslationScale;
        }
        break;
      }
      case EnumTransform::Euler3D:{
        // Euler3DTransform
        typedef itk::Euler3DTransform<TParametersValueType> Euler3DTransformType;
        Euler3DTransformType::Pointer Euler3D = Euler3DTransformType::New();
        Euler3D->SetIdentity();
        Euler3D->SetCenter(m_Center);
        Euler3D->SetTranslation(m_Translation);
        m_Transform = Euler3D;

        //define scales
        m_Scales.set_size(m_Transform->GetNumberOfParameters());
        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
        }
        break;
      }

      case EnumTransform::Similarity:{
        // Similarity3DTransform
        typedef itk::Similarity3DTransform<TParametersValueType> Similarity3DTransformType;
        Similarity3DTransformType::Pointer similarity3D = Similarity3DTransformType::New();
        similarity3D->SetIdentity();
        similarity3D->SetCenter(m_Center);
        similarity3D->SetTranslation(m_Translation);
        m_Transform = similarity3D;

        //define scales
        m_Scales.set_size(m_Transform->GetNumberOfParameters());
        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
        }

        m_Scales[count] = m_ScalingScale;
        break;
      }

      case EnumTransform::ScaleSkewVersor3D:{
        typedef itk::ScaleSkewVersor3DTransform<TParametersValueType> ScaleSkewVersor3DTransformType;
        ScaleSkewVersor3DTransformType::Pointer scaleskewversor3D = ScaleSkewVersor3DTransformType::New();
        scaleskewversor3D->SetIdentity();
        scaleskewversor3D->SetCenter(m_Center);
        scaleskewversor3D->SetTranslation(m_Translation);
        m_Transform = scaleskewversor3D;

        //define scales
        const size_t numberOfRotationComponents = 3;
        const size_t numberOfTranslationComponents = 3;
        const size_t numberOfScalingComponents = 3;
        const size_t numberOfSkewComponents = 6;

        m_Scales.set_size(m_Transform->GetNumberOfParameters());
        size_t count = 0;

        for (size_t i = 0; i < numberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
        }

        for (size_t i = 0; i < numberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = 0;
        }

        for (size_t i = 0; i < numberOfScalingComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingScale;
        }

        for (size_t i = 0; i < numberOfSkewComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingSkew;
        }

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
      os << "translation    " << m_Translation << std::endl;
      os << "center         " << m_Center << std::endl;
      os << "scales         " << m_Scales << std::endl;
      os << "parameters     " << m_Transform->GetParameters() << ", " << m_Transform->GetNumberOfParameters() << std::endl;
      os << std::endl;
    }

  protected:
    EnumTransform m_TypeOfTransform;
    typename TransformType::Pointer m_Transform;
    InputPointType m_Center;
    OutputVectorType m_Translation;
    itk::Array<double> m_Scales;
    double m_TranslationScale = 1;
    double m_RotationScale = 0.1;
    double m_ScalingScale = 0.1;
    double m_ScalingSkew = 0.1;

    size_t m_NumberOfTranslationComponents = PointDimension;
    size_t m_NumberOfRotationComponents = PointDimension;
    size_t m_NumberOfScalingComponents = PointDimension;

    TransformInitializer()
    {
      this->SetNumberOfRequiredInputs(0);
      this->SetNumberOfRequiredOutputs(0);
      m_TypeOfTransform = EnumTransform::Euler3D;
      m_Transform = nullptr;
    }
    ~TransformInitializer() {}
  };
}
