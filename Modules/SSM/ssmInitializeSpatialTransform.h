#pragma once

#include <itkTranslationTransform.h>
#include <itkEuler3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkScaleSkewVersor3DTransform.h>
#include <itkLogger.h>

namespace ssm
{
  template<typename TParametersValueType = double>
  class InitializeSpatialTransform : public itk::Object
  {
  public:
    /** Standard class typedefs. */
    typedef InitializeSpatialTransform     Self;
    typedef itk::Object                    Superclass;
    typedef itk::SmartPointer<Self>        Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);
    itkTypeMacro(InitializeSpatialTransform, itk::Object);

    enum class Transform
    {
      Translation,
      Euler3D,
      Similarity,
      ScaleSkewVersor3D
    };

    /** typedefs */
    itkStaticConstMacro(Dimension, unsigned int, 3U);

    typedef typename itk::Transform<TParametersValueType, Dimension> TransformType;
    typedef typename TransformType::InputPointType InputPointType;
    typedef typename TransformType::OutputVectorType OutputVectorType;
    typedef itk::Array<double> ParametersType;
    typedef itk::Array<unsigned int> ModeBoundsType;

    // Set logger
    itkSetObjectMacro(Logger, itk::Logger);

    // Get transform
    itkGetObjectMacro(Transform, TransformType);

    // Set/Get type of transform
    itkSetEnumMacro(TransformType, Transform);
    itkGetEnumMacro(TransformType, Transform);
    void SetTransformType(const size_t & type) { this->SetTransformType(static_cast<Transform>(type)); }

    itkSetMacro(RotationScale, double);
    itkGetMacro(RotationScale, double);

    itkSetMacro(TranslationScale, double);
    itkGetMacro(TranslationScale, double);

    itkSetMacro(ScalingScale, double);
    itkGetMacro(ScalingScale, double);

    itkSetMacro(SkewScale, double);
    itkGetMacro(SkewScale, double);

    // Get scales and bounds
    itkGetMacro(Scales, ParametersType);
    itkGetMacro(ModeBounds, ModeBoundsType);
    itkGetMacro(LowerBounds, ParametersType);
    itkGetMacro(UpperBounds, ParametersType);

    // Set/Get center and translation
    itkGetMacro(Center, InputPointType);
    itkSetMacro(Center, InputPointType);
    itkGetMacro(Translation, OutputVectorType);
    itkSetMacro(Translation, OutputVectorType);

    void Initialize()
    {
      switch (m_TransformType) {
      case Transform::Translation: {
        // Translation transform
        typedef typename itk::TranslationTransform<TParametersValueType, Dimension> TranslationTransformType;
        auto transform = TranslationTransformType::New();
        transform->Translate(m_Translation);
        m_Transform = transform;

        this->Allocate();

        // define scales
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }
      case Transform::Euler3D:{
        // Euler3DTransform
        typedef itk::Euler3DTransform<TParametersValueType> Euler3DTransformType;
        auto transform = Euler3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);
        m_Transform = transform;

        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        break;
      }

      case Transform::Similarity:{
        // Similarity3DTransform
        typedef itk::Similarity3DTransform<TParametersValueType> Similarity3DTransformType;
        auto transform = Similarity3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);
        m_Transform = transform;

        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 1;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        for (size_t i = 0; i < m_NumberOfScalingComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingScale;
          m_ModeBounds[count] = 2;
        }

        break;
      }

      case Transform::ScaleSkewVersor3D:{
        typedef itk::ScaleSkewVersor3DTransform<TParametersValueType> ScaleSkewVersor3DTransformType;
        auto transform = ScaleSkewVersor3DTransformType::New();
        transform->SetIdentity();
        transform->SetCenter(m_Center);
        transform->SetTranslation(m_Translation);
        m_Transform = transform;

        this->Allocate();

        // define scales
        m_NumberOfRotationComponents = 3;
        m_NumberOfTranslationComponents = 3;
        m_NumberOfScalingComponents = 3;
        m_NumberOfSkewComponents = 6;

        size_t count = 0;

        for (size_t i = 0; i < m_NumberOfRotationComponents; ++i, ++count) {
          m_Scales[count] = m_RotationScale;
          m_ModeBounds[count] = 2;
        }

        for (size_t i = 0; i < m_NumberOfTranslationComponents; ++i, ++count) {
          m_Scales[count] = m_TranslationScale;
          m_ModeBounds[count] = 0;
        }

        for (size_t i = 0; i < m_NumberOfScalingComponents; ++i, ++count) {
          m_Scales[count] = m_ScalingScale;
          m_ModeBounds[count] = 2;
        }

        for (size_t i = 0; i < m_NumberOfSkewComponents; ++i, ++count) {
          m_Scales[count] = m_SkewScale;
          m_ModeBounds[count] = 2;
        }

        break;
      }
      }

      for (size_t n = 0; n < m_NumberOfParameters; ++n) {
        if (m_ModeBounds[n] == 1 || m_ModeBounds[n] == 2) {
          m_LowerBounds[n] = m_Transform->GetParameters()[n] - m_Scales[n];
        }
        if (m_ModeBounds[n] == 2 || m_ModeBounds[n] == 3) {
          m_UpperBounds[n] = m_Transform->GetParameters()[n] + m_Scales[n];
        }
      }
    }

    void PrintReport() const
    {
      std::cout << this->GetNameOfClass() << std::endl;
      std::cout << "transform        " << m_Transform->GetTransformTypeAsString() << std::endl;
      std::cout << "center           " << m_Center << std::endl;
      std::cout << "translation      " << m_Translation << std::endl;
      std::cout << "fixed parameters " << m_Transform->GetFixedParameters() << ", " << m_Transform->GetNumberOfFixedParameters() << std::endl;
      std::cout << "parameters       " << m_Transform->GetParameters() << ", " << m_Transform->GetNumberOfParameters() << std::endl;
      std::cout << std::endl;
      std::cout << "scales       " << m_Scales << std::endl;
      std::cout << "mode bounds  " << m_ModeBounds << std::endl;
      std::cout << "lower bounds " << m_LowerBounds << std::endl;
      std::cout << "upper bounds " << m_UpperBounds << std::endl;
      std::cout << std::endl;
    }

  private:
    typename TransformType::Pointer m_Transform = nullptr;
    Transform m_TransformType = Transform::Similarity;
    InputPointType m_Center;
    OutputVectorType m_Translation;

    /** Set the boundary condition for each variable, where
    * select[i] = 0 if x[i] is unbounded,
    *           = 1 if x[i] has only a lower bound,
    *           = 2 if x[i] has both lower and upper bounds, and
    *           = 3 if x[1] has only an upper bound */
    ModeBoundsType m_ModeBounds;
    ParametersType m_LowerBounds;
    ParametersType m_UpperBounds;
    ParametersType m_Scales;

    size_t m_NumberOfParameters;
    size_t m_NumberOfTranslationComponents;
    size_t m_NumberOfRotationComponents;
    size_t m_NumberOfScalingComponents;
    size_t m_NumberOfSkewComponents;

    double m_TranslationScale = 1;
    double m_RotationScale = 0.2;
    double m_ScalingScale = 0.2;
    double m_SkewScale = 0.2;

  protected:
    InitializeSpatialTransform() {}
    ~InitializeSpatialTransform() {}

    void Allocate()
    {
      m_NumberOfParameters = m_Transform->GetNumberOfParameters();

      m_Scales.set_size(m_NumberOfParameters);
      m_Scales.Fill(NAN);

      m_ModeBounds.set_size(m_NumberOfParameters);
      m_ModeBounds.Fill(NAN);

      m_LowerBounds.set_size(m_NumberOfParameters);
      m_LowerBounds.Fill(NAN);

      m_UpperBounds.set_size(m_NumberOfParameters);
      m_UpperBounds.Fill(NAN);
    }

    itk::Logger::Pointer m_Logger = itk::Logger::New();
    mutable std::ostringstream m_Message;
  };
}
