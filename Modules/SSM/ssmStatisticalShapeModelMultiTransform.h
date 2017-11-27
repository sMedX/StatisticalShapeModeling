#pragma once

#include <itkTransform.h>
#include <itkStatisticalModel.h>

namespace ssm
{
template<typename TDataset, typename TParametersValueType = double>
class StatisticalShapeModelMultiTransform : public itk::Transform<TParametersValueType, TDataset::PointDimension>
{
public:
  /** Standard class typedefs. */
  typedef StatisticalShapeModelMultiTransform                            Self;
  typedef itk::Transform<TParametersValueType, TDataset::PointDimension> Superclass;
  typedef itk::SmartPointer<Self>                                        Pointer;
  typedef itk::SmartPointer<const Self>                                  ConstPointer;

  /** New macro for creation of through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StatisticalShapeModelMultiTransform, Transform);

  /** Dimension of the domain space. */
  itkStaticConstMacro(Dimension, unsigned int, TDataset::PointDimension);

  /** Standard scalar type for this class. */
  typedef TDataset                                                    DatasetType;
  typedef itk::StatisticalModel<DatasetType>                          ShapeModelType;
  typedef itk::Transform<TParametersValueType, Dimension, Dimension>  SpatialTransformType;

  /** Standard parameters container. */
  typedef TParametersValueType                      ParametersValueType;
  typedef typename Superclass::ParametersType       ParametersType;
  typedef typename Superclass::FixedParametersType  FixedParametersType;
  typedef typename ShapeModelType::VectorType       ShapeModelParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;

  /** The number of parameters defining this transform. */
  typedef typename Superclass::NumberOfParametersType NumberOfParametersType;

  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType  InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  /** This method sets the parameters for the transform value specified by the user. */
  virtual void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;

  /** Get the Transformation Parameters. */
  virtual const ParametersType & GetParameters() const ITK_OVERRIDE;

  /** Transform point. */
  OutputPointType TransformPoint(const InputPointType  & point) const ITK_OVERRIDE;
  OutputPointType TransformPoint(const size_t  & index) const;

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual void ComputeJacobianWithRespectToParameters(const InputPointType & point, JacobianType & j) const ITK_OVERRIDE;

  /**
  * Expanded interface to Compute the Jacobian with respect to the parameters for the composite transform using Jacobian rule. This version takes in temporary
  * variables to avoid excessive constructions. NOTE: outJacobian and jacobianWithRespectToPosition MUST be sized
  * prior to the call; outJacobian's size should be [NDimensions, this->GetNumberOfLocalParameters() ]
  * jacobianWithRespectToPosition size == [ NDimensions, NDimensions ]
  */
  virtual void ComputeJacobianWithRespectToParametersCachedTemporaries(const InputPointType & p, JacobianType & outJacobian, JacobianType & jacobianWithRespectToPosition) const ITK_OVERRIDE;
  virtual void ComputeJacobianWithRespectToParametersCachedTemporaries(const size_t & index, JacobianType & outJacobian, JacobianType & modelJacobian, JacobianType & spatialJacobian, JacobianType & jacobianWithRespectToPosition) const;

  /** Return the number of parameters that completely define the Transform  */
  virtual NumberOfParametersType GetNumberOfParameters() const ITK_OVERRIDE
  {
    return m_NumberOfParameters;
  }

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters(const FixedParametersType &) ITK_OVERRIDE
  {
    itkExceptionMacro(<< "method SetFixedParametersType is not implemented");
  }

  /** Get the Fixed Parameters. */
  virtual const FixedParametersType & GetFixedParameters() const ITK_OVERRIDE
  {
    itkExceptionMacro(<<"method GetFixedParametersType is not implemented");
  }

  itkGetConstReferenceMacro(NumberOfUsedComponents, NumberOfParametersType);
  itkGetConstReferenceMacro(NumberOfComponents, NumberOfParametersType);

  virtual void SetNumberOfUsedComponents(const NumberOfParametersType arg)
  {
    if (!m_ShapeModel) {
      itkExceptionMacro(<< "Number of used components must be initialized after shape model.")
    }

    if (m_NumberOfUsedComponents != arg) {
      m_NumberOfUsedComponents = arg;
      this->Modified();
    }
  }

  /**Set/Get methods for spatial transform */
  itkGetConstObjectMacro(SpatialTransform, SpatialTransformType);
  void SetSpatialTransform(SpatialTransformType * arg)
  {
    itkDebugMacro("setting m_SpatialTransform to " << arg);
    if (m_SpatialTransform != arg) {
      m_SpatialTransform = arg;
      m_NumberOfParameters = this->m_SpatialTransform->GetNumberOfParameters() + m_NumberOfComponents;
      this->Modified();
    }
  }

  /**Set/Get methods for shape model*/
  itkGetConstObjectMacro(ShapeModel, ShapeModelType);
  void SetShapeModel(const ShapeModelType * arg)
  {
    itkDebugMacro("setting m_ShapeModel to " << arg);
    if (m_ShapeModel != arg) {
      m_ShapeModel = arg;
      m_NumberOfComponents = m_ShapeModel->GetNumberOfPrincipalComponents();
      m_NumberOfParameters = m_SpatialTransform->GetNumberOfParameters() + m_NumberOfComponents;

      m_NumberOfUsedComponents = m_NumberOfComponents;
      m_ShapeModelParameters.set_size(m_NumberOfUsedComponents);
      m_ShapeModelParameters.fill(0);

      this->Modified();
    }
  }

protected:
  StatisticalShapeModelMultiTransform();
  ~StatisticalShapeModelMultiTransform() {};
  /** Print contents of an TranslationTransform. */
  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

private:
  StatisticalShapeModelMultiTransform(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  typename ShapeModelType::ConstPointer m_ShapeModel;
  typename SpatialTransformType::Pointer m_SpatialTransform;

  ShapeModelParametersType m_ShapeModelParameters;
  NumberOfParametersType m_NumberOfComponents;
  NumberOfParametersType m_NumberOfParameters;
  NumberOfParametersType m_NumberOfUsedComponents;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmStatisticalShapeModelMultiTransform.hxx"
#endif
