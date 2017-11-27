#pragma once

#include <itkIdentityTransform.h>

#include "ssmStatisticalShapeModelMultiTransform.h"

namespace ssm
{

template<typename TDataset, typename TParametersValueType>
StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>
::StatisticalShapeModelMultiTransform():Superclass(0)
{
  m_ShapeModel = nullptr;
  m_SpatialTransform = itk::IdentityTransform<TParametersValueType, Dimension>::New();

  m_NumberOfComponents = 0;
  m_NumberOfParameters = 0;
  m_NumberOfUsedComponents = 0;
}

template<typename TDataset, typename TParametersValueType>
void StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::SetParameters(const ParametersType & parameters)
{
  /* Verify proper input size. */
  if (parameters.Size() != this->GetNumberOfParameters()) {
    itkExceptionMacro(<< "Input parameter list size is not expected size. " << parameters.Size() << " instead of " << this->GetNumberOfParameters() << ".");
  }

  for (size_t i = 0; i < m_NumberOfUsedComponents; ++i) {
    m_ShapeModelParameters[i] = parameters[i];
  }

  for (size_t i = m_NumberOfUsedComponents; i < m_NumberOfComponents; ++i) {
    m_ShapeModelParameters[i] = 0;
  }

  /* set parameters to spatial transform */
  const auto parameterSize = m_SpatialTransform->GetParameters().Size();
  m_SpatialTransform->CopyInParameters(&(parameters.data_block())[m_NumberOfComponents],  &(parameters.data_block())[m_NumberOfComponents] + parameterSize);

  this->Modified();
}


template<typename TDataset, typename TParametersValueType>
const typename StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::ParametersType &
StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::GetParameters() const
{
  this->m_Parameters.SetSize(this->GetNumberOfParameters());

  /* use vnl_vector data_block() to get data ptr */
  std::copy(m_ShapeModelParameters.data_block(), m_ShapeModelParameters.data_block() + m_NumberOfUsedComponents, &(this->m_Parameters.data_block())[0]);

  for (size_t i = m_NumberOfUsedComponents; i < m_NumberOfComponents; ++i) {
    this->m_Parameters[i] = 0;
  }

  /* use vnl_vector data_block() to get data ptr */
  const ParametersType & subParameters = m_SpatialTransform->GetParameters();
  std::copy(subParameters.data_block(), subParameters.data_block() + subParameters.Size(), &(this->m_Parameters.data_block())[m_NumberOfComponents]);

  return this->m_Parameters;
}


template<typename TDataset, typename TParametersValueType>
void StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template<typename TDataset, typename TParametersValueType>
typename StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::OutputPointType
StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::TransformPoint(const size_t & index) const
{
  try {
    return m_SpatialTransform->TransformPoint(m_ShapeModel->DrawSampleAtPoint(m_ShapeModelParameters, index));
  }
  catch (itk::ExceptionObject &excep) {
    std::cout << "exception occurred at point " << index << std::endl;
    itkExceptionMacro(<< excep);
  }
}

template<typename TDataset, typename TParametersValueType>
typename StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::OutputPointType
StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>::TransformPoint(const InputPointType & p) const
{
  try {
    return m_SpatialTransform->TransformPoint(m_ShapeModel->DrawSampleAtPoint(m_ShapeModelParameters, p));
  }
  catch (itk::ExceptionObject &excep) {
    std::cout << "exception occurred at point " << p << std::endl;
    itkExceptionMacro(<< excep);
  }
}


template<typename TDataset, typename TParametersValueType>
void StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>
::ComputeJacobianWithRespectToParameters( const InputPointType & p, JacobianType & outJacobian) const
{
  /* Returns a concatenated MxN array, holding the Jacobian of each sub
  * transform that is selected for optimization. 
  * M rows = dimensionality of the transforms
  * N cols = total number of parameters in the selected sub transforms. */
  outJacobian.SetSize(Dimension, m_NumberOfParameters);
  JacobianType jacobianWithRespectToPosition(Dimension, Dimension);
  this->ComputeJacobianWithRespectToParametersCachedTemporaries(p, outJacobian, jacobianWithRespectToPosition);
}

template<typename TDataset, typename TParametersValueType>
void StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>
::ComputeJacobianWithRespectToParametersCachedTemporaries(const InputPointType & p, JacobianType & outJacobian, JacobianType & jacobianWithRespectToPosition) const
{
  outJacobian.set_size(Dimension, m_NumberOfParameters);
  JacobianType modelJacobian(Dimension, m_NumberOfComponents);
  JacobianType spatialJacobian(Dimension, this->m_SpatialTransform->GetNumberOfParameters());

  const auto & index = m_ShapeModel->GetRepresenter()->GetPointIdForPoint(p);

  this->ComputeJacobianWithRespectToParametersCachedTemporaries(index, outJacobian, modelJacobian, spatialJacobian, jacobianWithRespectToPosition);
}

template<typename TDataset, typename TParametersValueType>
void StatisticalShapeModelMultiTransform<TDataset, TParametersValueType>
::ComputeJacobianWithRespectToParametersCachedTemporaries(const size_t & index, JacobianType & outJacobian, JacobianType & modelJacobian, JacobianType & spatialJacobian, JacobianType & jacobianWithRespectToPosition) const
{
  const auto & jacobian = m_ShapeModel->GetJacobian(index);

  for (size_t i = 0; i < Dimension; ++i) {
    for (size_t k = 0; k < m_NumberOfUsedComponents; ++k) {
      modelJacobian[i][k] = jacobian[i][k];
    }

    for (size_t k = m_NumberOfUsedComponents; k < m_NumberOfComponents; ++k) {
      modelJacobian[i][k] = 0;
    }
  }

  const auto & modelTransformedPoint = m_ShapeModel->DrawSampleAtPoint(m_ShapeModelParameters, index);

  m_SpatialTransform->ComputeJacobianWithRespectToParameters(modelTransformedPoint, spatialJacobian);
  outJacobian.update(spatialJacobian, 0, m_NumberOfComponents);

  m_SpatialTransform->ComputeJacobianWithRespectToPosition(modelTransformedPoint, jacobianWithRespectToPosition);
  outJacobian.update(jacobianWithRespectToPosition * modelJacobian, 0, 0);
}
}
