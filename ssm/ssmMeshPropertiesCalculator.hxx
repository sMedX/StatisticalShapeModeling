#pragma once

#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkTriangleMeshToBinaryImageFilter.h>

#include "ssmMeshPropertiesCalculator.h"

namespace ssm
{
class InvalidMeshMomentsError:public itk::ExceptionObject
{
public:
  /**
   * Constructor. Needed to ensure the exception object can be copied.
   */
  InvalidMeshMomentsError(const char *file, unsigned int lineNumber):ExceptionObject(file,
                                                                                      lineNumber) { this->
                                                                                                    SetDescription(
                                                                                                      "No valid image moments are available."); }

  /**
   * Constructor. Needed to ensure the exception object can be copied.
   */
  InvalidMeshMomentsError(const std::string & file, unsigned int lineNumber):ExceptionObject(file,
                                                                                              lineNumber) { this->
                                                                                                            SetDescription(
                                                                                                              "No valid image moments are available."); }

  itkTypeMacro(InvalidMeshMomentsError, ExceptionObject);
};

//----------------------------------------------------------------------
// Construct without computing moments
template< typename TMesh >
MeshPropertiesCalculator< TMesh >::MeshPropertiesCalculator(void)
{
  m_Valid = false;
  m_Mesh = ITK_NULLPTR;
  m_SpatialObjectMask = ITK_NULLPTR;
}

//----------------------------------------------------------------------
template< typename TInputImage >
void MeshPropertiesCalculator< TInputImage >::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Image: " << m_Mesh.GetPointer() << std::endl;
  os << indent << "Valid: " << m_Valid << std::endl;
}

//----------------------------------------------------------------------
// Compute moments for a new or modified image
template< typename TMesh >
void MeshPropertiesCalculator< TMesh >::Compute()
{
  if ( m_Mesh==nullptr ) {
    return;
  }

  // compute binary mask
  // Compute a bounding box of the input mesh
  typename MeshType::BoundingBoxType::ConstPointer boundingBox = m_Mesh->GetBoundingBox();
  typename BinaryImageType::SpacingType spacing(1);
  typename BinaryImageType::PointType origin = boundingBox->GetMinimum();
  typename BinaryImageType::SizeType size;

  for (size_t n = 0; n < Dimension; ++n) {
    size[n] = (boundingBox->GetMaximum()[n] - boundingBox->GetMinimum()[n]) / spacing[n];
  }

  typedef itk::TriangleMeshToBinaryImageFilter<MeshType, BinaryImageType> ShapeToBinaryImageFilterType;
  typename ShapeToBinaryImageFilterType::Pointer shapeToImage = ShapeToBinaryImageFilterType::New();
  shapeToImage->SetInput(const_cast<MeshType*>(m_Mesh.GetPointer()));
  shapeToImage->SetSize(size);
  shapeToImage->SetOrigin(origin);
  shapeToImage->SetSpacing(spacing);
  shapeToImage->SetOutsideValue(0);
  shapeToImage->SetInsideValue(1);
  try {
    shapeToImage->Update();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }
  m_Mask = shapeToImage->GetOutput();

  // moment image calculator
  typedef itk::ImageMomentsCalculator<BinaryImageType>  ImageCalculatorType;
  typename ImageCalculatorType::Pointer m_ImageCalculator = ImageCalculatorType::New();
  m_ImageCalculator->SetImage(m_Mask);
  try {
    m_ImageCalculator->Compute();
  }
  catch (itk::ExceptionObject& excep) {
    std::cout << excep << std::endl;
    itkExceptionMacro(<< excep);
  }

  // set center of mask 
  m_CenterOfMaskGravity = m_ImageCalculator->GetCenterOfGravity();

  // compute radius
  m_Radius = itk::NumericTraits< ScalarType >::ZeroValue();

  for (auto it = m_Mesh->GetPoints()->Begin(); it != m_Mesh->GetPoints()->End(); ++it) {
    typename MeshType::PointType point = it.Value();
    VectorType vector;
    for (size_t i = 0; i < Dimension; ++i) {
      vector[i] = point[i] - m_CenterOfMaskGravity[i];
    }
    m_Radius += vector.GetNorm();
  }

  m_Radius /= m_Mesh->GetNumberOfPoints();

  /* Remember that the moments are valid */
  m_Valid = 1;
}

//--------------------------------------------------------------------
// Get center of mask gravity, in physical coordinates
template< typename TMesh >
typename MeshPropertiesCalculator< TMesh >::VectorType MeshPropertiesCalculator< TMesh >::GetCenterOfMaskGravity() const
{
  if (!m_Valid) {
    itkExceptionMacro(<< "GetCenterOfGravity() invoked, but the moments have not been computed. Call Compute() first.");
  }
  return m_CenterOfMaskGravity;
}
//--------------------------------------------------------------------
// Get center of radius
template< typename TMesh >
typename MeshPropertiesCalculator< TMesh >::ScalarType MeshPropertiesCalculator< TMesh >::GetRadius() const
{
  if (!m_Valid) {
    itkExceptionMacro(<< "GetCenterOfGravity() invoked, but the moments have not been computed. Call Compute() first.");
  }
  return m_Radius;
}
}
