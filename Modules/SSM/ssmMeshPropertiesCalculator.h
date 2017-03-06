#pragma once

#include "itkAffineTransform.h"
#include "itkImage.h"
#include "itkSpatialObject.h"

#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_diag_matrix.h"

namespace ssm
{
template< typename TMesh >
class MeshPropertiesCalculator:public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef MeshPropertiesCalculator< TMesh >      Self;
  typedef itk::Object                         Superclass;
  typedef itk::SmartPointer< Self >           Pointer;
  typedef itk::SmartPointer< const Self >     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshPropertiesCalculator, itk::Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(Dimension, unsigned int, TMesh::PointDimension);

  /** Standard scalar type within this class. */
  typedef double ScalarType;

  /** Standard vector type within this class. */
  typedef itk::Vector< ScalarType, Dimension> VectorType;

  /** Spatial Object type within this class. */
  typedef itk::SpatialObject<Dimension> SpatialObjectType;

  /** Spatial Object member types used within this class. */
  typedef typename SpatialObjectType::Pointer       SpatialObjectPointer;
  typedef typename SpatialObjectType::ConstPointer  SpatialObjectConstPointer;

  /** Standard matrix type within this class. */
  typedef itk::Matrix< ScalarType, Dimension, Dimension>   MatrixType;

  /** Standard image type and pointer within this class. */
  typedef TMesh MeshType;
  typedef typename MeshType::ConstPointer          MeshConstPointer;
  typedef itk::Image<unsigned char, Dimension>     BinaryImageType;

  /** Affine transform for mapping to and from principal axis */
  typedef itk::AffineTransform< double, itkGetStaticConstMacro(Dimension) > AffineTransformType;
  typedef typename AffineTransformType::Pointer                             AffineTransformPointer;

  /** Set the input image. */
  virtual void SetMesh(const MeshType *mesh)
  {
    if ( m_Mesh != mesh ) {
      m_Mesh = mesh;
      this->Modified();
      m_Valid = false;
    }
  }

  /** Set the spatial object mask. */
  virtual void SetSpatialObjectMask(const itk::SpatialObject< itkGetStaticConstMacro(Dimension) > *so)
  {
    if ( m_SpatialObjectMask != so )
      {
      m_SpatialObjectMask = so;
      this->Modified();
      m_Valid = false;
      }
  }

  void Compute();

  /** Get center of mask gravity, in physical coordinates.*/
  VectorType GetCenterOfMaskGravity() const;

  /** Get radius.*/
  ScalarType GetRadius() const;

protected:
  MeshPropertiesCalculator();
  virtual ~MeshPropertiesCalculator() {};
  virtual void PrintSelf(std::ostream & os, itk::Indent indent) const ITK_OVERRIDE;

private:
  MeshPropertiesCalculator(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  bool       m_Valid;
  VectorType m_CenterOfMaskGravity;
  ScalarType m_Radius;

  MeshConstPointer                     m_Mesh;
  SpatialObjectConstPointer            m_SpatialObjectMask;
  typename BinaryImageType::Pointer    m_Mask;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmMeshPropertiesCalculator.hxx"
#endif
