#pragma once

#include <itkImage.h>
#include <itkMesh.h>
#include <itkStatisticalModel.h>
#include <itkStandardMeshRepresenter.h>

const unsigned int Dimension = 3;

typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;

typedef itk::Point<float, Dimension> PointType;
typedef itk::PointSet<PointType, Dimension> PointSetType;
typedef itk::Mesh<float, Dimension> MeshType;

typedef itk::StatisticalModel<MeshType> ShapeModelType;
typedef itk::StandardMeshRepresenter<float, Dimension> RepresenterType;
typedef std::vector<MeshType::Pointer> MeshVectorType;

typedef std::vector<std::string> StringVector;

