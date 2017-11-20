#pragma once

#include <itkImage.h>
#include <itkMesh.h>

const unsigned int Dimension = 3;
typedef itk::Image<unsigned char, Dimension> BinaryImageType;
typedef itk::Image<float, Dimension> FloatImageType;
typedef itk::Mesh<float, Dimension> MeshType;

typedef std::list<std::string> StringList;

