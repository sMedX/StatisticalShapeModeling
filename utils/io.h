#ifndef __utils_h
#define __utils_h

#include <string>
#include <iostream>

#include <itkObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMeshFileWriter.h>
#include <itkMeshFileReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include "itkTransformFileWriter.h"

//! Reads a templated image from a file via ITK ImageFileReader
template <typename TImage>
bool readImage(typename TImage::Pointer image, const std::string& fileName)
{
  typedef itk::ImageFileReader<TImage> Reader;

  typename Reader::Pointer reader = Reader::New();

  reader->SetFileName(fileName);

  try {
    reader->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to read image from file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  image->Graft(reader->GetOutput());
  return true;
}

//! Writes a templated image to a file via ITK ImageFileWriter
template <typename TImage>
bool writeImage(const TImage* image, const std::string& fileName)
{
  typedef itk::ImageFileWriter<TImage> Writer;

  typename Writer::Pointer writer = Writer::New();

  writer->SetInput(image);
  writer->SetFileName(fileName);
  writer->SetUseCompression(true);

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to write image to file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  return true;
}

//! Reads a mesh from a file
template <typename TMesh>
bool readMesh(typename TMesh::Pointer mesh, const std::string& fileName)
{
  typedef itk::MeshFileReader<TMesh> MeshFileReader;
  typename MeshFileReader::Pointer reader = MeshFileReader::New();
  reader->SetFileName(fileName);

  try {
    reader->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to write mesh to file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  mesh->Graft(reader->GetOutput());
  return true;
}

//! Writes a mesh to a file
template <typename TMesh>
bool writeMesh(const TMesh* mesh, const std::string& fileName)
{
  typedef itk::MeshFileWriter<TMesh> MeshFileWriter;
  typename MeshFileWriter::Pointer writer = MeshFileWriter::New();

  writer->SetFileName(fileName);
  writer->SetInput(mesh);
  writer->SetUseCompression(true);
  writer->SetFileTypeAsBINARY();

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to write mesh to file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  return true;
}

bool writeVTKPolydata(vtkPolyData* surface, const std::string& fileName)
{
  typedef vtkSmartPointer<vtkPolyDataWriter> Writer;

  Writer writer = Writer::New();
  writer->SetInputData(surface);
  writer->SetFileName(fileName.c_str());
  writer->SetFileTypeToBinary();
  bool result = static_cast<bool>(writer->Write());

  if (!result) {
    std::cerr << "Error: Unable to write surface to file '" << fileName << "'" << std::endl;
  }

  return result;
}

//! Writes a mesh to a file
template <typename TransformType>
bool writeTransform(const TransformType* transform, const std::string& fileName)
{
  typedef TransformType::ScalarType ScalarType;

  itk::TransformFileWriterTemplate<ScalarType>::Pointer writer = itk::TransformFileWriterTemplate<ScalarType>::New();
  writer->SetInput(transform);
  writer->SetFileName(fileName);

  try {
    writer->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to write transform to file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }

  return true;
}


#endif
