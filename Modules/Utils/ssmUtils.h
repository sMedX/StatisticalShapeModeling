#pragma once
#include <boost/filesystem.hpp>

#include <itkObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMeshFileReader.h>
#include <itkMeshFileWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <itkMesh.h>

#include "ssmTypes.h"

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

//! print information about image
template <typename TImage>
void printImageInfo(const TImage* image, const std::string &info = "")
{
  if (info.size() > 0) {
    std::cout << info << std::endl;
  }
  std::cout << "    size " << image->GetLargestPossibleRegion().GetSize() << ", " << image->GetNumberOfComponentsPerPixel() << std::endl;
  std::cout << "  origin " << image->GetOrigin() << std::endl;
  std::cout << "spacing  " << image->GetSpacing() << std::endl;
  std::cout << "direction" << std::endl << image->GetDirection() << std::endl;
}

//! print information about mesh
template <typename TMesh>
void printMeshInfo(const TMesh* surface, const std::string &info = "")
{
  if (info.size() > 0) {
    std::cout << info << std::endl;
  }
  std::cout << "number of cells  " << surface->GetNumberOfCells() << std::endl;
  std::cout << "number of points " << surface->GetNumberOfPoints() << std::endl;
  std::cout << std::endl;
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
    std::cerr << "Unable to read mesh to file '" << fileName << "'" << std::endl;
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

bool readVTKPolydata(vtkPolyData* surface, const std::string& filename)
{
  typedef vtkSmartPointer<vtkPolyDataReader> Reader;
  Reader reader = Reader::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  surface->ShallowCopy(reader->GetOutput());
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

//! Writes a transform to a file
template <typename TransformType>
bool writeTransform(const TransformType* transform, const std::string& fileName)
{
  typedef typename TransformType::ScalarType ScalarType;

  typename itk::TransformFileWriterTemplate<ScalarType>::Pointer writer = itk::TransformFileWriterTemplate<ScalarType>::New();
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

//! Reads a transform from a file
typedef itk::TransformFileReader::TransformListType * TransformListType;
bool readTransform(TransformListType transforms, const std::string& fileName)
{
  itk::TransformFactoryBase::RegisterDefaultTransforms();
  itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
  reader->SetFileName(fileName);
  try {
    reader->Update();
  }
  catch (itk::ExceptionObject& err) {
    std::cerr << "Unable to read transform from file '" << fileName << "'" << std::endl;
    std::cerr << "Error: " << err << std::endl;
    return false;
  }
  transforms = reader->GetTransformList();
  return true;
}

std::string getDirectoryFromPath(const std::string& fileName)
{
  boost::filesystem::path path(fileName);
  return path.parent_path().string();
}

std::string getFileNameFromPath(const std::string& fileName)
{
  boost::filesystem::path path(fileName);
  return path.filename().string();
}

std::string addFileNameSuffix(const std::string& fileName, const::std::string& suffix)
{
  boost::filesystem::path path(fileName);
  path = path.parent_path() / boost::filesystem::path(path.stem().string() + suffix + path.extension().string());
  return path.string();
}

std::string getBaseNameFromPath(const std::string& fileName)
{
  boost::filesystem::path path(fileName);
  return path.stem().string();
}

StringVector readListFromFile(const std::string& fileName)
{
  StringVector list;

  std::ifstream file;
  try {
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    file.open(fileName.c_str(), std::ifstream::in);
    std::string line;
    while (getline(file, line)) {
      if (line != "") {
        //reading files with windows EOL on Linux results in the \r not being removed from the line ending
        if (*line.rbegin() == '\r') {
          line.erase(line.length() - 1, 1);
        }
        list.push_back(line);
      }
    }
  }
  catch (std::ifstream::failure e) {
    if (file.eof() == false) {
      throw std::ifstream::failure("Failed to read list from the file: " + fileName);
    }
  }

  return list;
}

void writeListToFile(const std::string & fileName, const StringVector & list)
{
  std::ofstream file(fileName, std::ofstream::out);
  if (!file.is_open()) {
    throw std::ofstream::failure("Failed to write list to the file : " + fileName);
  }

  for (const auto & string : list) {
    file << string << std::endl;
  }
  file.close();

  return;
}

bool checkFileName(const std::string & fileName)
{
  std::ofstream file(fileName, std::ofstream::out);
  if (!file.is_open()) {
    std::cerr << "Failed to open file " << fileName << std::endl;
    return false;
  }

  file.close();
  boost::filesystem::remove(fileName);

  return true;
}
