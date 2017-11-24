#pragma once

#include <itkProcessObject.h>
#include <vtkPolyData.h>

namespace ssm
{
template< typename TInputImage, typename TOutputMesh >
class BinaryMask3DMeshSource:public itk::ProcessObject
{
public:
  /** Standard "Self" typedef. */
  typedef BinaryMask3DMeshSource            Self;
  typedef itk::ProcessObject                Superclass;
  typedef itk::SmartPointer< Self >         Pointer;
  typedef itk::SmartPointer< const Self >   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BinaryMask3DMeshSource, itk::ProcessObject);

  /** Hold on to the type information specified by the template parameters. */
  typedef TOutputMesh                   OutputMeshType;
  typedef vtkSmartPointer<TOutputMesh>  OutputMeshPointer;

  /** Input Image Type Definition. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::SpacingType  SpacingType;
  typedef typename InputImageType::PointType    OriginType;
  typedef typename InputImageType::RegionType   RegionType;
  typedef typename InputImageType::SizeType     SizeType;

  void SetInput(const TInputImage *image);
  OutputMeshType * GetOutput();
  void Update() { this->GenerateData(); };

  itkSetMacro(LevelValue, double);
  itkGetMacro(LevelValue, double);

  itkSetMacro(ComputeLevelValue, bool);
  itkGetMacro(ComputeLevelValue, bool);

  itkSetMacro(Sigma, double);
  itkGetMacro(Sigma, double);

  itkSetMacro(NumberOfIterations, size_t);
  itkGetMacro(NumberOfIterations, size_t);

  itkSetMacro(NumberOfPoints, size_t);
  itkGetMacro(NumberOfPoints, size_t);

  itkSetMacro(RelaxationFactor, double);
  itkGetMacro(RelaxationFactor, double);

  void PrintReport() const;

protected:
  BinaryMask3DMeshSource();
  ~BinaryMask3DMeshSource() {};

  void GenerateData() ITK_OVERRIDE;
  typename TInputImage::ConstPointer GetInput();

private:
  BinaryMask3DMeshSource(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

  typedef typename InputImageType::SizeType InputImageSizeType;

  void CreateMesh();

  vtkSmartPointer<OutputMeshType> m_Output;
  double m_LevelValue;
  bool m_ComputeLevelValue;

  /** temporary variables used in CreateMesh to avoid thousands of
   *  calls to GetInput() and GetOutput()
   */
  OutputMeshPointer m_OutputMesh;

  double m_Sigma;
  size_t m_NumberOfIterations;
  double m_RelaxationFactor;
  double m_Reduction;
  size_t m_NumberOfPoints;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "ssmBinaryMask3DMeshSource.hxx"
#endif
