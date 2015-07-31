#ifndef __itkPeakIntensityFilter_h
#define __itkPeakIntensityFilter_h

#include "itkNeighborhoodOperatorImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

namespace itk
{

template <class TImage, class TLabelImage, class TInterpolator = NearestNeighborInterpolateImageFunction<TImage> >
class ITK_EXPORT PeakIntensityFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef PeakIntensityFilter      Self;
  typedef ProcessObject            Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  /** Useful class typedefs*/
  typedef TImage                           ImageType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::PixelType    PixelType;
  typedef typename ImageType::PointType    PointType;
  typedef typename ImageType::SpacingType  SpacingType;
  typedef typename ImageType::SizeType     SizeType;
  typedef typename ImageType::IndexType    IndexType;

  typedef TLabelImage                           LabelImageType;
  typedef typename LabelImageType::Pointer      LabelImagePointer;
  typedef typename LabelImageType::ConstPointer LabelImageConstPointer;
  typedef typename LabelImageType::PixelType    LabelPixelType;
  itkNewMacro( Self );
  
  /** Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int, ImageType::ImageDimension);
  //typedef typename itk::Image<double, ImageType::ImageDimension> InternalImageType;
  typedef typename itk::Neighborhood<double, ImageType::ImageDimension> NeighborhoodType;
  typedef typename itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef typename itk::NeighborhoodOperatorImageFunction<ImageType, double> NeighborhoodOperatorImageFunctionType;
  typedef typename itk::NeighborhoodOperatorImageFunction<LabelImageType, int> LabelNeighborhoodOperatorImageFunctionType;
  

  /** Run-time type information (and related methods). */
  itkTypeMacro(PeakIntensityFilter, ProcessObject);
  
  //Set functions for the images, in const and non-const form
  void SetInputImage( const ImageType* input );
  void SetInputImage( ImageType* input );
  void SetInputLabelImage( const LabelImageType* input );
  void SetInputLabelImage( LabelImageType* input );

  ImagePointer GetInputImage() const;
  LabelImagePointer GetInputLabelImage() const;

  //Set and Get macros for the various values
  itkSetMacro(CurrentLabel, LabelPixelType);
  itkGetMacro(CurrentLabel, LabelPixelType);

  itkGetMacro(MaximumValue, double);
  itkGetMacro(AverageValue, double);
  itkGetMacro(MinimumValue, double);
  itkGetMacro(PeakValue, double);
  itkGetMacro(PeakIndex, IndexType);
  itkGetMacro(SegmentedVolume, double);
  itkGetMacro(PeakLocation, PointType);
  itkGetMacro(NaivePeak, double);
  itkGetMacro(NaivePeakIndex, IndexType);
  itkGetMacro(NaivePeakLocation, PointType);
  
  itkSetMacro(SphereSpacing, SpacingType);
  itkGetMacro(SphereSpacing, SpacingType);
  itkSetMacro(SphereRadius, PointType);
  void SetSphereRadius(float r);
  itkGetMacro(SphereRadius, PointType);
  //itkSetMacro(SphereVolume, double);
  void SetSphereVolume(double volume);
  itkGetMacro(SphereVolume, double);
  itkSetMacro(SamplingFactor, int);
  itkSetMacro(UseInteriorOnly, bool);
  itkSetMacro(UseSourceSpacing, bool);

  void CalculatePeak();
  double GetKernelVolume();

protected:
  PeakIntensityFilter();
  ~PeakIntensityFilter();
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData();
  void BuildPeakKernel();
  void BuildIsotropicKernel();
  void BuildKernel();
  void MakeKernelOperators( NeighborhoodOperatorImageFunctionType* neighborhoodOperator,
                            LabelNeighborhoodOperatorImageFunctionType* labelNeighborhoodOperator );
  
  typedef typename itk::Image<double, ImageType::ImageDimension> InternalImageType;
  void ExtractLabelRegion();
  void CalculateSphereRadius();
  
private:
  PeakIntensityFilter(const PeakIntensityFilter&); //purposely not implemented
  void operator=(const PeakIntensityFilter&); //purposely not implemented

  /** The label used to calculate indices. */
  LabelPixelType m_CurrentLabel;
  /** The maximum segmented value.  */
  double m_MaximumValue;
  /** The average segmented value.  */
  double m_AverageValue;
  /** The minimum segmented value.  */
  double m_MinimumValue;
  /** The peak segmented value.  */
  double m_PeakValue;
  /** The index of the peak.  */
  IndexType m_PeakIndex;
  /** The location of the peak.  */
  PointType m_PeakLocation;
  /** The segmented volume.  */
  float m_SegmentedVolume;
  /** Naive peak value. */
  double m_NaivePeak;
  /** The location of the naive peak.  */
  PointType m_NaivePeakLocation;
  
  SpacingType m_SphereSpacing;
  PointType m_SphereRadius;
  double m_SphereVolume;
  double m_KernelSum;
  int m_MaskCount;
  SizeType m_KernelRadius;
  typename NeighborhoodOperatorImageFunctionType::Pointer m_PeakKernelOperator;
  typename LabelNeighborhoodOperatorImageFunctionType::Pointer m_MaskKernelOperator;
  int m_SamplingFactor;
  typename ImageType::Pointer m_CroppedInputImage;
  typename LabelImageType::Pointer m_CroppedLabelImage;
  bool m_UseInteriorOnly;
  IndexType m_NaivePeakIndex;
  bool m_UseSourceSpacing;
  typename InternalImageType::Pointer m_KernelImage;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPeakIntensityFilter.cxx"
#endif

#endif
