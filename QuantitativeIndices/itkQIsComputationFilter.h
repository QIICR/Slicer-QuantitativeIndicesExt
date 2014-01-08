#ifndef __itkQIsComputationFilter_h
#define __itkQIsComputationFilter_h

#include "itkMeshSource.h"

namespace itk
{

template <class TImage, class TLabelImage>
class ITK_EXPORT QIsComputationFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef QIsComputationFilter     Self;
  typedef ProcessObject  Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;
  
  typedef TImage                  ImageType;
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;
  typedef typename ImageType::PixelType     PixelType;
  

  typedef TLabelImage                  LabelImageType;
  typedef typename LabelImageType::Pointer      LabelImagePointer;
  typedef typename LabelImageType::ConstPointer LabelImageConstPointer;
  typedef typename LabelImageType::PixelType     LabelType;
  itkNewMacro( Self );

  typedef typename ImageType::PointType PointType;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(QIsComputationFilter, ProcessObject);
  
  //Set functions for the images, in const and non-const form
  void SetInputImage( const ImageType* input );
  void SetInputImage( ImageType* input );
  void SetInputLabelImage( const LabelImageType* input );
  void SetInputLabelImage( LabelImageType* input );

  ImagePointer GetInputImage() const;
  LabelImagePointer GetInputLabelImage() const;

  //Set and Get macros for the various values
  itkSetMacro(CurrentLabel, LabelType);
  itkGetMacro(CurrentLabel, LabelType);

  itkGetMacro(MaximumValue, PixelType);
  itkGetMacro(AverageValue, PixelType);
  itkGetMacro(RMSValue, PixelType);
  itkGetMacro(MinimumValue, PixelType);
  itkGetMacro(PeakValue, PixelType);
  itkGetMacro(MetabolicVolume, PixelType);
  itkGetMacro(SegmentedVolume, float);
  itkGetMacro(PeakLocation, typename ImageType::PointType);
  itkGetMacro(MedianValue, PixelType);
  itkGetMacro(FirstQuartileValue, PixelType);
  itkGetMacro(ThirdQuartileValue, PixelType);
  itkGetMacro(UpperAdjacentValue, PixelType);
  itkGetMacro(StandardDeviation, PixelType);
  itkGetMacro(SAMValue, PixelType);

  void CalculateMean();
  void CalculateQuartiles();
  void CalculatePeak();
  void CalculateSAM();

protected:
  QIsComputationFilter();
  ~QIsComputationFilter();
  virtual void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData();

  
private:
  QIsComputationFilter(const QIsComputationFilter&); //purposely not implemented
  void operator=(const QIsComputationFilter&); //purposely not implemented

  /** The label to calculate indices for. */
  LabelType m_CurrentLabel;

  /** The maximum segmented value.  */
  PixelType m_MaximumValue;
  /** The average segmented value.  */
  PixelType m_AverageValue;
  /** The root-mean-square segmented value */
  PixelType m_RMSValue;
  /** The median segmented value.  */
  PixelType m_MedianValue;
  /** The minimum segmented value.  */
  PixelType m_MinimumValue;
  /** The peak segmented value.  */
  PixelType m_PeakValue;
  /** The location of the peak.  */
  PointType m_PeakLocation;
  /** The segmented volume.  */
  float m_SegmentedVolume;
  /** The metabolic volume segmented.  */
  PixelType m_MetabolicVolume;
  /** The first quartile segmented value.  */
  PixelType m_FirstQuartileValue;
  /** The third quartile segmented value.  */
  PixelType m_ThirdQuartileValue;
  /** The upper adjacent value */
  PixelType m_UpperAdjacentValue;
  /** The standard deviation of the segmented values. */
  PixelType m_StandardDeviation;
  /** The standard added metabolic activity. */
  PixelType m_SAMValue;

  /** Flag indicating if list has been generated */
  bool m_ListGenerated;
  /** List of values in region of interest */
  std::list<double> m_SegmentedValues;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQIsComputationFilter.txx"
#endif

#endif
