#ifndef _itkPeakIntensityFilter_cxx
#define _itkPeakIntensityFilter_cxx

#include "itkPeakIntensityFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <itkResampleImageFilter.h>
#include <itkConstShapedNeighborhoodIterator.h>
#include <itkRegionOfInterestImageFilter.h>
#include "itkBinaryBallStructuringElement.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkExtractImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include <math.h>

#define QI_PEAK_RADIUS 6.204//2.5
#define QI_PEAK_RADIUS_SPACING_RATE 4.0
#define PI 3.14159265359

#define INITIAL_MARGIN 3

#include "itkImageFileWriter.h"


namespace itk
{

//----------------------------------------------------------------------------
template <class TImage, class TLabelImage, class TInterpolator>
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::PeakIntensityFilter()
{
  this->ProcessObject::SetNumberOfRequiredInputs(5);
  this->ProcessObject::SetNumberOfRequiredOutputs(0);
  SpacingType sphereSpacing;
  sphereSpacing.Fill(1);
  m_SphereSpacing = sphereSpacing;
  m_SphereVolume = 1000;  // 1cc sphere
  //this->CalculateSphereRadius();
  double r = std::pow(1000*0.75/PI,1.0/3.0);
  m_SphereRadius.Fill(r); // approx. 1cc sphere
  m_UseInteriorOnly = true;
  m_UseSourceSpacing = false;
}

//----------------------------------------------------------------------------
template <class TImage, class TLabelImage, class TInterpolator>
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::~PeakIntensityFilter()
{}


//----------------------------------------------------------------------------
/*
SetSphereRadius
Sets the radii for the sphere.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetSphereVolume(double volume)
{
  m_SphereVolume = volume;
  this->CalculateSphereRadius();
}


//----------------------------------------------------------------------------
/*
SetSphereRadius
Sets the radii for the sphere.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetSphereRadius(float r)
{
  m_SphereRadius.Fill(r);
}


//----------------------------------------------------------------------------
/*
CalculateSphereRadius
Determines the radius based on sphere volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::CalculateSphereRadius()
{
  double r = pow(m_SphereVolume*0.75/PI,1.0/3.0);
  m_SphereRadius.Fill(r);
}

//----------------------------------------------------------------------------
/*
SetInputImage
Sets the input volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetInputImage( const ImageType* input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, input);
}

//----------------------------------------------------------------------------
/*
SetInputImage
Sets the input volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetInputImage( ImageType* input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, input);
}


//----------------------------------------------------------------------------
/*
GetInputImage
Returns the input image volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
typename PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::ImagePointer
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::GetInputImage() const
{
  if ( this->GetNumberOfInputs() < 1 )
  { return NULL;  }
  else
  {
    ImagePointer image = ImageType::New();
    image->Graft(this->ProcessObject::GetInput(0));
    return image;
  }
}


//----------------------------------------------------------------------------
/*
SetInputLabelImage
Sets the label volume for the iput volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetInputLabelImage( const LabelImageType* input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(1, input);//this->ProcessObject::SetNthInput( 0, const_cast< HistogramType * >( input ) );
}

//----------------------------------------------------------------------------
/*
SetInputLabelImage
Sets the label volume for the input volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::SetInputLabelImage( LabelImageType*  input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(1, /*const_cast< const LabelImageType* >*/(input));//this->ProcessObject::SetNthInput( 0, const_cast< HistogramType * >( input ) );
}


//----------------------------------------------------------------------------
/*
GetInputLabelImage
Returns the label image volume.

*/
template <class TImage, class TLabelImage, class TInterpolator>
typename PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::LabelImagePointer
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::GetInputLabelImage() const
{
  if ( this->GetNumberOfInputs() < 2 )
  { return NULL;  }
  else
  {
    LabelImagePointer labelImage = LabelImageType::New();
    labelImage->Graft(this->ProcessObject::GetInput(1));
    return labelImage;
  }
  //{ return static_cast< LabelImageConstPointer >( this->ProcessObject::GetInput(1) );  }
}

//----------------------------------------------------------------------------
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::GenerateData()
{
  std::cout << "    GenerateData()\n";
  this->BuildPeakKernel();
  this->CalculatePeak();
  this->CalculateNaivePeak();
}


//----------------------------------------------------------------------------
/*
ExtractLabelRegion
Crops the input images to the area near the specified label.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::ExtractLabelRegion()
{
std::cout << "  ExtractLabelRegion()\n";
  ImagePointer inputImage = this->GetInputImage();
  LabelImagePointer inputLabel = this->GetInputLabelImage();
  m_CroppedInputImage = ImageType::New();
  m_CroppedLabelImage = LabelImageType::New();
  
  // determine extent of label
  typedef typename itk::ImageRegionConstIteratorWithIndex<LabelImageType> LabelIteratorType;
  LabelIteratorType lit(inputLabel,inputLabel->GetLargestPossibleRegion());
  lit.GoToBegin();
  IndexType lowerIndex;
  IndexType upperIndex;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    lowerIndex[i] = itk::NumericTraits<int>::max();
    upperIndex[i] = itk::NumericTraits<int>::min();
  }
  bool labelFound = false;
  while(!lit.IsAtEnd())
  {
    if(lit.Get() == m_CurrentLabel)
    {
      labelFound = true;
      typename LabelImageType::IndexType idx = lit.GetIndex();
      for(unsigned int i=0; i<ImageDimension; ++i)
      {
        if(idx[i]<lowerIndex[i]) lowerIndex[i]=idx[i];
        if(idx[i]>upperIndex[i]) upperIndex[i]=idx[i];
      }
    }
    ++lit;
  }
  if(!labelFound)
  {
    m_PeakValue = std::numeric_limits<double>::quiet_NaN();
    m_CroppedInputImage = NULL;
    m_CroppedLabelImage = NULL;
    return;
  }
  
  // define new region (pad region by 1.5x radius)
  SpacingType voxelSize = inputImage->GetSpacing();
  SizeType imageSize = inputImage->GetLargestPossibleRegion().GetSize();
  IndexType imageIndex = inputImage->GetLargestPossibleRegion().GetIndex();
  SizeType pad;
  SizeType size;
  IndexType minIndex;
  IndexType maxIndex;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    pad[i] = ceil(1.5*m_SphereRadius[i]/voxelSize[i]);
    minIndex[i] = std::floor(lowerIndex[i]-pad[i]);
    if(minIndex[i]<imageIndex[i]) minIndex[i]=imageIndex[i];
    maxIndex[i] = std::floor(upperIndex[i]+pad[i]);
    if(maxIndex[i]>imageSize[i]-1) maxIndex[i]=imageSize[i];
    size[i] = maxIndex[i]-minIndex[i]+1;
  }
  
  typename ImageType::RegionType region(minIndex, size);
  
  typedef typename itk::ExtractImageFilter<ImageType,ImageType> ImageExtractorType;
  typedef typename itk::ExtractImageFilter<LabelImageType,LabelImageType> LabelExtractorType;
  typename ImageExtractorType::Pointer imageExtractor = ImageExtractorType::New();
  typename LabelExtractorType::Pointer labelExtractor = LabelExtractorType::New();
  imageExtractor->SetInput(inputImage);
  labelExtractor->SetInput(inputLabel);
  imageExtractor->SetExtractionRegion(region);
  labelExtractor->SetExtractionRegion(region);
#if ITK_VERSION_MAJOR >= 4 // This is required.
  imageExtractor->SetDirectionCollapseToIdentity();
  labelExtractor->SetDirectionCollapseToIdentity();
#endif
  imageExtractor->Update();
  labelExtractor->Update();
  
  m_CroppedInputImage = imageExtractor->GetOutput();
  m_CroppedLabelImage = labelExtractor->GetOutput();

}


//----------------------------------------------------------------------------
/*
BuildKernel
Creates the NeighborhoodOperatorImageFunction for the peak kernel.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::BuildKernel()
{
std::cout << "  BuildKernel()\n";

}


//----------------------------------------------------------------------------
/*
BuildPeakKernel
Creates the NeighborhoodOperatorImageFunction for the peak kernel.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::BuildPeakKernel()
{
std::cout << "  BuildPeakKernel()\n";
  ImagePointer inputImage = this->GetInputImage();
  LabelImagePointer labelImage = this->GetInputLabelImage();
  SpacingType voxelSize = inputImage->GetSpacing();
  
  // build a higher-resolution image of the kernel
  SpacingType upsampledKernelSpacing;
  SizeType upsampledKernelSize;
  double voxelFraction = 1.0;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    upsampledKernelSpacing[i] = voxelSize[i]/m_SamplingFactor;
    upsampledKernelSize[i] = (std::ceil(m_SphereRadius[i]/voxelSize[i])*2+1)*m_SamplingFactor;
    voxelFraction /= m_SamplingFactor; // volume fraction of upsampled voxel
  }
  PointType origin; origin.Fill(0);
  IndexType startIndex; startIndex.Fill(0);
  typename InternalImageType::RegionType upsampledKernelRegion(startIndex, upsampledKernelSize);
  typename InternalImageType::Pointer upsampledKernel = InternalImageType::New();
  upsampledKernel->SetOrigin(origin);
  upsampledKernel->SetRegions(upsampledKernelRegion);
  upsampledKernel->SetSpacing(upsampledKernelSpacing);
  upsampledKernel->Allocate();
  
  itk::ContinuousIndex<double,ImageDimension> centerIndex;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    centerIndex[i] = (upsampledKernelSize[i]-1)*0.5;
  }
  PointType centerPoint;
  upsampledKernel->TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

  typedef typename itk::ImageRegionIteratorWithIndex<InternalImageType> KernelIteratorType;
  KernelIteratorType ukit(upsampledKernel,upsampledKernelRegion);
  ukit.GoToBegin();
  while(!ukit.IsAtEnd())
  {
    IndexType currentIndex = ukit.GetIndex();
    PointType currentPoint;
    upsampledKernel->TransformIndexToPhysicalPoint(currentIndex, currentPoint);
    double r = 0.0;
    for(unsigned int i=0; i<ImageDimension; ++i)
    {
      double psq = (currentPoint[i]-centerPoint[i])*(currentPoint[i]-centerPoint[i]);
      r += psq/(m_SphereRadius[i]*m_SphereRadius[i]);
    }
    if(r <= 1.0)
    {
      ukit.Set(1.0);
    }
    else{
      ukit.Set(0.0);
    }
    ++ukit;
  }
  
  // build the full-resolution image of the kernel
  SizeType kernelSize;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    kernelSize[i] = ceil(m_SphereRadius[i]/voxelSize[i])*2+1;
    m_KernelRadius[i] = (kernelSize[i]-1)*0.5;
  }

  m_KernelImage = InternalImageType::New();
  m_KernelImage->SetOrigin(origin);
  typename InternalImageType::RegionType kernelRegion(startIndex, kernelSize);
  m_KernelImage->SetRegions(kernelRegion);
  m_KernelImage->SetSpacing(voxelSize);
  m_KernelImage->Allocate();
  m_KernelImage->FillBuffer(0.0);
  KernelIteratorType kit(m_KernelImage, kernelRegion);
  kit.GoToBegin();
  const double inv_samplingFactor = 1/(double)m_SamplingFactor;
  ukit.GoToBegin();
  while(!ukit.IsAtEnd())
  {
    IndexType currentIndex = ukit.GetIndex();
    double val = ukit.Get();
    IndexType newIndex;
    for(unsigned int i=0; i<ImageDimension; ++i)
    {
      newIndex[i] = std::floor((double)currentIndex[i]*inv_samplingFactor);
    }
    kit.SetIndex( newIndex );
    double old_val = kit.Get();
    kit.Set( old_val + val );
    ++ukit;
  }

  // apply the voxel fraction to get the final kernel values
  kit.GoToBegin();
  while(!kit.IsAtEnd())
  {
    kit.Set( kit.Get()*voxelFraction );
    ++kit;
  }
  
}


//----------------------------------------------------------------------------
/*
BuildIsotropicKernel
Creates the NeighborhoodOperatorImageFunction for the peak kernel.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::BuildIsotropicKernel()
{
std::cout << "  BuildIsotropicKernel()\n";
  ImagePointer inputImage = this->GetInputImage();
  LabelImagePointer labelImage = this->GetInputLabelImage();
  SpacingType voxelSize = inputImage->GetSpacing();
  
  // build a higher-resolution image of the kernel
  SpacingType upsampledKernelSpacing;
  upsampledKernelSpacing.Fill(m_SphereSpacing[0]/m_SamplingFactor);
  SizeType upsampledKernelSize;
  upsampledKernelSize.Fill((std::ceil(m_SphereRadius[0]/m_SphereSpacing[0])*2+1)*m_SamplingFactor);
  //for(unsigned int i=0; i<ImageDimension; ++i)
  //{
    //upsampledKernelSpacing[i] = voxelSize[i]/m_SamplingFactor;
    //upsampledKernelSize[i] = (std::ceil(m_SphereRadius[i]/voxelSize[i])*2+1)*m_SamplingFactor;
  //}
  PointType origin; origin.Fill(0);
  IndexType startIndex; startIndex.Fill(0);
  PointType upsampledKernelOrigin;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    upsampledKernelOrigin[i] =  0.5*((std::ceil(m_SphereRadius[0]/voxelSize[i])*2+1)*voxelSize[i] - (upsampledKernelSize[i])*m_SphereSpacing[0]/(double)m_SamplingFactor) + 0.5*m_SphereSpacing[0]/(double)m_SamplingFactor - 0.5*voxelSize[i];
  }
  typename InternalImageType::RegionType upsampledKernelRegion(startIndex, upsampledKernelSize);
  typename InternalImageType::Pointer upsampledKernel = InternalImageType::New();
  upsampledKernel->SetOrigin(upsampledKernelOrigin);
  upsampledKernel->SetRegions(upsampledKernelRegion);
  upsampledKernel->SetSpacing(upsampledKernelSpacing);
  upsampledKernel->Allocate();
  
  // draw a sphere using physical coordinates
  itk::ContinuousIndex<double,ImageDimension> centerIndex;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    centerIndex[i] = (upsampledKernelSize[i]-1)*0.5;
  }
  PointType centerPoint;
  upsampledKernel->TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

  typedef typename itk::ImageRegionIteratorWithIndex<InternalImageType> KernelIteratorType;
  KernelIteratorType ukit(upsampledKernel,upsampledKernelRegion);
  ukit.GoToBegin();
  while(!ukit.IsAtEnd())
  {
    IndexType currentIndex = ukit.GetIndex();
    PointType currentPoint;
    upsampledKernel->TransformIndexToPhysicalPoint(currentIndex, currentPoint);
    double r = 0.0;
    for(unsigned int i=0; i<ImageDimension; ++i)
    {
      double psq = (currentPoint[i]-centerPoint[i])*(currentPoint[i]-centerPoint[i]);
      r += psq/(m_SphereRadius[i]*m_SphereRadius[i]);
    }
    if(r <= 1.0)
    {
      ukit.Set(1.0);
    }
    else{
      ukit.Set(0.0);
    }
    ++ukit;
  }
  
  // build the full-resolution image of the kernel
  SizeType kernelSize;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    kernelSize[i] = ceil(m_SphereRadius[i]/voxelSize[i])*2+1;
    m_KernelRadius[i] = (kernelSize[i]-1)*0.5;
  }
  typename InternalImageType::Pointer kernel = InternalImageType::New();
  kernel->SetOrigin(origin);
  typename InternalImageType::RegionType kernelRegion(startIndex, kernelSize);
  kernel->SetRegions(kernelRegion);
  kernel->SetSpacing(voxelSize);
  kernel->Allocate();
  kernel->FillBuffer(0.0);
  KernelIteratorType kit(kernel, kernelRegion);
  kit.GoToBegin();
  // image to keep track of counts within each kernel voxel
  typename InternalImageType::Pointer weightImage = InternalImageType::New();
  weightImage->SetOrigin(origin);
  weightImage->SetRegions(kernelRegion);
  weightImage->SetSpacing(voxelSize);
  weightImage->Allocate();
  weightImage->FillBuffer(0.0);
  KernelIteratorType wit(weightImage, kernelRegion);
  wit.GoToBegin();
  
  // iterate through the upsampled kernel and increment the full-resolution kernel
  ukit.GoToBegin();
  while(!ukit.IsAtEnd())
  {
    IndexType currentIndex = ukit.GetIndex();
    double val = ukit.Get();
    PointType currentPoint;
    upsampledKernel->TransformIndexToPhysicalPoint(currentIndex, currentPoint);
    IndexType newIndex;
    kernel->TransformPhysicalPointToIndex(currentPoint, newIndex);
    kit.SetIndex( newIndex );
    double old_val = kit.Get();
    kit.Set( old_val + val );
    wit.SetIndex( newIndex );
    wit.Set( wit.Get()+1 ); // may have different counts in each
    ++ukit;
  }
  kit.GoToBegin(); wit.GoToBegin();
  while(!kit.IsAtEnd())
  {
    double old_val = kit.Get();
    kit.Set( old_val*1/wit.Get() );
    ++kit; ++wit;
  }
  
  // create the mask kernel
  typedef typename itk::Neighborhood<LabelPixelType, ImageType::ImageDimension> LabelNeighborhoodType;
  LabelNeighborhoodType labelNeighborhood;
  labelNeighborhood.SetRadius(m_KernelRadius);
  
  // create the convolution kernel
  NeighborhoodType neighborhood;
  neighborhood.SetRadius(m_KernelRadius);
  kit.GoToBegin();
  
  typename NeighborhoodType::Iterator nit = neighborhood.Begin();
  typename LabelNeighborhoodType::Iterator lnit = labelNeighborhood.Begin();
  m_KernelSum = 0.0;
  m_MaskCount = 0.0;
  while(!kit.IsAtEnd())
  {
    double val = kit.Get();
    *nit = val;
    m_KernelSum += val;
    if( val > 0.0 )
    {
      *lnit = 1;
      ++m_MaskCount;
//std::cout << "valid: " << val << std::endl;
    }
    else{
      *lnit = 0;
    }
    ++nit; ++kit; ++lnit;
  }
  
  this->m_MaskKernelOperator = LabelNeighborhoodOperatorImageFunctionType::New();
  this->m_MaskKernelOperator->SetOperator(labelNeighborhood);
  this->m_MaskKernelOperator->SetInputImage(labelImage);
  
  this->m_PeakKernelOperator = NeighborhoodOperatorImageFunctionType::New();
  this->m_PeakKernelOperator->SetOperator(neighborhood);
  this->m_PeakKernelOperator->SetInputImage(inputImage);

}


//----------------------------------------------------------------------------
/*
CalculatePeak
Actually runs the calculations to determine:
Peak, Peak Location

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::CalculatePeak()
{

  if(m_UseSourceSpacing)
  {
    this->BuildPeakKernel();
  }
  else{
    this->BuildIsotropicKernel();
  }
  this->ExtractLabelRegion();
  
  // create the kernel operators
  typename NeighborhoodOperatorImageFunctionType::Pointer peakOperator = NeighborhoodOperatorImageFunctionType::New();
  peakOperator->SetInputImage(this->GetInputImage());
  typename LabelNeighborhoodOperatorImageFunctionType::Pointer maskOperator = LabelNeighborhoodOperatorImageFunctionType::New();
  maskOperator->SetInputImage(this->GetInputLabelImage());
  this->MakeKernelOperators(peakOperator,maskOperator);
  
std::cout << "  CalculatePeak()\n";
  
  // convolve the kernel and evaluate at valid indices
  typedef typename itk::ImageRegionIterator<ImageType> IteratorType;
  IteratorType it(m_CroppedInputImage,m_CroppedInputImage->GetRequestedRegion());
std::cout << "Mask Count: " << m_MaskCount << std::endl;
  typedef typename itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
  LabelIteratorType lit(m_CroppedLabelImage,m_CroppedLabelImage->GetRequestedRegion());
  it.GoToBegin(); lit.GoToBegin(); 
  double peak = itk::NumericTraits<double>::min();
  double max_center_val = itk::NumericTraits<double>::min();
  IndexType peakIndex;
  while(!it.IsAtEnd())
  {
    if(lit.Get() == m_CurrentLabel)
    {
      IndexType currentIndex = lit.GetIndex();
      int labelSum = maskOperator->EvaluateAtIndex(currentIndex);
      if( !m_UseInteriorOnly )
      {
        labelSum = m_MaskCount;
      }
      if( labelSum == m_MaskCount ) // valid kernel placement
      {
        double center_val = it.Get();
        double val = peakOperator->EvaluateAtIndex(currentIndex);
        if( (float)val>(float)peak )
        {
          peak = val;
          max_center_val = center_val;
          peakIndex = currentIndex;
        }
        if((float)val==(float)peak)
        {
          if(center_val > max_center_val)
          {
            max_center_val = center_val;
            peakIndex = currentIndex;
          }
        }
      }
    }
    ++it; ++lit;
  }
  m_PeakValue = peak;
  m_PeakIndex = peakIndex;
std::cout << "Max Center Value: " << max_center_val << std::endl;
std::cout << "Kernel Volume: " << this->GetKernelVolume() << std::endl;

}


//----------------------------------------------------------------------------
/*
MakeKernelOperators
Builds the NeighborhoodOperatorImageFunctions for the peak kernel.
Requires a pointer for the weighted version of the kernel and a pointer for
the binary version of the kernel.

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::MakeKernelOperators( NeighborhoodOperatorImageFunctionType* neighborhoodOperator,
                      LabelNeighborhoodOperatorImageFunctionType* labelNeighborhoodOperator)
{
std::cout << "  MakeKernelOperators()" << std::endl;

  // create the mask kernel
  typedef typename itk::Neighborhood<LabelPixelType, ImageType::ImageDimension> LabelNeighborhoodType;
  LabelNeighborhoodType labelNeighborhood;
  labelNeighborhood.SetRadius(m_KernelRadius);
  
  // create the convolution kernel
  NeighborhoodType neighborhood;
  neighborhood.SetRadius(m_KernelRadius);
  
  typedef typename itk::ImageRegionConstIterator<InternalImageType> KernelIteratorType;
  KernelIteratorType kit(this->m_KernelImage,this->m_KernelImage->GetLargestPossibleRegion());
  kit.GoToBegin();
  
  typename NeighborhoodType::Iterator nit = neighborhood.Begin();
  typename LabelNeighborhoodType::Iterator lnit = labelNeighborhood.Begin();
  double kernelSum = 0.0;
  m_MaskCount = 0.0;
  while(!kit.IsAtEnd())
  {
    double val = kit.Get();
    *nit = val;
    kernelSum += val;
    if( val > 0.0 )
    {
      *lnit = 1;
      ++m_MaskCount;
    }
    else{
      *lnit = 0;
    }
    ++nit; ++kit; ++lnit;
  }
  nit = neighborhood.Begin();
  while( nit != neighborhood.End() )
  {
    double val = *nit;
    *nit = val/kernelSum;
    ++nit;
  }
  
  labelNeighborhoodOperator->SetOperator(labelNeighborhood); 
  neighborhoodOperator->SetOperator(neighborhood);
  
  typedef typename itk::ImageFileWriter<InternalImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "mask.nrrd" );
  writer->SetInput( m_KernelImage );
  writer->Update();

}


//----------------------------------------------------------------------------
/*
GetKernelVolume
Determines the total volume of the kernel based on the voxel size and weights

*/
template <class TImage, class TLabelImage, class TInterpolator>
double
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::GetKernelVolume()
{
  // determine voxel volume
  typename InternalImageType::SpacingType spacing = this->m_KernelImage->GetSpacing();
  double voxelVolume = 1.0;
  for(unsigned int i=0; i<ImageDimension; ++i)
  {
    voxelVolume *= spacing[i];
  }
  
  // sum up the weights of the kernel image
  double weightSum = 0.0;
  typedef typename itk::ImageRegionConstIterator<InternalImageType> KernelIteratorType;
  KernelIteratorType kit(this->m_KernelImage,this->m_KernelImage->GetLargestPossibleRegion());
  kit.GoToBegin();
  while(!kit.IsAtEnd())
  {
    weightSum += kit.Get();
    ++kit;
  }
  
  return voxelVolume*weightSum;
  
}


//----------------------------------------------------------------------------
/*
CalculateNaivePeak
Actually runs the calculations to determine:
NaivePeak, NaivePeak Location

*/
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::CalculateNaivePeak()
{
  // determine bounding box
  typedef typename itk::ImageRegionConstIteratorWithIndex<LabelImageType> LabelIteratorType;
  ImagePointer inputImage = this->GetInputImage();
  LabelImagePointer inputLabel = this->GetInputLabelImage();

  LabelIteratorType lit(inputLabel,inputLabel->GetLargestPossibleRegion());
  lit.GoToBegin();
  typename LabelImageType::IndexType lowerIndex;
  int xmin = itk::NumericTraits<int>::max(); int xmax = itk::NumericTraits<int>::min();
  int ymin = itk::NumericTraits<int>::max(); int ymax = itk::NumericTraits<int>::min();
  int zmin = itk::NumericTraits<int>::max(); int zmax = itk::NumericTraits<int>::min();
  bool labelFound = false;
  while(!lit.IsAtEnd())
  {
    if(lit.Get()==m_CurrentLabel)
    {
      labelFound = true;
      typename LabelImageType::IndexType idx = lit.GetIndex();
      if(idx[0]<xmin) xmin=idx[0];
      if(idx[0]>xmax) xmax=idx[0];
      if(idx[1]<ymin) ymin=idx[1];
      if(idx[1]>ymax) ymax=idx[1];
      if(idx[2]<zmin) zmin=idx[2];
      if(idx[2]>zmax) zmax=idx[2];
    }
    ++lit;
  }
  if(!labelFound)
  {
    m_NaivePeak = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  
  // pad region by approx 1 cm
  typename LabelImageType::SpacingType voxelSize = inputLabel->GetSpacing();
  int xpad = std::ceil(10/voxelSize[0]);
  int ypad = std::ceil(10/voxelSize[1]);
  int zpad = std::ceil(10/voxelSize[2]);
  
  // resample regions to 1mm isotropic
  typedef typename itk::ResampleImageFilter<ImageType,ImageType> ResampleImageType;
  typedef typename itk::ResampleImageFilter<LabelImageType,LabelImageType> ResampleLabelType;
  typename ResampleImageType::Pointer resampler = ResampleImageType::New();
  typename ResampleLabelType::Pointer labelResampler = ResampleLabelType::New();
  typedef typename itk::NearestNeighborInterpolateImageFunction<ImageType,double> InterpolatorType;
  typedef typename itk::NearestNeighborInterpolateImageFunction<LabelImageType,double> LabelInterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename LabelInterpolatorType::Pointer labelInterpolator = LabelInterpolatorType::New();

  typename ImageType::SizeType size;
  size[0] = std::ceil((xmax-xmin+1+2*xpad)*voxelSize[0]); 
  size[1] = std::ceil((ymax-ymin+1+2*ypad)*voxelSize[1]); 
  size[2] = std::ceil((zmax-zmin+1+2*zpad)*voxelSize[2]);

  lowerIndex[0]=std::floor((xmin-xpad)*voxelSize[0]); 
  lowerIndex[1]=std::floor((ymin-ypad)*voxelSize[1]); 
  lowerIndex[2]=std::floor((zmin-zpad)*voxelSize[2]);
  
  resampler->SetInterpolator(interpolator);
  typename ImageType::SpacingType spacing;
  spacing.Fill(1);
  resampler->SetOutputSpacing(spacing);
  resampler->SetOutputOrigin(inputImage->GetOrigin());
  resampler->SetOutputDirection(inputImage->GetDirection());
  resampler->SetOutputStartIndex(lowerIndex);
  resampler->SetSize(size);
  resampler->SetInput(inputImage);
  resampler->Update();
  
  labelResampler->SetInterpolator(labelInterpolator);
  labelResampler->SetOutputSpacing(spacing);
  labelResampler->SetOutputOrigin(inputLabel->GetOrigin());
  labelResampler->SetOutputDirection(inputLabel->GetDirection());
  labelResampler->SetOutputStartIndex(lowerIndex);
  labelResampler->SetSize(size);
  labelResampler->SetInput(inputLabel);
  labelResampler->Update();
  
  // determine highest average, provided center voxel is still within label
  typedef typename itk::BinaryBallStructuringElement<PixelType, 3> StructuringElementType;
  StructuringElementType structuringElement;
  structuringElement.SetRadius(QI_PEAK_RADIUS); // approx 1 cm^3 sphere
  structuringElement.CreateStructuringElement();
  typedef typename itk::ShapedNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  NeighborhoodIteratorType nit(structuringElement.GetRadius(),resampler->GetOutput(),resampler->GetOutput()->GetRequestedRegion());
  nit.CreateActiveListFromNeighborhood(structuringElement);
  nit.GoToBegin();
  LabelIteratorType Lit(labelResampler->GetOutput(),labelResampler->GetOutput()->GetRequestedRegion());
  Lit.GoToBegin();
  double peak = itk::NumericTraits<double>::min();
  typename ImageType::IndexType peakIndex;
  while(!nit.IsAtEnd())
  {
    if(Lit.Get()==m_CurrentLabel)
    {
      typename NeighborhoodIteratorType::Iterator sphereIt = nit.Begin();
      double sphereMean = 0.0;
      double activeIndices = nit.GetActiveIndexListSize();
      while(!sphereIt.IsAtEnd())
      {
        sphereMean += sphereIt.Get();
        ++sphereIt;
      }
      sphereMean /= activeIndices;
      if(sphereMean>peak)
      {
        peak = sphereMean;
        peakIndex = nit.GetIndex();
      }
    }
    ++nit; ++Lit;
  }
  typename ImageType::PointType peakPoint;
  resampler->GetOutput()->TransformIndexToPhysicalPoint(peakIndex, peakPoint);
  
  m_NaivePeak = peak;
  IndexType truePeakIndex;
  this->GetInputImage()->TransformPhysicalPointToIndex(peakPoint, truePeakIndex);
  m_NaivePeakIndex = truePeakIndex;
  m_NaivePeakLocation = peakPoint;
}


//----------------------------------------------------------------------------
template <class TImage, class TLabelImage, class TInterpolator>
void
PeakIntensityFilter<TImage, TLabelImage, TInterpolator>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // namespace

#endif
