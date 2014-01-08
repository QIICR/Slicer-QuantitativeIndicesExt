#include "QuantitativeIndicesCLP.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <iostream>

#include "itkQIsComputationFilter.h"

using namespace std;

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef  float  PixelType;
	const unsigned int Dimension = 3;

  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::Image< int, Dimension > LabelImageType;

	//image reader
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileReader< LabelImageType > LabelReaderType;
	ReaderType::Pointer ptImage = ReaderType::New();
  LabelReaderType::Pointer labelImage = LabelReaderType::New();

  ptImage->SetFileName( PET_Image );
  labelImage->SetFileName( Label_Image );
  ptImage->Update();
  labelImage->Update();

  ofstream writeFile;
  writeFile.open( returnParameterFile.c_str() );
  if(!Mean){writeFile << "Mean_s = --" << endl;};
  if(!RMS){writeFile << "RMS_s = --" << endl;};
  if(!Quart1){writeFile << "Quart1_s = --" << endl;};
  if(!Median){writeFile << "Median_s = --" << endl;};
  if(!Quart3){writeFile << "Quart3_s = --" << endl;};
  if(!Adj){writeFile << "Adj_s = --" << endl;};
  if(!SAM){writeFile << "SAM_s = --" << endl;};
  if(!Peak){writeFile << "Peak_s = --" << endl;};

  typedef itk::QIsComputationFilter<ImageType,LabelImageType> QIFilterType;
  QIFilterType::Pointer qiCompute = QIFilterType::New();
  qiCompute->SetInputImage(ptImage->GetOutput());
  qiCompute->SetInputLabelImage(labelImage->GetOutput());
  qiCompute->SetCurrentLabel( (int)Label_Value );
  qiCompute->Update();


  if(Mean || RMS)
    {
      qiCompute->CalculateMean();
      if(Mean){
        writeFile << "Mean_s = " << (double) qiCompute->GetAverageValue() << endl;
        cout << "Mean: " << (double) qiCompute->GetAverageValue() << endl;
      }
      if(RMS){
        writeFile << "RMS_s = " << (double) qiCompute->GetRMSValue() << endl;
        cout << "RMS: " << (double) qiCompute->GetRMSValue() << endl;
      }
    }

  if(Quart1 || Median || Quart3 || Adj)
    {
      qiCompute->CalculateQuartiles();
      if(Quart1){
        writeFile << "Quart1_s = " << (double) qiCompute->GetFirstQuartileValue() << endl;
        cout << "1st Quartile: " << (double) qiCompute->GetFirstQuartileValue() << endl;
      }
      if(Median){
        writeFile << "Median_s = " << (double) qiCompute->GetMedianValue() << endl;
        cout << "Median: " << (double) qiCompute->GetMedianValue() << endl;
      }
      if(Quart3){
        writeFile << "Quart3_s = " << (double) qiCompute->GetThirdQuartileValue() << endl;
        cout << "3rd Quartile: " << (double) qiCompute->GetThirdQuartileValue() << endl;
      }
      if(Adj){
        writeFile << "Adj_s = " << (double) qiCompute->GetUpperAdjacentValue() << endl;
        cout << "Upper Adjacent: " << (double) qiCompute->GetUpperAdjacentValue() << endl;
      }
    }

  if(SAM)
    {
      qiCompute->CalculateSAM();
      writeFile << "SAM_s = " << (double) qiCompute->GetSAMValue() << endl;
      cout << "SAM: " << (double) qiCompute->GetSAMValue() << endl;
    }

  if(Peak)
    {
      qiCompute->CalculatePeak();
      writeFile << "Peak_s = " << (double) qiCompute->GetPeakValue() << endl;
      cout << "Peak: " << (double) qiCompute->GetPeakValue() << endl;
    }

  writeFile.close();

  return EXIT_SUCCESS;
}
