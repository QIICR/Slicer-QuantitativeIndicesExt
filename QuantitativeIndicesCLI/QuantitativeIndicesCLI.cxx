#include "QuantitativeIndicesCLICLP.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <iostream>

#include "itkQuantitativeIndicesComputationFilter.h"

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

  ptImage->SetFileName( Grayscale_Image );
  labelImage->SetFileName( Label_Image );
  ptImage->Update();
  labelImage->Update();

  ofstream writeFile;
  writeFile.open( returnParameterFile.c_str() );
  if(!Mean){writeFile << "Mean_s = --" << endl;};
  if(!Variance){writeFile << "Variance_s = --" << endl;};
  if(!RMS){writeFile << "RMS_s = --" << endl;};
  if(!Max){writeFile << "Max_s = --" << endl;};
  if(!Min){writeFile << "Min_s = --" << endl;};
  if(!Volume){writeFile << "Volume_s = --" << endl;};
  if(!Quart1){writeFile << "Quart1_s = --" << endl;};
  if(!Median){writeFile << "Median_s = --" << endl;};
  if(!Quart3){writeFile << "Quart3_s = --" << endl;};
  if(!Adj){writeFile << "Adj_s = --" << endl;};
  if(!MTV){writeFile << "MTV_s = --" << endl;};
  if(!Gly1){writeFile << "Gly1_s = --" << endl;};
  if(!Gly2){writeFile << "Gly2_s = --" << endl;};
  if(!Gly3){writeFile << "Gly3_s = --" << endl;};
  if(!Gly4){writeFile << "Gly4_s = --" << endl;};
  if(!Q1){writeFile << "Q1_s = --" << endl;};
  if(!Q2){writeFile << "Q2_s = --" << endl;};
  if(!Q3){writeFile << "Q3_s = --" << endl;};
  if(!Q4){writeFile << "Q4_s = --" << endl;};
  if(!SAM){writeFile << "SAM_s = --" << endl;};
  if(!SAMBG){writeFile << "SAMBG_s = --" << endl;};
  if(!Peak){writeFile << "Peak_s = --" << endl;};

  typedef itk::QuantitativeIndicesComputationFilter<ImageType,LabelImageType> QIFilterType;
  QIFilterType::Pointer qiCompute = QIFilterType::New();
  qiCompute->SetInputImage(ptImage->GetOutput());
  qiCompute->SetInputLabelImage(labelImage->GetOutput());
  qiCompute->SetCurrentLabel( (int)Label_Value );
  qiCompute->Update();


  if(Mean||RMS||Variance||Max||Min||Volume||MTV||Gly1||Gly2||Gly3||Gly4||Q1||Q2||Q3||Q4)
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
      if(Variance){
        double var = (double) qiCompute->GetStandardDeviation();
        writeFile << "Variance_s = " << var*var << endl;
        cout << "Variance: " << var*var << endl;
      }
      if(Max){
        writeFile << "Max_s = " << (double) qiCompute->GetMaximumValue() << endl;
        cout << "Max: " << (double) qiCompute->GetMaximumValue() << endl;
      }
      if(Min){
        writeFile << "Min_s = " << (double) qiCompute->GetMinimumValue() << endl;
        cout << "Min: " << (double) qiCompute->GetMinimumValue() << endl;
      }
      if(Volume){
        writeFile << "Volume_s = " << (double) qiCompute->GetSegmentedVolume() << endl;
        cout << "Volume: " << (double) qiCompute->GetSegmentedVolume() << endl;
      }
      if(MTV){
        writeFile << "MTV_s = " << (double) qiCompute->GetMetabolicVolume() << endl;
        cout << "MTV: " << (double) qiCompute->GetMetabolicVolume() << endl;
      }
      if(Gly1){
        writeFile << "Gly1_s = " << (double) qiCompute->GetGly1() << endl;
        cout << "Gly1: " << (double) qiCompute->GetGly1() << endl;
      }
      if(Gly2){
        writeFile << "Gly2_s = " << (double) qiCompute->GetGly2() << endl;
        cout << "Gly2: " << (double) qiCompute->GetGly2() << endl;
      }
      if(Gly3){
        writeFile << "Gly3_s = " << (double) qiCompute->GetGly3() << endl;
        cout << "Gly3: " << (double) qiCompute->GetGly3() << endl;
      }
      if(Gly4){
        writeFile << "Gly4_s = " << (double) qiCompute->GetGly4() << endl;
        cout << "Gly4: " << (double) qiCompute->GetGly4() << endl;
      }
      if(Q1){
        writeFile << "Q1_s = " << (double) qiCompute->GetQ1() << endl;
        cout << "Q1: " << (double) qiCompute->GetQ1() << endl;
      }
      if(Q2){
        writeFile << "Q2_s = " << (double) qiCompute->GetQ2() << endl;
        cout << "Q2: " << (double) qiCompute->GetQ2() << endl;
      }
      if(Q3){
        writeFile << "Q3_s = " << (double) qiCompute->GetQ3() << endl;
        cout << "Q3: " << (double) qiCompute->GetQ3() << endl;
      }
      if(Q4){
        writeFile << "Q4_s = " << (double) qiCompute->GetQ4() << endl;
        cout << "Q4: " << (double) qiCompute->GetQ4() << endl;
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

  if(SAM||SAMBG)
    {
      qiCompute->CalculateSAM();
      if(SAM){
        writeFile << "SAM_s = " << (double) qiCompute->GetSAMValue() << endl;
        cout << "SAM: " << (double) qiCompute->GetSAMValue() << endl;
      }
      if(SAMBG){
        writeFile << "SAMBG_s = " << (double) qiCompute->GetSAMBackground() << endl;
        cout << "SAM mean background: " << (double) qiCompute->GetSAMBackground() << endl;
      }
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
