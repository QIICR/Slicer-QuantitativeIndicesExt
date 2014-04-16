#include "QuantitativeIndicesCLICLP.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <iostream>

#include "itkQuantitativeIndicesComputationFilter.h"

//versioning info
#include "vtkQuantitativeIndicesExtVersionConfigure.h"

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
  if(!First_Quartile){writeFile << "First_Quartile_s = --" << endl;};
  if(!Median){writeFile << "Median_s = --" << endl;};
  if(!Third_Quartile){writeFile << "Third_Quartile_s = --" << endl;};
  if(!Upper_Adjacent){writeFile << "Upper_Adjacent_s = --" << endl;};
  if(!Eightieth){writeFile << "Eightieth_s = --" << endl;};
  if(!NinetyFifth){writeFile << "NinetyFifth_s = --" << endl;};
  if(!TLG){writeFile << "TLG_s = --" << endl;};
  if(!Glycolysis_Q1){writeFile << "Glycolysis_Q1_s = --" << endl;};
  if(!Glycolysis_Q2){writeFile << "Glycolysis_Q2_s = --" << endl;};
  if(!Glycolysis_Q3){writeFile << "Glycolysis_Q3_s = --" << endl;};
  if(!Glycolysis_Q4){writeFile << "Glycolysis_Q4_s = --" << endl;};
  if(!Q1_Distribution){writeFile << "Q1_Distribution_s = --" << endl;};
  if(!Q2_Distribution){writeFile << "Q2_Distribution_s = --" << endl;};
  if(!Q3_Distribution){writeFile << "Q3_Distribution_s = --" << endl;};
  if(!Q4_Distribution){writeFile << "Q4_Distribution_s = --" << endl;};
  if(!SAM){writeFile << "SAM_s = --" << endl;};
  if(!SAM_Background){writeFile << "SAM_Background_s = --" << endl;};
  if(!Peak){writeFile << "Peak_s = --" << endl;};

  typedef itk::QuantitativeIndicesComputationFilter<ImageType,LabelImageType> QIFilterType;
  QIFilterType::Pointer qiCompute = QIFilterType::New();
  qiCompute->SetInputImage(ptImage->GetOutput());
  qiCompute->SetInputLabelImage(labelImage->GetOutput());
  qiCompute->SetCurrentLabel( (int)Label_Value );
  qiCompute->Update();


  if(Mean||RMS||Variance||Max||Min||Volume||TLG||Glycolysis_Q1||Glycolysis_Q2||Glycolysis_Q3||Glycolysis_Q4||Q1_Distribution||Q2_Distribution||Q3_Distribution||Q4_Distribution)
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
      if(TLG){
        writeFile << "TLG_s = " << (double) qiCompute->GetTotalLesionGlycolysis() << endl;
        cout << "TLG: " << (double) qiCompute->GetTotalLesionGlycolysis() << endl;
      }
      if(Glycolysis_Q1){
        writeFile << "Glycolysis_Q1_s = " << (double) qiCompute->GetGly1() << endl;
        cout << "Glycolysis Q1: " << (double) qiCompute->GetGly1() << endl;
      }
      if(Glycolysis_Q2){
        writeFile << "Glycolysis_Q2_s = " << (double) qiCompute->GetGly2() << endl;
        cout << "Glycolysis Q2: " << (double) qiCompute->GetGly2() << endl;
      }
      if(Glycolysis_Q3){
        writeFile << "Glycolysis_Q3_s = " << (double) qiCompute->GetGly3() << endl;
        cout << "Glycolysis Q3: " << (double) qiCompute->GetGly3() << endl;
      }
      if(Glycolysis_Q4){
        writeFile << "Glycolysis_Q4_s = " << (double) qiCompute->GetGly4() << endl;
        cout << "Glycolysis Q4: " << (double) qiCompute->GetGly4() << endl;
      }
      if(Q1_Distribution){
        writeFile << "Q1_Distribution_s = " << (double) qiCompute->GetQ1() << endl;
        cout << "Q1 Distribution: " << (double) qiCompute->GetQ1() << endl;
      }
      if(Q2_Distribution){
        writeFile << "Q2_Distribution_s = " << (double) qiCompute->GetQ2() << endl;
        cout << "Q2 Distribution: " << (double) qiCompute->GetQ2() << endl;
      }
      if(Q3_Distribution){
        writeFile << "Q3_Distribution_s = " << (double) qiCompute->GetQ3() << endl;
        cout << "Q3 Distribution: " << (double) qiCompute->GetQ3() << endl;
      }
      if(Q4_Distribution){
        writeFile << "Q4_Distribution_s = " << (double) qiCompute->GetQ4() << endl;
        cout << "Q4 Distribution: " << (double) qiCompute->GetQ4() << endl;
      }
    }

  if(First_Quartile || Median || Third_Quartile || Upper_Adjacent || Eightieth || NinetyFifth)
    {
      qiCompute->CalculateQuartiles();
      if(First_Quartile){
        writeFile << "First_Quartile_s = " << (double) qiCompute->GetFirstQuartileValue() << endl;
        cout << "1st Quartile: " << (double) qiCompute->GetFirstQuartileValue() << endl;
      }
      if(Median){
        writeFile << "Median_s = " << (double) qiCompute->GetMedianValue() << endl;
        cout << "Median: " << (double) qiCompute->GetMedianValue() << endl;
      }
      if(Third_Quartile){
        writeFile << "Third_Quartile_s = " << (double) qiCompute->GetThirdQuartileValue() << endl;
        cout << "3rd Quartile: " << (double) qiCompute->GetThirdQuartileValue() << endl;
      }
      if(Upper_Adjacent){
        writeFile << "Upper_Adjacent_s = " << (double) qiCompute->GetUpperAdjacentValue() << endl;
        cout << "Upper Adjacent: " << (double) qiCompute->GetUpperAdjacentValue() << endl;
      }
      if(Eightieth){
        writeFile << "Eightieth_s = " << (double) qiCompute->GetPercentile80Value() << endl;
        cout << "80th Percentile: " << (double) qiCompute->GetPercentile80Value() << endl;
      }
      if(NinetyFifth){
        writeFile << "NinetyFifth_s = " << (double) qiCompute->GetPercentile95Value() << endl;
        cout << "95th Percentile: " << (double) qiCompute->GetPercentile95Value() << endl;
      }
    }

  if(SAM||SAM_Background)
    {
      qiCompute->CalculateSAM();
      if(SAM){
        writeFile << "SAM_s = " << (double) qiCompute->GetSAMValue() << endl;
        cout << "SAM: " << (double) qiCompute->GetSAMValue() << endl;
      }
      if(SAM_Background){
        writeFile << "SAM_Background_s = " << (double) qiCompute->GetSAMBackground() << endl;
        cout << "SAM mean background: " << (double) qiCompute->GetSAMBackground() << endl;
      }
    }

  if(Peak)
    {
      qiCompute->CalculatePeak();
      writeFile << "Peak_s = " << (double) qiCompute->GetPeakValue() << endl;
      cout << "Peak: " << (double) qiCompute->GetPeakValue() << endl;
    }
    
  writeFile << "Software_Version = " << QuantitativeIndicesExt_WC_REVISION << endl;

  writeFile.close();

  return EXIT_SUCCESS;
}
