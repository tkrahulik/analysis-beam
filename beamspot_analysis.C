/*---------------Beam Spot Analysis-----------------
 * Project: A Compact Magnetic Field Cloaking Device
 * Author: Thomas Krahulik
 * Date: June 2, 2016
 * Purpose: To develop a program that locates the
 * center of a beam spot and performs an analysis
 * of any shift of the beam.
 * To compile and execute run the command:
 * root -L beamspot_analysis.C++
 ---------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TROOT.h>
#include <TFile.h>
#include <TF1.h>
#include <TImage.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TStyle.h>

using namespace std;

void beamspot_analysis(
		       TString in_image = "../Images/Test_2016_07_08/IMG_0017_edit.png"
		       )
{
  //Read in Configuration File of Conversions and Measurements
  const char* fconfig = "beamspot_config.txt";
  cout<< "processing file " << fconfig <<endl;
  std::ifstream infile(fconfig);

  std::string label, FullLine;
  double convert, pipe;

  getline(infile, FullLine);

  while ( !infile.eof() )
    {
      if ( FullLine[0]=='#' || FullLine.substr(0, 2) == "//" )
	{
	  getline(infile,FullLine);
	  continue;
	}

      istringstream line( FullLine.c_str() );
      line >> label;
      if (label == "conversion"){
	line >> convert;
	std::cout << "Conversion: " << convert << endl;
      }
      else if (label == "pipe_distance"){
	line >> pipe;
	std::cout << "Beam Pipe Length: " << pipe << endl;
      }

      getline(infile, FullLine);

    }

  gStyle->SetOptStat(kFALSE);
  TImage *image = TImage::Open(in_image);
  UInt_t xPixels = image->GetWidth();
  UInt_t yPixels = image->GetHeight();
  UInt_t *argb   = image->GetArgbArray();

  std::cout << "Pixels along x-axis: " << xPixels << endl;
  std::cout << "Pixels along y-axis: " << yPixels << endl;

  /*
   *----------Intensity Histogram----------
   * Fill intensity histogram to determine 
   * where to make cuts to exclude noise
   */
  TH1D* h_I = new TH1D("h_I", "Intensity Distribution", 100, 0, 1);

  for (UInt_t row=0; row<xPixels; ++row) {
    for (UInt_t col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      float gray = float(argb[index]&0xff)/256;
      h_I->Fill(gray);
    }
  }

  TCanvas *c_h_I = new TCanvas();
  c_h_I->Draw();
  h_I->Draw();

  float I_lowcut = 0.8;
  float I_highcut = 1.0;

  /*
   *---------Beam Spot Histograms-----------
   * Turn beam spot images into two-dimensional histograms with intensity on z-axis
   * One histogram is entire image without intensity cut
   * Second histogram makes intensity cuts
   * and removes lower right corner to exclude time stamp
   */

  TH2D* h_rep = new TH2D("h_rep","Beam Spot",xPixels,-1,1,yPixels,-1,1);
  TH2D* h_cut = new TH2D("h_cut","Beam Spot with Intensity Cut",xPixels,-1,1,yPixels,-1,1);

  float weight_sumx = 0;
  float weight_sumy = 0;
  float gray_sum = 0;

  float SQweight_sumx = 0;
  float SQweight_sumy = 0;

  for (UInt_t row=0; row<xPixels; ++row) {
    for (UInt_t col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      float gray = float(argb[index]&0xff)/256.;

      h_rep->SetBinContent(row+1,yPixels-col,gray);
      
      //Intensity Cut
      if ((gray >= I_lowcut) & (gray <= I_highcut)){
	//Cut out time stamp in lower right corner
	if ((row < yPixels-100.) & (col < xPixels-100.)){

	  h_cut->SetBinContent(row+1,yPixels-col,gray);

	  weight_sumx = weight_sumx + (row+1)*gray;
	  weight_sumy = weight_sumy + (col)*gray;
	  gray_sum = gray_sum + gray;
	  
	  //Calculate for Variance of Center Coordinates
	  SQweight_sumx = SQweight_sumx + (pow((row+1),2))*gray;
	  SQweight_sumy = SQweight_sumy + (pow(col,2))*gray;

	}      
      }
    }
  }
  //Center of Intensity
  float pix_x = weight_sumx / gray_sum;
  float pix_y = weight_sumy / gray_sum;
  //Mean of the Squares
  float SQpix_x = SQweight_sumx / gray_sum;
  float SQpix_y = SQweight_sumy / gray_sum;
  //Standard Deviation
  float stdx = sqrt(SQpix_x - pow(pix_x,2));
  float stdy = sqrt(SQpix_y - pow(pix_y,2));
  //Convert Pixels to Physical Distances
  float centerx = pix_x * convert;
  float centery = pix_y * convert;

  std::cout << "Mean x: " << centerx << endl;
  std::cout << "std x: " << stdx << endl;
  std::cout << "Mean y: " << centery << endl;
  std::cout << "std y: " << stdy << endl;

  gStyle->SetPalette(55);
  gStyle->SetCanvasDefH(yPixels);
  gStyle->SetCanvasDefW(xPixels);
  TCanvas *c_rep = new TCanvas();
  c_rep->Draw();
  h_rep->Draw("colz");
  TCanvas *c_cut = new TCanvas();
  c_cut->Draw();
  h_cut->Draw("colz");

  TCanvas *c_image = new TCanvas();
  c_image->Draw();
  image->Draw();
  image->DrawEllips(pix_x, pix_y, stdx, stdy, 1, "#FF0000", 3);
}
