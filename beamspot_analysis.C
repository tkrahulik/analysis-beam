/*---------------Beam Spot Analysis-----------------
 * Project: A Compact Magnetic Field Cloaking Device
 * Author: Thomas Krahulik
 * Date: June 2, 2016
 * Purpose: To develop a program that locates the
 * center of a beam spot and performs an analysis
 * of any shift of the beam.
 ---------------------------------------------------*/
void beamspot_analysis(
		       const TString in_image = "../Images/beam_2014_10_30_001.jpg"
		       )
{
  eicStyle->SetOptStat(kFALSE);
  TImage *image = TImage::Open(in_image);
  UInt_t xPixels = image->GetWidth();
  UInt_t yPixels = image->GetHeight();
  UInt_t *argb   = image->GetArgbArray();

  std::cout << xPixels << endl;
  std::cout << yPixels << endl;

  /*
   *----------Intensity Histogram----------
   * Fill intensity histogram to determine 
   * where to make cuts to exclude noise
   */
  TH1D* h_I = new TH1D("h_I", "Intensity Distribution", 100, 0, 1);

  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      float gray = float(argb[index]&0xff)/256;
      h_I->Fill(gray);
    }
  }

  TCanvas *c_h_I = new TCanvas();
  h_I->Draw();
  h_I->Fit("gaus", "", "", 0.7, 1);

  /*
   * Calculate cuts on intensity using parameters from Gaussian fit
   */
  float mean = h_I->GetFunction("gaus")->GetParameter(1);
  float sigma = h_I->GetFunction("gaus")->GetParameter(2);

  float I_lowcut = mean - 5*sigma;
  float I_highcut = mean + 5*sigma;

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

  for (int row=0; row<xPixels; ++row) {
    for (int col=0; col<yPixels; ++col) {
      int index = col*xPixels+row;
      float gray = float(argb[index]&0xff)/256;

      h_rep->SetBinContent(row+1,yPixels-col,gray);
      
      //Intensity Cut
      if (gray >= I_lowcut & gray <= I_highcut){
	//Cut out time stamp in lower right corner
	if (row < 700 & col < 1200){

	  h_cut->SetBinContent(row+1,yPixels-col,gray);

	  weight_sumx = weight_sumx + (row+1)*gray;
	  weight_sumy = weight_sumy + (col)*gray;
	  gray_sum = gray_sum + gray;
	  
	  //Calculate for Variance of Center Coordinates
	  SQweight_sumx = SQweight_sumx + ((row+1)**2)*gray;
	  SQweight_sumy = SQweight_sumy + ((col)**2)*gray;

	}      
      }
    }
  }
  //Center of Intensity
  float centerx = weight_sumx / gray_sum;
  float centery = weight_sumy / gray_sum;
  //Mean of the Squares
  float SQcenterx = SQweight_sumx / gray_sum;
  float SQcentery = SQweight_sumy / gray_sum;
  //Standard Deviation
  float stdx = sqrt(SQcenterx - centerx**2);
  float stdy = sqrt(SQcentery - centery**2);

  std::cout << "mean x: " << centerx << endl;
  std::cout << "std x: " << stdx << endl;
  std::cout << "mean y: " << centery << endl;
  std::cout << "std y: " << stdy << endl;

  eicStyle->SetPalette(55);//grayscale
  eicStyle->SetCanvasDefH(yPixels);
  eicStyle->SetCanvasDefW(xPixels);
  TCanvas *c_rep = new TCanvas();
  h_rep->Draw("colz");
  TCanvas *c_cut = new TCanvas();
  h_cut->Draw("colz");

  TCanvas *c_image = new TCanvas();
  image->Draw();
  image->DrawEllips(centerx, centery, stdx, stdy, 1, "#FF0000", 3);
}
