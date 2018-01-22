void FitProfile(TH2F *h);

void PlotProfileLKL()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  using namespace RooFit;

  TChain *t = new TChain("t");
  t->Add("../../analysis/ProfileLikelihood/bscale_bb2n/FitResult_*.root");

  TH2F *h = new TH2F("h","h",20,4200,5000,11,0.94,1.06);

  RooFitResult *result = 0;
  double bScale;
  double numbb2n;

  t->SetBranchAddress("result",&result);
  t->SetBranchAddress("bScale",&bScale);
  t->SetBranchAddress("numbb2n",&numbb2n);

  int nentries = t->GetEntries();
  for (int i = 0; i < nentries; i++) {
    t->GetEntry(i);

    double m = result->minNll();
    int bin = h->FindBin(numbb2n,bScale);
    if (m > 0) {m = -2639;}
    h->SetBinContent(bin,m);
  }

  h->SetMaximum(2630);
  
  double ymin = 0.94;
  double ymax = 1.06;
  double xmin = 4200;
  double xmax = 5000;
  double xstp = 8;
  double ystp = 0.0012;
  int count = 0;
  
  TGraph2D *grTMP = new TGraph2D(h);
  
  TGraphDelaunay *gr_Delaunay = new TGraphDelaunay(grTMP);
  
  TGraph2D *gr = new TGraph2D(10000);
  for (double j = xmin; j <= xmax; j+=xstp) {
    for (double k = ymin; k <= ymax; k+=ystp) {
      gr->SetPoint(count,j,k,gr_Delaunay->ComputeZ(j,k));
      cout << gr_Delaunay->Interpolate(j,k) << "\t" << gr_Delaunay->ComputeZ(j,k) << endl;
      //gr->SetPointError(count,xstp,ystp,h->GetBinError(h->FindBin(j,k)));
      
      count++;
    }
  }
  
  gr->SetMaximum(-2630);
  
  //gr->Draw("colz");
  
  //TGraph2D *gr2D = new TGraph2D(h);
  
  TF2 *chifit = new TF2("chifit","[5]+(([0]-x)*(-([0]-x)/[2]/[2]/(-1.+[4]*[4])+[4]*([1]-y)/[2]/[3]/(-1.+[4]*[4]))+([1]-y)*([4]*([0]-x)/[2]/[3]/(-1.+[4]*[4])-([1]-y)/[3]/[3]/(-1.+[4]*[4])))",4200,5000,0.96,1.04);
  
  double xmin = 4220;
  double ymin = 0.96;
  double xmax = 4980;
  double ymax = 1.05;
  
  double xbest = 4550;
  double ybest = 1.01;
  double xerr = 200;
  double yerr = 0.005;
  double corcoef = 0.1;
  double minval = -2639;
  
  chifit->SetRange(xmin,ymin,xmax,ymax);
  chifit->SetParameters(xbest,ybest,xerr,yerr,corcoef,minval);
  chifit->SetParLimits(0,xmin,xmax);
  chifit->SetParLimits(1,ymin,ymax);
  chifit->SetParLimits(2,xmin*0.001,xmin);
  chifit->SetParLimits(3,ymin*0.001,ymin);
  chifit->SetParLimits(4,-0.5,0.5);
  chifit->SetParLimits(5,-2640,-2630);
  
  gr->Fit("chifit","r");
  
  double NLL_min = chifit->GetParameter(5);
  double cont_levels[4] = {NLL_min + 2.41, NLL_min + 4.60, NLL_min + 5.99, NLL_min + 9.21};
  
  chifit->SetContour(4,cont_levels);
  
  chifit->SetLineWidth(2);
  chifit->SetLineStyle(2);
  
  TCanvas *c1 = new TCanvas();
  gr->Draw("colz");
  chifit->Draw("cont3 same");
  
  TCanvas *c2 = new TCanvas();
  chifit->Draw();
  
  //FitProfile(h);

  //h->Draw("cont4");

  return;
}

void FitProfile(TH1F *h)
{
  gROOT->SetStyle("Plain");
  
  TF2 *chifit = new TF2("chifit","[5]+(([0]-x)*(-([0]-x)/[2]/[2]/(-1.+[4]*[4])+[4]*([1]-y)/[2]/[3]/(-1.+[4]*[4]))+([1]-y)*([4]*([0]-x)/[2]/[3]/(-1.+[4]*[4])-([1]-y)/[3]/[3]/(-1.+[4]*[4])))",4200,5000,0.94,1.06);
  
  double xmin = 4200;
  double ymin = 0.94;
  double xmax = 5000;
  double ymax = 1.06;
  
  double xbest = 4400;
  double ybest = 1.01;
  double xerr = 200;
  double yerr = 0.01;
  double corcoef = -0.2;
  double minval = 2630;
  
  chifit->SetRange(xmin,ymin,xmax,ymax);
  chifit->SetParameters(xbest,ybest,xerr,yerr,corcoef,minval);
  chifit->SetParLimits(0,xmin,xmax);
  chifit->SetParLimits(1,ymin,ymax);
  chifit->SetParLimits(2,xmin*0.001,xmin);
  chifit->SetParLimits(3,ymin*0.001,ymin);
  chifit->SetParLimits(4,-1,1);
  chifit->SetParLimits(5,2500,2700);
  
  h->Fit("chifit","r");
  
  h->Draw("cont4");
  chifit->Draw("same");
  
  //chifit->Draw();
  
  return;
}

void Test()
{
  TF2 *f2 = new TF2("f2","sin(x)*sin(y)/(x*y)",0,5,0,5);
  f2->Draw();
}