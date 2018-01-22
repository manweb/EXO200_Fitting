double fitFunctionTh228(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  
  return gauss1 + erf1;
}

void EvaluateFitFunctionTh228()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  using namespace RooFit;
  
  gSystem->Load("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/lib/libEXOFitting.so");
  
  EXOFitter *fitter = new EXOFitter(true);
  fitter->LoadWorkspace("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/analysis/EXO_Workspace_floatEres.root");
  //fitter->BuildRotationModel();
  
  RooWorkspace *w = fitter->GetWorkspace();
  
  RooRealVar *energy_ss = (RooRealVar*)w->var("energy_ss");
  RooRealVar *energy_ms = (RooRealVar*)w->var("energy_ms");
  
  // open file with MC histograms
  TFile *f = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/MC_Source_Histos.root","READ");
  
  // get the histograms
  TH1F *hTh228_ss = (TH1F*)f->Get("SourceP4_px_228Th_ss");
  TH1F *hTh228_ms = (TH1F*)f->Get("SourceP4_px_228Th_ms");
  
  // generate RooDataHist
  RooDataHist hTh228_d_ss("hTh228_d_ss","hTh228_d_ss",*energy_ss,hTh228_ss);
  RooDataHist hTh228_d_ms("hTh228_d_ms","hTh228_d_ms",*energy_ms,hTh228_ms);
  
  // generate RooHistPdf
  RooHistPdf Th228_Pdf_ss("Th228_Pdf_ss","Th228_Pdf_ss",*energy_ss,hTh228_d_ss,0);
  RooHistPdf Th228_Pdf_ms("Th228_Pdf_ms","Th228_Pdf_ms",*energy_ms,hTh228_d_ms,0);
  
  //RooHistPdf *Th228_Pdf_ss = (RooHistPdf*)w->pdf("SourceP4_px_228Th_Pdf_ss");
  //RooHistPdf *Th228_Pdf_ms = (RooHistPdf*)w->pdf("SourceP4_px_228Th_Pdf_ms");
  
  // smear PDF
  RooAbsPdf *Th228_Pdf_s_ss = fitter->SmearPdf(&Th228_Pdf_ss,energy_ss);
  RooAbsPdf *Th228_Pdf_s_ms = fitter->SmearPdf(&Th228_Pdf_ms,energy_ms);
  
  // initialize arrays
  const int n = 10;
  double *resolution_ss = new double[n];
  double *resolution_ms = new double[n];
  double *Th228_Peak_ss = new double[n];
  double *Th228_Peak_ms = new double[n];
  
  double *Th228_Peak_err_ss = new double[n];
  double *Th228_Peak_err_ms = new double[n];
  
  double *sigma_ss = new double[n];
  double *sigma_ms = new double[n];
  double *Th228_Sigma_ss = new double[n];
  double *Th228_Sigma_ms = new double[n];
  
  double *Th228_Sigma_err_ss = new double[n];
  double *Th228_Sigma_err_ms = new double[n];
  
  // generat fit functions
  TF1 *Th228_fit_ss = new TF1("Th228_fit_ss",fitFunctionTh228,2400,3000,4);
  TF1 *Th228_fit_ms = new TF1("Th228_fit_ms",fitFunctionTh228,2400,2800,4);
  
  TF1 *erf_ss = new TF1("erf_ss","[3]*[0]*0.5*TMath::Erfc((x - [1]) / (TMath::Sqrt(2) * [2]))",2400,3000);
  TF1 *gaus_ss = new TF1("gaus_ss","gaus",2400,3000);
  
  TF1 *erf_ms = new TF1("erf_ms","[3]*[0]*0.5*TMath::Erfc((x - [1]) / (TMath::Sqrt(2) * [2]))",2400,2800);
  TF1 *gaus_ms = new TF1("gaus_ms","gaus",2400,2800);
  
  Th228_fit_ss->SetLineColor(kBlue);
  Th228_fit_ms->SetLineColor(kBlue);
  
  Th228_fit_ss->SetLineWidth(2);
  Th228_fit_ms->SetLineWidth(2);
  
  erf_ss->SetLineColor(kBlue);
  gaus_ss->SetLineColor(kBlue);
  
  erf_ss->SetLineWidth(1);
  gaus_ss->SetLineWidth(1);
  
  erf_ss->SetLineStyle(2);
  gaus_ss->SetLineStyle(2);
  
  erf_ms->SetLineColor(kBlue);
  gaus_ms->SetLineColor(kBlue);
  
  erf_ms->SetLineWidth(1);
  gaus_ms->SetLineWidth(1);
  
  erf_ms->SetLineStyle(2);
  gaus_ms->SetLineStyle(2);
  
  Th228_fit_ss->SetParameters(1000,2615,40,0.2);
  Th228_fit_ms->SetParameters(1000,2615,40,0.2);
  
  Th228_fit_ss->SetParLimits(0,0.0,1e6);
  Th228_fit_ss->SetParLimits(1,2500,2700);
  Th228_fit_ss->SetParLimits(2,0.0,200.0);
  Th228_fit_ss->SetParLimits(3,0.01,0.5);
  
  Th228_fit_ms->SetParLimits(0,0.0,1e6);
  Th228_fit_ms->SetParLimits(1,2500,2700);
  Th228_fit_ms->SetParLimits(2,0.0,200.0);
  Th228_fit_ms->SetParLimits(3,0.01,0.5);
  
  RooRealVar res0_ss("r0_ss","r0_ss",0.0);
  RooRealVar res1_ss("r1_ss","r1_ss",0.0);
  RooRealVar res2_ss("r2_ss","r2_ss",0.0);
  
  RooRealVar res0_ms("r0_ms","r0_ms",0.0);
  RooRealVar res1_ms("r1_ms","r1_ms",0.0);
  RooRealVar res2_ms("r2_ms","r2_ms",0.0);
  
  RooArgSet* comp_ss = Th228_Pdf_s_ss->getComponents();
  EXOresolution *res_ss = comp_ss->find("res_ss");
  ((RooRealVar*)res_ss->getA())->setVal(0.0);
  ((RooRealVar*)res_ss->getB())->setVal(0.0);
  ((RooRealVar*)res_ss->getC())->setVal(0.0);
  
  RooArgSet* comp_ms = Th228_Pdf_s_ms->getComponents();
  EXOresolution *res_ms = comp_ms->find("res_ms");
  ((RooRealVar*)res_ms->getA())->setVal(0.0);
  ((RooRealVar*)res_ms->getB())->setVal(0.0);
  ((RooRealVar*)res_ms->getC())->setVal(0.0);
  
  TH1F *h_ss = new TH1F("h_ss","h_ss",100,2000,3000);
  TH1F *h_ms = new TH1F("h_ms","h_ms",100,2000,3000);
  
  double rStart = 0.012;
  double rEnd = 0.045;
  double stp = (rEnd - rStart) / double(n);
  
  double fitRange1_ss = 2350.0;
  double fitRange1_ms = 2300.0;
  double fitRange2_ss = 2900.0;
  double fitRange2_ms = 2800.0;
  
  RooPlot *frame_ss;
  RooPlot *frame_ms;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,2);
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionTh228.ps[");
  
  for (int i = 0; i < n; i++) {
    double r = i*stp + rStart;
    ((RooRealVar*)res_ss->getA())->setVal(r*TMath::Sqrt(2615.0));
    ((RooRealVar*)res_ms->getA())->setVal(r*TMath::Sqrt(2615.0));
    
    fitRange1_ss = -4040.0*r + 2468.0;
    fitRange1_ms = -6734.0*r + 2521.0;
    
    Th228_fit_ss->SetRange(fitRange1_ss,fitRange2_ss);
    erf_ss->SetRange(fitRange1_ss,fitRange2_ss);
    gaus_ss->SetRange(fitRange1_ss,fitRange2_ss);
    
    Th228_fit_ms->SetRange(fitRange1_ms,fitRange2_ms);
    erf_ms->SetRange(fitRange1_ms,fitRange2_ms);
    gaus_ms->SetRange(fitRange1_ms,fitRange2_ms);
    
    // generate toy MC
    RooDataSet *Th228_ToyMC_ss = Th228_Pdf_s_ss->generate(*energy_ss,100000);
    RooDataSet *Th228_ToyMC_ms = Th228_Pdf_s_ms->generate(*energy_ms,100000);
    
    ((TTree*)Th228_ToyMC_ss->tree())->Draw("energy_ss>>h_ss","","goff");
    h_ss->Fit("Th228_fit_ss","rn");
    
    ((TTree*)Th228_ToyMC_ms->tree())->Draw("energy_ms>>h_ms","","goff");
    h_ms->Fit("Th228_fit_ms","rn");
    
    h_ss->SetTitle(Form("resolution = %.4f",r));
    h_ms->SetTitle(Form("resolution = %.4f",r));
    
    resolution_ss[i] = r;
    resolution_ms[i] = r;
    Th228_Peak_ss[i] = 2615.0 - Th228_fit_ss->GetParameter(1);
    Th228_Peak_ms[i] = 2615.0 - Th228_fit_ms->GetParameter(1);
    
    Th228_Peak_err_ss[i] = Th228_fit_ss->GetParError(1);
    Th228_Peak_err_ms[i] = Th228_fit_ms->GetParError(1);
    
    sigma_ss[i] = r*2615.0;
    sigma_ms[i] = r*2615.0;
    Th228_Sigma_ss[i] = Th228_fit_ss->GetParameter(2);
    Th228_Sigma_ms[i] = Th228_fit_ms->GetParameter(2);
    
    Th228_Sigma_err_ss[i] = Th228_fit_ss->GetParError(2);
    Th228_Sigma_err_ms[i] = Th228_fit_ss->GetParError(2);
    
    frame_ss = energy_ss->frame();
    Th228_Pdf_ss.plotOn(frame_ss,LineStyle(2),LineWidth(1));
    Th228_Pdf_s_ss->plotOn(frame_ss);
    frame_ss->SetTitle(Form("Single Site, resolution = %.4f",r));
    
    frame_ms = energy_ms->frame();
    Th228_Pdf_ms.plotOn(frame_ms,LineStyle(2),LineWidth(1));
    Th228_Pdf_s_ms->plotOn(frame_ms);
    frame_ms->SetTitle(Form("Multi Site, resolution = %.4f",r));
    
    erf_ss->SetParameters(Th228_fit_ss->GetParameters());
    gaus_ss->SetParameters(Th228_fit_ss->GetParameters());
    
    erf_ms->SetParameters(Th228_fit_ms->GetParameters());
    gaus_ms->SetParameters(Th228_fit_ms->GetParameters());
    
    c1->cd(1);
    frame_ss->Draw();
    
    c1->cd(2);
    frame_ms->Draw();
    
    c1->cd(3);
    h_ss->Draw("EZP");
    Th228_fit_ss->Draw("same");
    erf_ss->Draw("same");
    gaus_ss->Draw("same");
    
    c1->cd(4);
    h_ms->Draw("EZP");
    Th228_fit_ms->Draw("same");
    erf_ms->Draw("same");
    gaus_ms->Draw("same");
    
    c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionTh228.ps");
    
    //h_ss->Delete();
    //h_ms->Delete();
  }
  
  TGraphErrors *grDev_ss = new TGraphErrors(n,resolution_ss,Th228_Peak_ss,0,Th228_Peak_err_ss);
  TGraphErrors *grDev_ms = new TGraphErrors(n,resolution_ms,Th228_Peak_ms,0,Th228_Peak_err_ms);
  TGraphErrors *grSigma_ss = new TGraphErrors(n,sigma_ss,Th228_Sigma_ss,0,Th228_Sigma_err_ss);
  TGraphErrors *grSigma_ms = new TGraphErrors(n,sigma_ms,Th228_Sigma_ms,0,Th228_Sigma_err_ms);
  
  grDev_ss->SetTitle("Peak bias vs resolution, single site");
  grDev_ms->SetTitle("Peak bias vs resolution, multi site");
  
  grSigma_ss->SetTitle("Fit sigma vs true sigma, single site");
  grSigma_ms->SetTitle("Fit sigma vs true sigma, multi site");
  
  grDev_ss->GetXaxis()->SetTitle("resolution");
  grDev_ss->GetYaxis()->SetTitle("bias (keV)");
  
  grDev_ms->GetXaxis()->SetTitle("resolution");
  grDev_ms->GetYaxis()->SetTitle("bias (keV)");
  
  grSigma_ss->GetXaxis()->SetTitle("#sigma_{true}");
  grSigma_ss->GetYaxis()->SetTitle("#sigma_{fit}");
  
  grSigma_ms->GetXaxis()->SetTitle("#sigma_{true}");
  grSigma_ms->GetYaxis()->SetTitle("#sigma_{fit}");
  
  grDev_ss->SetMarkerStyle(20);
  grDev_ss->SetMarkerSize(0.8);
  
  grDev_ms->SetMarkerStyle(20);
  grDev_ms->SetMarkerSize(0.8);
  
  grSigma_ss->SetMarkerStyle(20);
  grSigma_ss->SetMarkerSize(0.8);
  
  grSigma_ms->SetMarkerStyle(20);
  grSigma_ms->SetMarkerSize(0.8);
  
  c1->cd(1);
  grDev_ss->Draw("AZP");
  
  c1->cd(2);
  grDev_ms->Draw("AZP");
  
  c1->cd(3);
  grSigma_ss->Draw("AZP");
  
  c1->cd(4);
  grSigma_ms->Draw("AZP");
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionTh228.ps");
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionTh228.ps]");
  
  TFile *fOut = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvalFitFunctionTh228.root","RECREATE");
  grDev_ss->Write();
  grDev_ms->Write();
  grSigma_ss->Write();
  grSigma_ms->Write();
  fOut->Close();
  
  return;
}