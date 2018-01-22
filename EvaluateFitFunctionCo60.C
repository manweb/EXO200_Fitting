double fitFunctionCo60(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = A1_gaus*par[3];
  double A2_gaus = A1_gaus*par[4];
  double E2 = par[5];
  double sigma2 = par[6];
  double A2_erf = A2_gaus*par[3];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  double gauss2 = A2_gaus * TMath::Gaus(x[0],E2,sigma2);
  double erf2 = A2_erf * 0.5 * TMath::Erfc((x[0] - E2) / (TMath::Sqrt(2)*sigma2));
  
  return gauss1 + erf1 + gauss2 + erf2;
}

void EvaluateFitFunctionCo60()
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
  TFile *f = new TFile("MC_Source_Histos.root","READ");
  
  // get the histograms
  TH1F *hCo60_ss = (TH1F*)f->Get("SourceP4_px_60Co_ss");
  TH1F *hCo60_ms = (TH1F*)f->Get("SourceP4_px_60Co_ms");
  
  // generate RooDataHist
  RooDataHist hCo60_d_ss("hCo60_d_ss","hCo60_d_ss",*energy_ss,hCo60_ss);
  RooDataHist hCo60_d_ms("hCo60_d_ms","hCo60_d_ms",*energy_ms,hCo60_ms);
  
  // generate RooHistPdf
  RooHistPdf Co60_Pdf_ss("Co60_Pdf_ss","Co60_Pdf_ss",*energy_ss,hCo60_d_ss,0);
  RooHistPdf Co60_Pdf_ms("Co60_Pdf_ms","Co60_Pdf_ms",*energy_ms,hCo60_d_ms,0);
  
  //RooHistPdf *Co60_Pdf_ss = (RooHistPdf*)w->pdf("SourceP4_px_60Co_Pdf_ss");
  //RooHistPdf *Co60_Pdf_ms = (RooHistPdf*)w->pdf("SourceP4_px_60Co_Pdf_ms");
  
  // smear PDF
  RooAbsPdf *Co60_Pdf_s_ss = fitter->SmearPdf(&Co60_Pdf_ss,energy_ss);
  RooAbsPdf *Co60_Pdf_s_ms = fitter->SmearPdf(&Co60_Pdf_ms,energy_ms);
  
  // initialize arrays
  const int n = 10;
  double *resolution_ss = new double[n];
  double *resolution_ms = new double[n];
  double *Co60_Peak_ss = new double[n];
  double *Co60_Peak_ms = new double[n];
  
  double *Co60_Peak_err_ss = new double[n];
  double *Co60_Peak_err_ms = new double[n];
  
  double *sigma_ss = new double[n];
  double *sigma_ms = new double[n];
  double *Co60_Sigma_ss = new double[n];
  double *Co60_Sigma_ms = new double[n];
  
  double *Co60_Sigma_err_ss = new double[n];
  double *Co60_Sigma_err_ms = new double[n];
  
  // generat fit functions
  TF1 *Co60_fit_ss = new TF1("Co60_fit_ss",fitFunctionCo60,1000,1500,7);
  TF1 *Co60_fit_ms = new TF1("Co60_fit_ms",fitFunctionCo60,1000,1500,7);
  
  TF1 *erf1_ss = new TF1("erf1_ss","[3]*[0]*0.5*TMath::Erfc((x - [1]) / (TMath::Sqrt(2) * [2]))",1000,1500);
  TF1 *gaus1_ss = new TF1("gaus1_ss","gaus",1000,1500);
  
  TF1 *erf1_ms = new TF1("erf1_ms","[3]*[0]*0.5*TMath::Erfc((x - [1]) / (TMath::Sqrt(2) * [2]))",1000,1500);
  TF1 *gaus1_ms = new TF1("gaus1_ms","gaus",1000,1500);
  
  TF1 *erf2_ss = new TF1("erf2_ss","[3]*[0]*[4]*0.5*TMath::Erfc((x - [5]) / (TMath::Sqrt(2) * [6]))",1000,1500);
  TF1 *gaus2_ss = new TF1("gaus2_ss","[0]*[4]*TMath::Gaus(x,[5],[6])",1000,1500);
  
  TF1 *erf2_ms = new TF1("erf2_ms","[3]*[0]*[4]*0.5*TMath::Erfc((x - [5]) / (TMath::Sqrt(2) * [6]))",1000,1500);
  TF1 *gaus2_ms = new TF1("gaus2_ms","[0]*[4]*TMath::Gaus(x,[5],[6])",1000,1500);
  
  Co60_fit_ss->SetLineColor(kBlue);
  Co60_fit_ms->SetLineColor(kBlue);
  
  Co60_fit_ss->SetLineWidth(2);
  Co60_fit_ms->SetLineWidth(2);
  
  erf1_ss->SetLineColor(kBlue);
  gaus1_ss->SetLineColor(kBlue);
  
  erf1_ss->SetLineWidth(1);
  gaus1_ss->SetLineWidth(1);
  
  erf1_ss->SetLineStyle(2);
  gaus1_ss->SetLineStyle(2);
  
  erf1_ms->SetLineColor(kBlue);
  gaus1_ms->SetLineColor(kBlue);
  
  erf1_ms->SetLineWidth(1);
  gaus1_ms->SetLineWidth(1);
  
  erf1_ms->SetLineStyle(2);
  gaus1_ms->SetLineStyle(2);
  
  erf2_ss->SetLineColor(kBlue);
  gaus2_ss->SetLineColor(kBlue);
  
  erf2_ss->SetLineWidth(1);
  gaus2_ss->SetLineWidth(1);
  
  erf2_ss->SetLineStyle(2);
  gaus2_ss->SetLineStyle(2);
  
  erf2_ms->SetLineColor(kBlue);
  gaus2_ms->SetLineColor(kBlue);
  
  erf2_ms->SetLineWidth(1);
  gaus2_ms->SetLineWidth(1);
  
  erf2_ms->SetLineStyle(2);
  gaus2_ms->SetLineStyle(2);
  
  Co60_fit_ss->SetParameters(1000,1173,40,0.2,0.8,1333,40);
  Co60_fit_ms->SetParameters(1000,1173,40,0.2,0.8,1333,40);
  
  Co60_fit_ss->SetParLimits(0,0.0,1e6);
  Co60_fit_ss->SetParLimits(1,1050,1200);
  Co60_fit_ss->SetParLimits(2,0.0,200.0);
  Co60_fit_ss->SetParLimits(3,0.01,0.8);
  Co60_fit_ss->SetParLimits(4,0.5,1.2);
  Co60_fit_ss->SetParLimits(5,1280,1400);
  Co60_fit_ss->SetParLimits(6,0.0,200.0);
  
  Co60_fit_ms->SetParLimits(0,0.0,1e6);
  Co60_fit_ms->SetParLimits(1,1050,1200);
  Co60_fit_ms->SetParLimits(2,0.0,200.0);
  Co60_fit_ms->SetParLimits(3,0.01,0.8);
  Co60_fit_ms->SetParLimits(4,0.5,1.2);
  Co60_fit_ms->SetParLimits(5,1280,1400);
  Co60_fit_ms->SetParLimits(6,0.0,200.0);
  
  RooRealVar res0_ss("r0_ss","r0_ss",0.0);
  RooRealVar res1_ss("r1_ss","r1_ss",0.0);
  RooRealVar res2_ss("r2_ss","r2_ss",0.0);
  
  RooRealVar res0_ms("r0_ms","r0_ms",0.0);
  RooRealVar res1_ms("r1_ms","r1_ms",0.0);
  RooRealVar res2_ms("r2_ms","r2_ms",0.0);
  
  RooArgSet* comp_ss = Co60_Pdf_s_ss->getComponents();
  EXOresolution *res_ss = comp_ss->find("res_ss");
  ((RooRealVar*)res_ss->getA())->setVal(0.0);
  ((RooRealVar*)res_ss->getB())->setVal(0.0);
  ((RooRealVar*)res_ss->getC())->setVal(0.0);
  
  RooArgSet* comp_ms = Co60_Pdf_s_ms->getComponents();
  EXOresolution *res_ms = comp_ms->find("res_ms");
  ((RooRealVar*)res_ms->getA())->setVal(0.0);
  ((RooRealVar*)res_ms->getB())->setVal(0.0);
  ((RooRealVar*)res_ms->getC())->setVal(0.0);
  
  TH1F *h_ss = new TH1F("h_ss","h_ss",200,500,2500);
  TH1F *h_ms = new TH1F("h_ms","h_ms",200,500,2500);
  
  double rStart = 0.012;
  double rEnd = 0.045;
  double stp = (rEnd - rStart) / double(n);
  
  RooPlot *frame_ss;
  RooPlot *frame_ms;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
  c1->Divide(2,2);
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionCo60.ps[");
  
  for (int i = 0; i < n; i++) {
    double r = i*stp + rStart;
    ((RooRealVar*)res_ss->getA())->setVal(r*TMath::Sqrt(2615.0));
    ((RooRealVar*)res_ms->getA())->setVal(r*TMath::Sqrt(2615.0));
    
    // generate toy MC
    RooDataSet *Co60_ToyMC_ss = Co60_Pdf_s_ss->generate(*energy_ss,100000);
    RooDataSet *Co60_ToyMC_ms = Co60_Pdf_s_ms->generate(*energy_ms,100000);
    
    ((TTree*)Co60_ToyMC_ss->tree())->Draw("energy_ss>>h_ss","","goff");
    h_ss->Fit("Co60_fit_ss","rn");
    
    ((TTree*)Co60_ToyMC_ms->tree())->Draw("energy_ms>>h_ms","","goff");
    h_ms->Fit("Co60_fit_ms","rn");
    
    h_ss->SetTitle(Form("resolution = %.4f",r));
    h_ms->SetTitle(Form("resolution = %.4f",r));
    
    resolution_ss[i] = r;
    resolution_ms[i] = r;
    Co60_Peak_ss[i] = 1333.0 - Co60_fit_ss->GetParameter(5);
    Co60_Peak_ms[i] = 1333.0 - Co60_fit_ms->GetParameter(5);
    
    Co60_Peak_err_ss[i] = Co60_fit_ss->GetParError(5);
    Co60_Peak_err_ms[i] = Co60_fit_ms->GetParError(5);
    
    sigma_ss[i] = r*TMath::Sqrt(2615.0*1333.0);
    sigma_ms[i] = r*TMath::Sqrt(2615.0*1333.0);
    Co60_Sigma_ss[i] = Co60_fit_ss->GetParameter(6);
    Co60_Sigma_ms[i] = Co60_fit_ms->GetParameter(6);
    
    Co60_Sigma_err_ss[i] = Co60_fit_ss->GetParError(6);
    Co60_Sigma_err_ms[i] = Co60_fit_ss->GetParError(6);
    
    frame_ss = energy_ss->frame();
    Co60_Pdf_ss.plotOn(frame_ss,LineStyle(2),LineWidth(1));
    Co60_Pdf_s_ss->plotOn(frame_ss);
    frame_ss->SetTitle(Form("Single Site, resolution = %.4f",r));
    
    frame_ms = energy_ms->frame();
    Co60_Pdf_ms.plotOn(frame_ms,LineStyle(2),LineWidth(1));
    Co60_Pdf_s_ms->plotOn(frame_ms);
    frame_ms->SetTitle(Form("Multi Site, resolution = %.4f",r));
    
    erf1_ss->SetParameters(Co60_fit_ss->GetParameters());
    gaus1_ss->SetParameters(Co60_fit_ss->GetParameters());
    
    erf1_ms->SetParameters(Co60_fit_ms->GetParameters());
    gaus1_ms->SetParameters(Co60_fit_ms->GetParameters());
    
    erf2_ss->SetParameters(Co60_fit_ss->GetParameters());
    gaus2_ss->SetParameters(Co60_fit_ss->GetParameters());
    
    erf2_ms->SetParameters(Co60_fit_ms->GetParameters());
    gaus2_ms->SetParameters(Co60_fit_ms->GetParameters());
    
    c1->cd(1);
    frame_ss->Draw();
    
    c1->cd(2);
    frame_ms->Draw();
    
    c1->cd(3);
    h_ss->Draw("EZP");
    Co60_fit_ss->Draw("same");
    erf1_ss->Draw("same");
    gaus1_ss->Draw("same");
    erf2_ss->Draw("same");
    gaus2_ss->Draw("same");
    
    c1->cd(4);
    h_ms->Draw("EZP");
    Co60_fit_ms->Draw("same");
    erf1_ms->Draw("same");
    gaus1_ms->Draw("same");
    erf2_ms->Draw("same");
    gaus2_ms->Draw("same");
    
    c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionCo60.ps");
    
    //h_ss->Delete();
    //h_ms->Delete();
  }
  
  TGraphErrors *grDev_ss = new TGraphErrors(n,resolution_ss,Co60_Peak_ss,0,Co60_Peak_err_ss);
  TGraphErrors *grDev_ms = new TGraphErrors(n,resolution_ms,Co60_Peak_ms,0,Co60_Peak_err_ms);
  TGraphErrors *grSigma_ss = new TGraphErrors(n,sigma_ss,Co60_Sigma_ss,0,Co60_Sigma_err_ss);
  TGraphErrors *grSigma_ms = new TGraphErrors(n,sigma_ms,Co60_Sigma_ms,0,Co60_Sigma_err_ms);
  
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
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionCo60.ps");
  
  c1->Print("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvaluateFitFunctionCo60.ps]");
  
  TFile *fOut = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/EvalFitFunctionCo60.root","RECREATE");
  grDev_ss->Write();
  grDev_ms->Write();
  grSigma_ss->Write();
  grSigma_ms->Write();
  fOut->Close();
  
  return;
}