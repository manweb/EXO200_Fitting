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

double fitFunctionCs137(double *x, double *par)
{
  double A1_gaus = par[0];
  double E1 = par[1];
  double sigma1 = par[2];
  double A1_erf = par[3]*par[0];
  
  double gauss1 = A1_gaus * TMath::Gaus(x[0],E1,sigma1);
  double erf1 = A1_erf * 0.5 * TMath::Erfc((x[0] - E1) / (TMath::Sqrt(2)*sigma1));
  
  return gauss1 + erf1;
}

void EvaluateFitFunction()
{
  using namespace RooFit;
  
  gSystem->Load("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/lib/libEXOFitting.so");
  
  EXOFitter *fitter = new EXOFitter(true);
  fitter->BuildWorkspace();
  fitter->BuildRotationModel();
  
  RooWorkspace *w = fitter->GetWorkspace();
  
  RooRealVar *energy_ss = (RooRealVar*)w->var("energy_ss");
  RooRealVar *energy_ms = (RooRealVar*)w->var("energy_ms");
  
  RooRealVar *res0_ss = (RooRealVar*)w->var("res0_ss");
  RooRealVar *res1_ss = (RooRealVar*)w->var("res1_ss");
  RooRealVar *res2_ss = (RooRealVar*)w->var("res2_ss");
  
  RooRealVar *res0_ms = (RooRealVar*)w->var("res0_ms");
  RooRealVar *res1_ms = (RooRealVar*)w->var("res1_ms");
  RooRealVar *res2_ms = (RooRealVar*)w->var("res2_ms");
  
  energy_ss->setRange(100,3500);
  energy_ms->setRange(100,3500);
  
  res0_ss->setVal(0.0);
  res1_ss->setVal(104.0);
  res2_ss->setVal(0.0);
  
  res0_ms->setVal(0.0);
  res1_ms->setVal(104.0);
  res2_ms->setVal(0.0);
  
  // open file with MC histograms
  TFile *f = new TFile("MC_Source_Histos.root","READ");
  
  // get the histograms
  TH1F *hTh228_ss = (TH1F*)f->Get("SourceP4_px_228Th_ss");
  TH1F *hCo60_ss = (TH1F*)f->Get("SourceP4_px_60Co_ss");
  TH1F *hCs137_ss = (TH1F*)f->Get("SourceP4_px_137Cs_ss");
  
  TH1F *hTh228_ms = (TH1F*)f->Get("SourceP4_px_228Th_ms");
  TH1F *hCo60_ms = (TH1F*)f->Get("SourceP4_px_60Co_ms");
  TH1F *hCs137_ms = (TH1F*)f->Get("SourceP4_px_137Cs_ms");
  
  // generate RooDataHist
  RooDataHist hTh228_d_ss("hTh228_d_ss","hTh228_d_ss",*energy_ss,hTh228_ss);
  RooDataHist hCo60_d_ss("hCo60_d_ss","hCo60_d_ss",*energy_ss,hCo60_ss);
  RooDataHist hCs137_d_ss("hCs137_d_ss","hCs137_d_ss",*energy_ss,hCs137_ss);
  
  RooDataHist hTh228_d_ms("hTh228_d_ms","hTh228_d_ms",*energy_ms,hTh228_ms);
  RooDataHist hCo60_d_ms("hCo60_d_ms","hCo60_d_ms",*energy_ms,hCo60_ms);
  RooDataHist hCs137_d_ms("hCs137_d_ms","hCs137_d_ms",*energy_ms,hCs137_ms);
  
  // generate RooHistPdf
  RooHistPdf Th228_Pdf_ss("Th228_Pdf_ss","Th228_Pdf_ss",*energy_ss,hTh228_d_ss,0);
  RooHistPdf Co60_Pdf_ss("Co60_Pdf_ss","Co60_Pdf_ss",*energy_ss,hCo60_d_ss,0);
  RooHistPdf Cs137_Pdf_ss("Cs137_Pdf_ss","Cs137_Pdf_ss",*energy_ss,hCs137_d_ss,0);
  
  RooHistPdf Th228_Pdf_ms("Th228_Pdf_ms","Th228_Pdf_ms",*energy_ms,hTh228_d_ms,0);
  RooHistPdf Co60_Pdf_ms("Co60_Pdf_ms","Co60_Pdf_ms",*energy_ms,hCo60_d_ms,0);
  RooHistPdf Cs137_Pdf_ms("Cs137_Pdf_ms","Cs137_Pdf_ms",*energy_ms,hCs137_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Th228_Pdf_s_ss = fitter->SmearPdf(&Th228_Pdf_ss,energy_ss);
  RooAbsPdf *Co60_Pdf_s_ss = fitter->SmearPdf(&Co60_Pdf_ss,energy_ss);
  RooAbsPdf *Cs137_Pdf_s_ss = fitter->SmearPdf(&Cs137_Pdf_ss,energy_ss);
  
  RooAbsPdf *Th228_Pdf_s_ms = fitter->SmearPdf(&Th228_Pdf_ms,energy_ms);
  RooAbsPdf *Co60_Pdf_s_ms = fitter->SmearPdf(&Co60_Pdf_ms,energy_ms);
  RooAbsPdf *Cs137_Pdf_s_ms = fitter->SmearPdf(&Cs137_Pdf_ms,energy_ms);
  
  // initialize arrays
  const int n = 10;
  double *resolution_ss = new double[n];
  double *resolution_ms = new double[n];
  double *Th228_Peak_ss = new double[n];
  double *Co60_Peak1_ss = new double[n];
  double *Co60_Peak2_ss = new double[n];
  double *Cs137_Peak_ss = new double[n];
  double *Th228_Peak_ms = new double[n];
  double *Co60_Peak1_ms = new double[n];
  double *Co60_Peak2_ms = new double[n];
  double *Cs137_Peak_ms = new double[n];
  
  double *Th228_Peak_err_ss = new double[n];
  double *Co60_Peak1_err_ss = new double[n];
  double *Co60_Peak2_err_ss = new double[n];
  double *Cs137_Peak_err_ss = new double[n];
  double *Th228_Peak_err_ms = new double[n];
  double *Co60_Peak1_err_ms = new double[n];
  double *Co60_Peak2_err_ms = new double[n];
  double *Cs137_Peak_err_ms = new double[n];
  
  // generat fit functions
  TF1 *Th228_fit_ss = new TF1("Th228_fit_ss",fitFunctionTh228,2400,3200,4);
  TF1 *Co60_fit_ss = new TF1("Co60_fit_ss",fitFunctionCo60,1000,1500,7);
  TF1 *Cs137_fit_ss = new TF1("Cs137_fit_ss",fitFunctionCs137,500,800,4);
  
  TF1 *Th228_fit_ms = new TF1("Th228_fit_ms",fitFunctionTh228,2400,3200,4);
  TF1 *Co60_fit_ms = new TF1("Co60_fit_ms",fitFunctionCo60,1000,1500,7);
  TF1 *Cs137_fit_ms = new TF1("Cs137_fit_ms",fitFunctionCs137,500,800,4);
  
  Th228_fit_ss->SetLineColor(kBlue);
  Co60_fit_ss->SetLineColor(kBlue);
  Cs137_fit_ss->SetLineColor(kBlue);
  
  Th228_fit_ms->SetLineColor(kBlue);
  Co60_fit_ms->SetLineColor(kBlue);
  Cs137_fit_ms->SetLineColor(kBlue);
  
  Th228_fit_ss->SetLineWidth(2);
  Co60_fit_ss->SetLineWidth(2);
  Cs137_fit_ss->SetLineWidth(2);
  
  Th228_fit_ms->SetLineWidth(2);
  Co60_fit_ms->SetLineWidth(2);
  Cs137_fit_ms->SetLineWidth(2);
  
  Th228_fit_ss->SetParameters(1000,2615,40,0.2);
  Co60_fit_ss->SetParameters(1000,1173,40,0.2,0.8,1333,40);
  Cs137_fit_ss->SetParameters(1000,662,40,0.2);
  
  Th228_fit_ms->SetParameters(1000,2615,40,0.2);
  Co60_fit_ms->SetParameters(1000,1173,40,0.2,0.8,1333,40);
  Cs137_fit_ms->SetParameters(1000,662,40,0.2);
  
  RooRealVar res0("res0","res0",0.0);
  RooRealVar res1("res1","res1",104.0);
  RooRealVar res2("res2","res2",0.0);
  RooArgSet* comp = Co60_Pdf_s_ss->getComponents();
  EXOresolution *res = comp->find("res_ss");
  res->setB(res1);
  
  //((EXOresolution*)((RooArgSet*)Co60_Pdf_s_ss->getComponents())->find("res_ss"))->setB(RooRealVar("res1_ss","res1_ss",104.0));
  
  // generate toy MC
  //RooDataSet *Th228_ToyMC_ss = Th228_Pdf_s_ss->generate(*energy_ss,100000);
  
  //((TTree*)Th228_ToyMC_ss->tree())->Draw("energy_ss>>h_ss(100,2000,3000)","","EZP");
  //h_ss->Fit("Th228_fit_ss","r");
  
  RooPlot *frame_ss = energy_ss->frame();
  Co60_Pdf_ss.plotOn(frame_ss,LineStyle(2),LineWidth(1));
  Co60_Pdf_s_ss->plotOn(frame_ss);
  
  RooPlot *frame_ms = energy_ms->frame();
  Co60_Pdf_ms.plotOn(frame_ms,LineStyle(2),LineWidth(1));
  Co60_Pdf_s_ms->plotOn(frame_ms);
  
  TCanvas *cPdf = new TCanvas("cPdf", "PDF", 1200, 600);
  cPdf->Divide(2, 1);
  
  cPdf->cd(1);
  frame_ss->Draw();
  
  cPdf->cd(2);
  frame_ms->Draw();
}