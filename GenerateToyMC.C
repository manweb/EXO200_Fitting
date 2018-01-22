void GenerateToyMC_Th228()
{
  using namespace RooFit;
  
  gSystem->Load("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/lib/libEXOFitting.so");
  
  EXOFitter *fitter = new EXOFitter(true);
  fitter->BuildWorkspace();
  fitter->BuildRotationModel();
  
  RooWorkspace *w = fitter->GetWorkspace();
  
  RooRealVar *energy_ss = (RooRealVar*)w->var("energy_ss");
  RooRealVar *energy_ms = (RooRealVar*)w->var("energy_ms");
  
  energy_ss->setRange(100,3500);
  energy_ms->setRange(100,3500);
  
  // open file with MC histograms
  TFile *f = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/MC_Source_Histos.root","READ");
  
// #### Th228 ####################################################################
  
  // get the histograms
  TH1F *hTh228_ss = (TH1F*)f->Get("SourceP4_px_228Th_ss");
  TH1F *hTh228_ms = (TH1F*)f->Get("SourceP4_px_228Th_ms");
  
  // generate RooDataHist
  RooDataHist hTh228_d_ss("hTh228_d_ss","hTh228_d_ss",*energy_ss,hTh228_ss);
  RooDataHist hTh228_d_ms("hTh228_d_ms","hTh228_d_ms",*energy_ms,hTh228_ms);
  
  // generate RooHistPdf
  RooHistPdf Th228_Pdf_ss("Th228_Pdf_ss","Th228_Pdf_ss",*energy_ss,hTh228_d_ss,0);
  RooHistPdf Th228_Pdf_ms("Th228_Pdf_ms","Th228_Pdf_ms",*energy_ms,hTh228_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Th228_Pdf_s_ss = fitter->SmearPdf(&Th228_Pdf_ss,energy_ss);
  RooAbsPdf *Th228_Pdf_s_ms = fitter->SmearPdf(&Th228_Pdf_ms,energy_ms);
  
  // generate toy MC
  RooDataSet *Th228_ToyMC_ss = Th228_Pdf_s_ss->generate(*energy_ss,100000);
  RooDataSet *Th228_ToyMC_ms = Th228_Pdf_s_ms->generate(*energy_ms,100000);
  
  TTree *tTh228_ss = Th228_ToyMC_ss->tree();
  tTh228_ss->SetName("tTh228_ss");
  
  TTree *tTh228_ms = Th228_ToyMC_ms->tree();
  tTh228_ms->SetName("tTh228_ms");
  
  TFile *fTh228_ss = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Th228_ToyMC_ss.root","RECREATE");
  tTh228_ss->Write();
  fTh228_ss->Close();
  
  TFile *fTh228_ms = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Th228_ToyMC_ms.root","RECREATE");
  tTh228_ms->Write();
  fTh228_ms->Close();
  
// ###############################################################################
  
// #### Co60 ####################################################################
  
  // get the histograms
  TH1F *hCo60_ss = (TH1F*)f->Get("SourceP4_px_60Co_ss");
  TH1F *hCo60_ms = (TH1F*)f->Get("SourceP4_px_60Co_ms");
  
  // generate RooDataHist
  RooDataHist hCo60_d_ss("hCo60_d_ss","hCo60_d_ss",*energy_ss,hCo60_ss);
  RooDataHist hCo60_d_ms("hCo60_d_ms","hCo60_d_ms",*energy_ms,hCo60_ms);
  
  // generate RooHistPdf
  RooHistPdf Co60_Pdf_ss("Co60_Pdf_ss","Co60_Pdf_ss",*energy_ss,hCo60_d_ss,0);
  RooHistPdf Co60_Pdf_ms("Co60_Pdf_ms","Co60_Pdf_ms",*energy_ms,hCo60_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Co60_Pdf_s_ss = fitter->SmearPdf(&Co60_Pdf_ss,energy_ss);
  RooAbsPdf *Co60_Pdf_s_ms = fitter->SmearPdf(&Co60_Pdf_ms,energy_ms);
  
  // generate toy MC
  RooDataSet *Co60_ToyMC_ss = Co60_Pdf_s_ss->generate(*energy_ss,100000);
  RooDataSet *Co60_ToyMC_ms = Co60_Pdf_s_ms->generate(*energy_ms,100000);
  
  TTree *tCo60_ss = Co60_ToyMC_ss->tree();
  tCo60_ss->SetName("tCo60_ss");
  
  TTree *tCo60_ms = Co60_ToyMC_ms->tree();
  tCo60_ms->SetName("tCo60_ms");
  
  TFile *fCo60_ss = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Co60_ToyMC_ss.root","RECREATE");
  tCo60_ss->Write();
  fCo60_ss->Close();
  
  TFile *fCo60_ms = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Co60_ToyMC_ms.root","RECREATE");
  tCo60_ms->Write();
  fCo60_ms->Close();
  
// ###############################################################################
  
// #### Cs137 ####################################################################
  
  // get the histograms
  TH1F *hCs137_ss = (TH1F*)f->Get("SourceP4_px_137Cs_ss");
  TH1F *hCs137_ms = (TH1F*)f->Get("SourceP4_px_137Cs_ms");
  
  // generate RooDataHist
  RooDataHist hCs137_d_ss("hCs137_d_ss","hCs137_d_ss",*energy_ss,hCs137_ss);
  RooDataHist hCs137_d_ms("hCs137_d_ms","hCs137_d_ms",*energy_ms,hCs137_ms);
  
  // generate RooHistPdf
  RooHistPdf Cs137_Pdf_ss("Cs137_Pdf_ss","Cs137_Pdf_ss",*energy_ss,hCs137_d_ss,0);
  RooHistPdf Cs137_Pdf_ms("Cs137_Pdf_ms","Cs137_Pdf_ms",*energy_ms,hCs137_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Cs137_Pdf_s_ss = fitter->SmearPdf(&Cs137_Pdf_ss,energy_ss);
  RooAbsPdf *Cs137_Pdf_s_ms = fitter->SmearPdf(&Cs137_Pdf_ms,energy_ms);
  
  // generate toy MC
  RooDataSet *Cs137_ToyMC_ss = Cs137_Pdf_s_ss->generate(*energy_ss,100000);
  RooDataSet *Cs137_ToyMC_ms = Cs137_Pdf_s_ms->generate(*energy_ms,100000);
  
  TTree *tCs137_ss = Cs137_ToyMC_ss->tree();
  tCs137_ss->SetName("tCs137_ss");
  
  TTree *tCs137_ms = Cs137_ToyMC_ms->tree();
  tCs137_ms->SetName("tCs137_ms");
  
  TFile *fCs137_ss = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Cs137_ToyMC_ss.root","RECREATE");
  tCs137_ss->Write();
  fCs137_ss->Close();
  
  TFile *fCs137_ms = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Cs137_ToyMC_ms.root","RECREATE");
  tCs137_ms->Write();
  fCs137_ms->Close();
  
// ###############################################################################
  
  return;
}