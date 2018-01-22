void GenerateToyMC_Co60()
{
  using namespace RooFit;
  
  gSystem->Load("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/lib/libEXOFitting.so");
  
  EXOFitter *fitter = new EXOFitter(true);
  fitter->LoadWorkspace("/nfs/slac/g/exo/maweber/EXO_Fitting_050412/analysis/EXO_Workspace_floatEres.root");
  //fitter->BuildRotationModel();
  
  RooWorkspace *w = fitter->GetWorkspace();
  
  RooRealVar *energy_ss = (RooRealVar*)w->var("energy_ss");
  RooRealVar *energy_ms = (RooRealVar*)w->var("energy_ms");
  
  //energy_ss->setRange(100,3500);
  //energy_ms->setRange(100,3500);
  
  // open file with MC histograms
  TFile *f = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/MC_Source_Histos.root","READ");
  
// #### Co60 ####################################################################
  
  // get the histograms
  TH1F *hCo60_ss = (TH1F*)f->Get("SourceP4_px_60Co_ss");
  TH1F *hCo60_ms = (TH1F*)f->Get("SourceP4_px_60Co_ms");
  
  // generate RooDataHist
  RooDataHist hCo60_d_ss("hCo60_d_ss","hCo60_d_ss",*energy_ss,hCo60_ss);
  RooDataHist hCo60_d_ms("hCo60_d_ms","hCo60_d_ms",*energy_ms,hCo60_ms);
  
  // generate RooHistPdf
  RooHistPdf *Co60_Pdf_ss = new RooHistPdf("Co60_Pdf_ss","Co60_Pdf_ss",*energy_ss,hCo60_d_ss,0);
  RooHistPdf *Co60_Pdf_ms = new RooHistPdf("Co60_Pdf_ms","Co60_Pdf_ms",*energy_ms,hCo60_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Co60_Pdf_s_ss = fitter->SmearPdf(Co60_Pdf_ss,energy_ss);
  RooAbsPdf *Co60_Pdf_s_ms = fitter->SmearPdf(Co60_Pdf_ms,energy_ms);
  
  // set the resolution parameters
  RooRealVar res0_ss("r0_ss","r0_ss",0.0);
  RooRealVar res1_ss("r1_ss","r1_ss",35.2857);
  RooRealVar res2_ss("r2_ss","r2_ss",8.60865e-03);
  
  RooRealVar res0_ms("r0_ms","r0_ms",0.0);
  RooRealVar res1_ms("r1_ms","r1_ms",39.7640);
  RooRealVar res2_ms("r2_ms","r2_ms",9.12999e-03);
  
  RooArgSet* comp_ss = Co60_Pdf_s_ss->getComponents();
  EXOresolution *res_ss = comp_ss->find("res_ss");
  //res_ss->setA(res0_ss);
  //res_ss->setB(res1_ss);
  //res_ss->setC(res2_ss);
  ((RooRealVar*)res_ss->getA())->setVal(0.0);
  ((RooRealVar*)res_ss->getB())->setVal(35.2857);
  ((RooRealVar*)res_ss->getC())->setVal(8.60865e-03);
  
  RooArgSet* comp_ms = Co60_Pdf_s_ms->getComponents();
  EXOresolution *res_ms = comp_ms->find("res_ms");
  //res_ms->setA(res0_ms);
  //res_ms->setB(res1_ms);
  //res_ms->setC(res2_ms);
  ((RooRealVar*)res_ms->getA())->setVal(0.0);
  ((RooRealVar*)res_ms->getB())->setVal(39.7640);
  ((RooRealVar*)res_ms->getC())->setVal(9.12999e-03);
  
  // generate toy MC
  RooDataSet *Co60_ToyMC_ss = Co60_Pdf_s_ss->generate(*energy_ss,1000000);
  RooDataSet *Co60_ToyMC_ms = Co60_Pdf_s_ms->generate(*energy_ms,1000000);
  
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
  
  return;
}
