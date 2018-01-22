void GenerateToyMC_Th228_Charge()
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
  
// #### Th228 ####################################################################
  
  // get the histograms
  TH1F *hTh228_ss = (TH1F*)f->Get("SourceP4_px_228Th_ss");
  TH1F *hTh228_ms = (TH1F*)f->Get("SourceP4_px_228Th_ms");
  
  // generate RooDataHist
  RooDataHist hTh228_d_ss("hTh228_d_ss","hTh228_d_ss",*energy_ss,hTh228_ss);
  RooDataHist hTh228_d_ms("hTh228_d_ms","hTh228_d_ms",*energy_ms,hTh228_ms);
  
  // generate RooHistPdf
  RooHistPdf *Th228_Pdf_ss = new RooHistPdf("Th228_Pdf_ss","Th228_Pdf_ss",*energy_ss,hTh228_d_ss,0);
  RooHistPdf *Th228_Pdf_ms = new RooHistPdf("Th228_Pdf_ms","Th228_Pdf_ms",*energy_ms,hTh228_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Th228_Pdf_s_ss = fitter->SmearPdf(Th228_Pdf_ss,energy_ss);
  RooAbsPdf *Th228_Pdf_s_ms = fitter->SmearPdf(Th228_Pdf_ms,energy_ms);
  
  // set the resolution parameters
  RooRealVar res0_ss("r0_ss","r0_ss",1.79413);
  RooRealVar res1_ss("r1_ss","r1_ss",9.97843);
  RooRealVar res2_ss("r2_ss","r2_ss",0.0);
  
  RooRealVar res0_ms("r0_ms","r0_ms",1.71945);
  RooRealVar res1_ms("r1_ms","r1_ms",2.49290);
  RooRealVar res2_ms("r2_ms","r2_ms",2.39957);
  
  RooArgSet* comp_ss = Th228_Pdf_s_ss->getComponents();
  EXOresolution *res_ss = comp_ss->find("res_ss");
  //res_ss->setA(res0_ss);
  //res_ss->setB(res1_ss);
  //res_ss->setC(res2_ss);
  ((RooRealVar*)res_ss->getA())->setVal(1.79413);
  ((RooRealVar*)res_ss->getB())->setVal(9.97843);
  ((RooRealVar*)res_ss->getC())->setVal(0.0);
  
  RooArgSet* comp_ms = Th228_Pdf_s_ms->getComponents();
  EXOresolution *res_ms = comp_ms->find("res_ms");
  //res_ms->setA(res0_ms);
  //res_ms->setB(res1_ms);
  //res_ms->setC(res2_ms);
  ((RooRealVar*)res_ms->getA())->setVal(1.71945);
  ((RooRealVar*)res_ms->getB())->setVal(2.49290);
  ((RooRealVar*)res_ms->getC())->setVal(2.39957);
  
  // generate toy MC
  RooDataSet *Th228_ToyMC_ss = Th228_Pdf_s_ss->generate(*energy_ss,1000000);
  RooDataSet *Th228_ToyMC_ms = Th228_Pdf_s_ms->generate(*energy_ms,1000000);
  
  TTree *tTh228_ss = Th228_ToyMC_ss->tree();
  tTh228_ss->SetName("tTh228_ss");
  
  TTree *tTh228_ms = Th228_ToyMC_ms->tree();
  tTh228_ms->SetName("tTh228_ms");
  
  TFile *fTh228_ss = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Th228_ToyMC_Charge_ss.root","RECREATE");
  tTh228_ss->Write();
  fTh228_ss->Close();
  
  TFile *fTh228_ms = new TFile("/nfs/slac/g/exo/maweber/EXO200Analysis/Fitting/scripts/Th228_ToyMC_Charge_ms.root","RECREATE");
  tTh228_ms->Write();
  fTh228_ms->Close();
  
// ###############################################################################
  
  return;
}
