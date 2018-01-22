void GenerateToyMC_Cs137()
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
  
// #### Cs137 ####################################################################
  
  // get the histograms
  TH1F *hCs137_ss = (TH1F*)f->Get("SourceP4_px_137Cs_ss");
  TH1F *hCs137_ms = (TH1F*)f->Get("SourceP4_px_137Cs_ms");
  
  // generate RooDataHist
  RooDataHist hCs137_d_ss("hCs137_d_ss","hCs137_d_ss",*energy_ss,hCs137_ss);
  RooDataHist hCs137_d_ms("hCs137_d_ms","hCs137_d_ms",*energy_ms,hCs137_ms);
  
  // generate RooHistPdf
  RooHistPdf *Cs137_Pdf_ss = new RooHistPdf("Cs137_Pdf_ss","Cs137_Pdf_ss",*energy_ss,hCs137_d_ss,0);
  RooHistPdf *Cs137_Pdf_ms = new RooHistPdf("Cs137_Pdf_ms","Cs137_Pdf_ms",*energy_ms,hCs137_d_ms,0);
  
  // smear PDF
  RooAbsPdf *Cs137_Pdf_s_ss = fitter->SmearPdf(Cs137_Pdf_ss,energy_ss);
  RooAbsPdf *Cs137_Pdf_s_ms = fitter->SmearPdf(Cs137_Pdf_ms,energy_ms);
  
  RooRealVar res0_ss("r0_ss","r0_ss",0.0);
  RooRealVar res1_ss("r1_ss","r1_ss",35.2857);
  RooRealVar res2_ss("r2_ss","r2_ss",8.60865e-03);
  
  RooRealVar res0_ms("r0_ms","r0_ms",0.0);
  RooRealVar res1_ms("r1_ms","r1_ms",39.7640);
  RooRealVar res2_ms("r2_ms","r2_ms",9.12999e-03);
  
  RooArgSet* comp_ss = Cs137_Pdf_s_ss->getComponents();
  EXOresolution *res_ss = comp_ss->find("res_ss");
  //res_ss->setA(res0_ss);
  //res_ss->setB(res1_ss);
  //res_ss->setC(res2_ss);
  ((RooRealVar*)res_ss->getA())->setVal(0.0);
  ((RooRealVar*)res_ss->getB())->setVal(35.2857);
  ((RooRealVar*)res_ss->getC())->setVal(8.60865e-03);
  
  RooArgSet* comp_ms = Cs137_Pdf_s_ms->getComponents();
  EXOresolution *res_ms = comp_ms->find("res_ms");
  //res_ms->setA(res0_ms);
  //res_ms->setB(res1_ms);
  //res_ms->setC(res2_ms);
  ((RooRealVar*)res_ms->getA())->setVal(0.0);
  ((RooRealVar*)res_ms->getB())->setVal(39.7640);
  ((RooRealVar*)res_ms->getC())->setVal(9.12999e-03);
  
  // generate toy MC
  RooDataSet *Cs137_ToyMC_ss = Cs137_Pdf_s_ss->generate(*energy_ss,1000000);
  RooDataSet *Cs137_ToyMC_ms = Cs137_Pdf_s_ms->generate(*energy_ms,1000000);
  
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
