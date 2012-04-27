analysis.C
==========
#define analysis_cxx
#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>

void analysis::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L analysis.C
//      Root > analysis t
//      Root > t.Loop();       // Loop on all entries
//

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  TLorentzVector par1,par2,par3,par4;
  TH1F *diMuonMass = new TH1F("diMuonMass","DiMuon Mass", 100,0.,100.);
  TH1F *diMuonDEta = new TH1F("diMuonDEta","DiMuon Eta",100,-10.,10.);
  TH1F *diMuonDPhi = new TH1F("diMuonDPhi","DiMuon Phi",20,-15.,15.);
  TH1F *diMuonDR   = new TH1F("diMuonDR","DiMuon Phi",20,-5.,5.);
  TH1F *diMuonN    = new TH1F("diMuonN","DiMuon N",20,-5.,5.);
  diMuonN->SetLineColor(9);
  //
  TH1F *diMuonDPhi_test = new TH1F("diMuonDPhi_test","DiMuon Phi", 20,-15.,15.);
  diMuonDPhi_test->SetLineColor(2);
  TH1F *diMuonDR_test = new TH1F("diMuonDR_test","DiMuon Eta", 20,-5.,5.);
  diMuonDR_test->SetLineColor(2);
  //
  //Electron
  TH1F *diElectronMass = new TH1F("diElectronMass","DiElectron Mass", 100,0.,100.);
  TH1F *diElectronDEta = new TH1F("diElectronDEta","DiElectron Eta", 100,-10.,10.);
  TH1F *diElectronDPhi = new TH1F("diElectronDPhi","DiElectron Phi", 20,-15.,15.);
  TH1F *diElectronDR   = new TH1F("diElectronDR","DiElectron Phi", 20,-5.,5.);
  TH1F *diElectronN    = new TH1F("diElectronN","DiElectron N", 20,-5.,5.);
  diElectronN->SetLineColor(9);
  //
  TH1F *diElectronDPhi_test = new TH1F("diElectronDPhi_test","DiElectron Phi", 20,-15.,15.);
  diElectronDPhi_test->SetLineColor(2);
  TH1F *diElectronDR_test   = new TH1F("diElectronDR_test","DiElectron Eta" ,20, -5.,5.);
  diElectronDR_test->SetLineColor(2);
  //
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //General cuts (Try to change to other physics objects)
    
    diMuonN->Fill(nMuon);

    if(nMuon!=2) continue;
    if(muonIsGlobal[0]!=true || muonIsGlobal[1]!=true) continue;


    //Set Four-Vectors
    par1.SetPxPyPzE(muonPx[0],muonPy[0],muonPz[0],muonE[0]);
    par2.SetPxPyPzE(muonPx[1],muonPy[1],muonPz[1],muonE[1]);
    par3.SetPxPyPzE(electronPx[0],electronPy[0],electronPz[0],electronE[0]);
    par4.SetPxPyPzE(electronPx[1],electronPy[1],electronPz[1],electronE[1]);
    //Fill Histogram
    diMuonDEta->Fill(muonEta[0]-muonEta[1]);
    diMuonDPhi->Fill(par1.DeltaPhi(par2));
    diMuonDR->Fill(par1.DeltaR(par2));
    diMuonMass->Fill((par1+par2).Mag());
    //Fill histogram electron
    diElectronDEta->Fill(electronEta[0]-electronEta[1]);
    diElectronDPhi->Fill(par3.DeltaPhi(par4));
    diElectronDR->Fill(par3.DeltaR(par4));
    diElectronMass->Fill((par3+par4).Mag());
    diElectronN->Fill(nElectron);

    //Fill testing histograms
    float Dphi = muonPhi[0]-muonPhi[1];
    if(fabs(Dphi) > TMath::Pi()){
      if(Dphi < TMath::Pi())
	Dphi=Dphi+2.0*TMath::Pi();
      else Dphi=Dphi-2.0*TMath::Pi();
      }
    //
    float Dphi2 = electronPhi[0]-electronPhi[1];
    if(fabs(Dphi2) > TMath::Pi()){
      if(Dphi2 < TMath::Pi())
	Dphi2=Dphi2+2.0*TMath::Pi();
      else Dphi2=Dphi2-2.0*TMath::Pi();
      }

    float DR_muon = TMath::Sqrt(TMath::Power(Dphi,2)+TMath::Power(muonEta[0]-muonEta[1],2));
    diMuonDPhi_test->Fill(Dphi);
    diMuonDR_test->Fill(DR_muon);
    
    float DR_electron = TMath::Sqrt(TMath::Power(Dphi2,2)+TMath::Power(electronEta[0]-electronEta[1],2));
    diElectronDPhi_test->Fill(Dphi2);
    diElectronDR_test->Fill(DR_electron);

    //OL
    /* int a;
    float muonOL[a];
    if(DR_muon < 0.1)
    */


    //

  }
  //Draw

  TCanvas *c1 = new TCanvas("c1","c1",1500,600);
  c1->Divide(5,2);
  c1->cd(1);
  diMuonDEta->Draw();
  c1->cd(2);
  diMuonDPhi->Draw();
  diMuonDPhi_test->Draw("same");
  c1->cd(3);
  diMuonDR->Draw();
  diMuonDR_test->Draw("same");
  c1->cd(4);
  diMuonMass->Draw();
  c1->cd(5);
  diMuonN->Draw();
  //electron
  c1->cd(6);
  diElectronDEta->Draw();
  c1->cd(7);
  diElectronDPhi->Draw();
  diElectronDPhi_test->Draw("same");
  c1->cd(8);
  diElectronDR->Draw();
  diElectronDR_test->Draw("same");
  c1->cd(9);
  diElectronMass->Draw();
  c1->cd(10);
  diElectronN->Draw();

  c1->Update();

   
   
}
