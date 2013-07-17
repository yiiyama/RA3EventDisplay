#include "SusyEvent.h"

#include "TString.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TText.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"

#include <iostream>
#include <map>
#include <set>

double const sizeNorm(2.);

class RA3EventDisplay {
public:
  RA3EventDisplay();
  ~RA3EventDisplay();

  int addPath(TString const&);
  void setPtThreshold(double _pt) { ptThreshold_ = _pt; }
  void setPFJetCollection(TString const&);

  bool showEvent(unsigned, unsigned);
  bool showNextEvent();
  void print(TString const& _fileName) { canvas_.Print(_fileName); }

private:
  void showEvent_();

  TChain mainTree_;
  TChain scanTree_;

  double ptThreshold_;
  TString pfJetCollection_;

  unsigned runNumber_;
  unsigned luminosityBlockNumber_;
  unsigned eventNumber_;
  long currentEntry_;

  susy::Event event_;

  TCanvas canvas_;
  TH2F etaPhiField_;
  TLegend sizeLegend_;
  TLegend colorLegend_;
  TPaveText eventInfo_;

  TObjArray garbageCollection_;
};

RA3EventDisplay::RA3EventDisplay() :
  mainTree_("susyTree"),
  scanTree_("susyTree"),
  ptThreshold_(3.),
  pfJetCollection_("ak5"),
  runNumber_(0),
  luminosityBlockNumber_(0),
  eventNumber_(0),
  currentEntry_(-1),
  event_(),
  canvas_("evdisp", "RA3 Event Display"),
  etaPhiField_("evdisp", ";#eta;#phi", 100, -4., 4., 100, -TMath::Pi(), TMath::Pi()),
  sizeLegend_(0.82, 0.35, 0.95, 0.5),
  colorLegend_(0.82, 0.05, 0.95, 0.3),
  eventInfo_(0.82, 0.8, 0.95, 0.95, "brNDC"),
  garbageCollection_()
{
  mainTree_.SetBranchStatus("*", 0);
  mainTree_.SetBranchStatus("pfParticles*", 1);
  mainTree_.SetBranchStatus("pfJets*", 1);
  mainTree_.SetBranchStatus("met_pfType01CorrectedMet*", 1);

  scanTree_.SetBranchStatus("*", 0);
  scanTree_.SetBranchStatus("runNumber", 1);
  scanTree_.SetBranchStatus("luminosityBlockNumber", 1);
  scanTree_.SetBranchStatus("eventNumber", 1);
  scanTree_.SetBranchAddress("runNumber", &runNumber_);
  scanTree_.SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber_);
  scanTree_.SetBranchAddress("eventNumber", &eventNumber_);

  canvas_.SetGrid(1, 1);
  canvas_.SetMargin(0.05, 0.2, 0.05, 0.05);

  etaPhiField_.SetStats(false);
  etaPhiField_.SetTitleSize(0.15);
  etaPhiField_.GetXaxis()->SetNdivisions(20, 1, 0);
  etaPhiField_.GetYaxis()->SetNdivisions(20, 1, 0);
  etaPhiField_.GetXaxis()->SetTitleOffset(-0.3);
  etaPhiField_.GetYaxis()->SetTitleOffset(-0.3);
  etaPhiField_.GetXaxis()->SetTitleSize(0.05);
  etaPhiField_.GetYaxis()->SetTitleSize(0.05);

  TMarker* tenGeV(new TMarker(0., 0., kFullCircle));
  TMarker* hundredGeV(new TMarker(0., 0., kFullCircle));
  tenGeV->SetMarkerSize(sizeNorm);
  hundredGeV->SetMarkerSize(sizeNorm * 2.);

  sizeLegend_.SetBorderSize(0);
  sizeLegend_.SetFillStyle(0);
  sizeLegend_.SetTextAlign(32); // right, center
  sizeLegend_.AddEntry(tenGeV, "10 GeV", "P");
  sizeLegend_.AddEntry(hundredGeV, "100 GeV", "P");

  TMarker* photon(new TMarker(0., 0., kFullCircle));
  TMarker* electron(new TMarker(0., 0., kFullCircle));
  TMarker* muon(new TMarker(0., 0., kFullCircle));
  TMarker* jetHadron(new TMarker(0., 0., kFullCircle));
  TMarker* strayHadron(new TMarker(0., 0., kFullCircle));
  photon->SetMarkerColor(kGreen);
  electron->SetMarkerColor(kBlue);
  muon->SetMarkerColor(kRed);
  jetHadron->SetMarkerColor(kOrange);
  strayHadron->SetMarkerColor(kBlack);

  colorLegend_.AddEntry(photon, "Photon", "P");
  colorLegend_.AddEntry(electron, "Electron", "P");
  colorLegend_.AddEntry(muon, "Muon", "P");
  colorLegend_.AddEntry(jetHadron, "Hadron in " + pfJetCollection_ + "PFJet", "P");
  colorLegend_.AddEntry(strayHadron, "Hadron not in " + pfJetCollection_ + "PFJet", "P");

  eventInfo_.SetBorderSize(1);

  garbageCollection_.SetOwner(true);
}

RA3EventDisplay::~RA3EventDisplay()
{
  event_.releaseTree(mainTree_);

  TList* colorList(colorLegend_.GetListOfPrimitives());
  for(int iL(0); iL != colorList->GetEntries(); ++iL){
    TLegendEntry* entry(static_cast<TLegendEntry*>(colorList->At(iL)));
    delete entry->GetObject();
  }

  TList* sizeList(sizeLegend_.GetListOfPrimitives());
  for(int iL(0); iL != sizeList->GetEntries(); ++iL){
    TLegendEntry* entry(static_cast<TLegendEntry*>(sizeList->At(iL)));
    delete entry->GetObject();
  }
}

int
RA3EventDisplay::addPath(TString const& _path)
{
  bool first(mainTree_.GetNtrees() == 0);

  int result(mainTree_.Add(_path));
  scanTree_.Add(_path);

  if(first) event_.setInput(mainTree_);

  return result;
}

void
RA3EventDisplay::setPFJetCollection(TString const& _pfJetCollection)
{
  TList* colorList(colorLegend_.GetListOfPrimitives());
  for(int iL(0); iL != colorList->GetEntries(); ++iL){
    TLegendEntry* entry(static_cast<TLegendEntry*>(colorList->At(iL)));
    TString label(entry->GetLabel());
    if(label.Contains("PFJet")){
      label.ReplaceAll(pfJetCollection_, _pfJetCollection);
      entry->SetLabel(label);
    }
  }

  pfJetCollection_ = _pfJetCollection;
}

bool
RA3EventDisplay::showEvent(unsigned _runNumber, unsigned _eventNumber)
{
  long start(currentEntry_);
  while(runNumber_ != _runNumber || eventNumber_ != _eventNumber){
    ++currentEntry_;
    if(currentEntry_ == start) break;
    int bytes(scanTree_.GetEntry(currentEntry_));
    if(bytes == 0){
      currentEntry_ = 0;
      bytes = scanTree_.GetEntry(currentEntry_);
    }
    if(bytes < 0) break;
  }

  if(_runNumber == runNumber_ && _eventNumber == eventNumber_){
    showEvent_();
    return true;
  }
  else{
    std::cerr << "Run " << _runNumber << " Event " << _eventNumber << " not found in input" << std::endl;
    return false;
  }
}

bool
RA3EventDisplay::showNextEvent()
{
  if(scanTree_.GetEntry(++currentEntry_) > 0){
    showEvent_();
    return true;
  }
  else{
    std::cerr << "End of input" << std::endl;
    return false;
  }
}

void
RA3EventDisplay::showEvent_()
{
  garbageCollection_.Clear();
  canvas_.Clear();

  etaPhiField_.SetTitle(TString::Format("Run %d, Lumi %d, Event %d", runNumber_, luminosityBlockNumber_, eventNumber_));
  etaPhiField_.Draw();

  sizeLegend_.Draw();
  colorLegend_.Draw();

  event_.getEntry(currentEntry_);

  susy::PFParticleCollection const& particles(event_.pfParticles);
  susy::PFJetCollection const& pfJets(event_.pfJets[pfJetCollection_]);
  susy::MET const& met(event_.metMap["pfType01CorrectedMet"]);

  unsigned nPF(particles.size());
  unsigned nJet(pfJets.size());

  double ht(0.);
  std::set<unsigned> jetConstituents;
  for(unsigned iJ(0); iJ != nJet; ++iJ){
    susy::PFJet const& jet(pfJets[iJ]);

    TLorentzVector const& p(jet.momentum);
    TLorentzVector corrP(jet.momentum * jet.jecScaleFactors.find("L1FastL2L3")->second);
    if(corrP.Pt() > 30.
       && TMath::Abs(corrP.Eta()) < 3.
       && jet.chargedHadronEnergy / p.E() > 0.
       && jet.neutralHadronEnergy / p.E() < 0.99
       && jet.chargedEmEnergy / p.E() < 0.99
       && jet.neutralEmEnergy / p.E() < 0.99
       && jet.nConstituents > 1
       && jet.chargedMultiplicity > 0)
      ht += corrP.Pt();

    std::vector<unsigned short> const& list(jet.pfParticleList);
    unsigned nC(list.size());
    for(unsigned iC(0); iC != nC; ++iC)
      jetConstituents.insert(list[iC]);

    TEllipse* circle(new TEllipse(jet.momentum.Eta(), jet.momentum.Phi(), TMath::Sqrt(jet.jetArea / TMath::Pi())));
    circle->SetLineColor(kOrange);
    circle->SetLineWidth(2);
    circle->SetLineStyle(kSolid);
    circle->SetFillStyle(0);

    circle->Draw();

    garbageCollection_.Add(circle);
  }

  std::set<std::pair<double, double> > uniqueEtaPhi;

  for(unsigned iP(0); iP != nPF; ++iP){
    susy::PFParticle const& part(particles[iP]);
    double pt(part.momentum.Pt());
    if(pt < ptThreshold_) continue;

    double eta(part.momentum.Eta());
    double phi(part.momentum.Phi());

    std::pair<double, double> etaPhiPair(eta, phi);
    if(uniqueEtaPhi.find(etaPhiPair) != uniqueEtaPhi.end()) continue;
    uniqueEtaPhi.insert(etaPhiPair);

    if(eta < etaPhiField_.GetXaxis()->GetXmin() || eta > etaPhiField_.GetXaxis()->GetXmax()) continue;
    if(phi < etaPhiField_.GetYaxis()->GetXmin() || phi > etaPhiField_.GetYaxis()->GetXmax()) continue;

    TMarker* marker(new TMarker(eta, phi, kFullCircle));
    switch(part.pdgId){
    case 22:
      marker->SetMarkerColor(kGreen);
      break;
    case 11:
    case -11:
      marker->SetMarkerColor(kBlue);
      break;
    case 13:
    case -13:
      marker->SetMarkerColor(kRed);
      break;
    default:
      if(jetConstituents.find(iP) != jetConstituents.end())
        marker->SetMarkerColor(kOrange);
      else
        marker->SetMarkerColor(kBlack);
      break;
    }
    marker->SetMarkerSize(sizeNorm * TMath::Log10(pt));

    marker->Draw();

    garbageCollection_.Add(marker);
  }

  TLine* metLine(new TLine(etaPhiField_.GetXaxis()->GetXmin(), TVector2::Phi_mpi_pi(met.mEt.Phi()), etaPhiField_.GetXaxis()->GetXmax(), TVector2::Phi_mpi_pi(met.mEt.Phi())));
  metLine->SetLineWidth(2);
  metLine->SetLineColor(kRed);
  metLine->SetLineStyle(kDashed);

  metLine->Draw();

  garbageCollection_.Add(metLine);

  eventInfo_.Clear();

  eventInfo_.AddText(TString::Format("MET = %.2f GeV", met.met()));
  eventInfo_.AddText(TString::Format("HT = %.2f GeV", ht));

  eventInfo_.Draw();
}

