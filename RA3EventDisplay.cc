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
#include "TVector2.h"

#include <iostream>
#include <map>
#include <utility>
#include <set>
#include <stdexcept>

double const sizeNorm(2.);
double const etaMax(4.);
TString const HtDescription("H_{T}: P_{T}^{j} > 30 GeV && |#eta^{j}| < 3");

class RA3EventDisplay {
public:
  RA3EventDisplay();
  ~RA3EventDisplay();

  int addPath(TString const&);
  void setPtThreshold(double);
  void setPFJetCollection(TString const&);
  void addMet(TString const&, bool = true, int = 0);
  bool setMainMet(TString const&);

  bool showEvent(unsigned, unsigned);
  bool showNextEvent();
  bool refresh();

  double mass(unsigned, unsigned = -1, unsigned = -1, unsigned = -1) const;
  double mt(unsigned, unsigned = -1, unsigned = -1, unsigned = -1) const;
  TString const& getMainMet() const { return mainMet_; }

  void print(TString const& _fileName) { canvas_.Print(_fileName); }

private:
  void showEvent_();
  void formLegend_(bool = false);

  TChain mainTree_;
  TChain scanTree_;

  double ptThreshold_;
  TString pfJetCollectionTag_;
  std::map<TString, int> metTags_;
  TString mainMet_;

  unsigned runNumber_;
  unsigned luminosityBlockNumber_;
  unsigned eventNumber_;
  long currentEntry_;

  susy::Event event_;

  TCanvas canvas_;
  TH2F etaPhiField_;
  TPaveText eventInfo_;
  TPaveText notes_;

  TLegend sizeLegend_;
  TLegend colorLegend_;

  TObjArray garbageCollection_;
};

RA3EventDisplay::RA3EventDisplay() :
  mainTree_("susyTree"),
  scanTree_("susyTree"),
  ptThreshold_(1.),
  pfJetCollectionTag_("ak5"),
  metTags_(),
  mainMet_("pfType01CorrectedMet"),
  runNumber_(0),
  luminosityBlockNumber_(0),
  eventNumber_(0),
  currentEntry_(-1),
  event_(),
  canvas_("evdisp", "RA3 Event Display", 800, 600),
  etaPhiField_("evdisp", ";#eta;#phi", 100, -etaMax, etaMax, 100, -TMath::Pi(), TMath::Pi()),
  eventInfo_(0.82, 0.8, 0.98, 0.95, "brNDC"),
  notes_(0.82, 0.7, 0.98, 0.75, "brNDC"),
  sizeLegend_(0.82, 0.4, 0.98, 0.55),
  colorLegend_(0.82, 0.05, 0.98, 0.35),
  garbageCollection_()
{
  mainTree_.SetBranchStatus("*", 0);
  mainTree_.SetBranchStatus("pfParticles*", 1);
  mainTree_.SetBranchStatus("pfJets*", 1);
  mainTree_.SetBranchStatus("met*", 1);

  scanTree_.SetBranchStatus("*", 0);
  scanTree_.SetBranchStatus("runNumber", 1);
  scanTree_.SetBranchStatus("luminosityBlockNumber", 1);
  scanTree_.SetBranchStatus("eventNumber", 1);
  scanTree_.SetBranchAddress("runNumber", &runNumber_);
  scanTree_.SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber_);
  scanTree_.SetBranchAddress("eventNumber", &eventNumber_);

  metTags_[mainMet_] = kRed;

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

  eventInfo_.SetBorderSize(1);

  notes_.SetBorderSize(1);
  notes_.SetTextAlign(12); // left, center
  notes_.SetFillStyle(0);
  notes_.AddText(HtDescription);

  garbageCollection_.SetOwner(true);

  formLegend_();
}

RA3EventDisplay::~RA3EventDisplay()
{
  event_.releaseTree(mainTree_);

  formLegend_(true);
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
RA3EventDisplay::setPtThreshold(double _ptThresh)
{
  if(_ptThresh < 1.){
    std::cerr << "Cannot set threshold to lower than 1 GeV" << std::endl;
    return;
  }
  ptThreshold_ = _ptThresh;
}

void
RA3EventDisplay::setPFJetCollection(TString const& _pfJetCollectionTag)
{
  pfJetCollectionTag_ = _pfJetCollectionTag;

  formLegend_();

  refresh();
}

void
RA3EventDisplay::addMet(TString const& _metName, bool _add/* = true*/, int _color/* = 0*/)
{
  if(_add){
    int color(_color);
    if(color == 0){
      color = kRed;
      while(true){
        std::map<TString, int>::iterator itr(metTags_.begin());
        for(; itr != metTags_.end(); ++itr)
          if(itr->second == color) break;

        if(itr != metTags_.end()) ++color;
        else break;
      }
    }

    metTags_[_metName] = color;

    colorLegend_.SetY2(colorLegend_.GetY2() + 0.05);
  }
  else{
    std::map<TString, int>::iterator itr(metTags_.find(_metName));
    if(itr == metTags_.end()) return;
    metTags_.erase(itr);
    colorLegend_.SetY2(colorLegend_.GetY2() - 0.05);

    if(mainMet_ == _metName){
      if(metTags_.size() != 0)
        mainMet_ = metTags_.begin()->first;
      else
        mainMet_ = "";
    }
  }

  formLegend_();

  setMainMet(mainMet_);

  refresh();
}

bool
RA3EventDisplay::setMainMet(TString const& _metTag)
{
  if(metTags_.find(_metTag) == metTags_.end()) return false;

  mainMet_ = _metTag;

  notes_.Clear();
  notes_.AddText(HtDescription);
  if(metTags_.size() != 1) notes_.AddText("MET: " + mainMet_);

  formLegend_();

  refresh();

  return true;
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
      currentEntry_ = -1;
      if(start == -1)
        break;
    }
    else if(bytes < 0)
      throw std::runtime_error("Input error");
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
  int bytes(scanTree_.GetEntry(++currentEntry_));
  if(bytes > 0){
    showEvent_();
    return true;
  }
  else if(bytes == 0){
    currentEntry_ = -1;
    std::cerr << "End of input" << std::endl;
    return false;
  }
  else
    throw std::runtime_error("Input error");
    
}

bool
RA3EventDisplay::refresh()
{
  if(currentEntry_ != -1){
    showEvent_();
    return true;
  }
  else
    return false;
}

double
RA3EventDisplay::mass(unsigned _p0, unsigned _p1/* = -1*/, unsigned _p2/* = -1*/, unsigned _p3/* = -1*/) const
{
  if(currentEntry_ == -1) return 0.;

  susy::PFParticleCollection const& particles(event_.pfParticles);

  TLorentzVector sumP;
  if(_p0 < particles.size()) sumP += particles[_p0].momentum;
  if(_p1 < particles.size()) sumP += particles[_p1].momentum;
  if(_p2 < particles.size()) sumP += particles[_p2].momentum;
  if(_p3 < particles.size()) sumP += particles[_p3].momentum;

  return sumP.M();
}

double
RA3EventDisplay::mt(unsigned _p0, unsigned _p1/* = -1*/, unsigned _p2/* = -1*/, unsigned _p3/* = -1*/) const
{
  if(currentEntry_ == -1) return 0.;

  susy::PFParticleCollection const& particles(event_.pfParticles);
  susy::METMap::const_iterator mItr(event_.metMap.find(mainMet_));
  if(mItr == event_.metMap.end()){
    std::cerr << mainMet_ << " not found" << std::endl;
    return 0.;
  }

  TVector2 const& met(mItr->second.mEt);

  TVector2 sumPt;
  if(_p0 < particles.size()) sumPt += TVector2(particles[_p0].momentum.X(), particles[_p0].momentum.Y());
  if(_p1 < particles.size()) sumPt += TVector2(particles[_p1].momentum.X(), particles[_p1].momentum.Y());
  if(_p2 < particles.size()) sumPt += TVector2(particles[_p2].momentum.X(), particles[_p2].momentum.Y());
  if(_p3 < particles.size()) sumPt += TVector2(particles[_p3].momentum.X(), particles[_p3].momentum.Y());

  return TMath::Sqrt(TMath::Power(sumPt.Mod() + met.Mod(), 2.) - (sumPt + met).Mod2());
}

void
RA3EventDisplay::showEvent_()
{
  garbageCollection_.Clear();
  canvas_.Clear();

  etaPhiField_.SetTitle(TString::Format("Run %d, Lumi %d, Event %d", runNumber_, luminosityBlockNumber_, eventNumber_));
  etaPhiField_.Draw();

  colorLegend_.Draw();
  sizeLegend_.Draw();
  notes_.Draw();

  event_.getEntry(currentEntry_);

  susy::PFParticleCollection const& particles(event_.pfParticles);
  susy::PFJetCollectionMap::const_iterator jetColItr(event_.pfJets.find(pfJetCollectionTag_));
  if(jetColItr == event_.pfJets.end()){
    std::cerr << pfJetCollectionTag_ << " not found" << std::endl;
    canvas_.Clear();
    return;
  }

  susy::PFJetCollection const& pfJets(jetColItr->second);

  unsigned nPF(particles.size());
  unsigned nJet(pfJets.size());

  double ht(0.);
  std::set<unsigned> jetConstituents;
  for(unsigned iJ(0); iJ != nJet; ++iJ){
    susy::PFJet const& jet(pfJets[iJ]);

    TLorentzVector const& p(jet.momentum);
    TLorentzVector corrP(jet.momentum * jet.jecScaleFactors.find("L1FastL2L3")->second);
    if(corrP.Pt() < 1.) continue; // never happens, but just for mathematical sanity..

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

    // Marker size 1 is absolute 8 pix. Ellipse radius is given in dR, so the conversion factor is:
    // radius = 8 * log10(jetEt) * sizeNorm * 2etaMax / (windowWidth * 0.75) 
    TEllipse* circle(new TEllipse(jet.momentum.Eta(), jet.momentum.Phi(), 8. * TMath::Log10(corrP.Pt()) * sizeNorm * 2. * etaMax / canvas_.GetWindowWidth() / 0.75));
    circle->SetLineColor(kOrange);
    circle->SetLineWidth(2);
    circle->SetLineStyle(kSolid);
    circle->SetFillStyle(0);

    circle->Draw();

    garbageCollection_.Add(circle);
  }

  std::map<std::pair<double, double>, unsigned> sortedAndCleaned;
  for(unsigned iP(0); iP != nPF; ++iP){
    susy::PFParticle const& part(particles[iP]);
    double pt(part.momentum.Pt());
    if(pt < ptThreshold_) continue;
    double eta(part.momentum.Eta());
    if(eta < -etaMax || eta > etaMax) continue;

    std::pair<double, double> key(pt, eta);
    std::map<std::pair<double, double>, unsigned>::iterator itr(sortedAndCleaned.find(key));
    if(itr != sortedAndCleaned.end()){
      if(jetConstituents.find(iP) != jetConstituents.end()) itr->second = iP;
    }
    else
      sortedAndCleaned.insert(std::pair<std::pair<double, double>, unsigned>(key, iP));
  }

  std::map<std::pair<double, double>, unsigned>::reverse_iterator pEnd(sortedAndCleaned.rend());
  for(std::map<std::pair<double, double>, unsigned>::reverse_iterator pItr(sortedAndCleaned.rbegin()); pItr != pEnd; ++pItr){
    unsigned iP(pItr->second);
    susy::PFParticle const& part(particles[iP]);

    double eta(part.momentum.Eta());
    double phi(part.momentum.Phi());

    TMarker* marker(new TMarker(eta, phi, kFullCircle));
    marker->SetUniqueID(iP);
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
    marker->SetMarkerSize(sizeNorm * TMath::Log10(part.momentum.Pt()));

    marker->Draw();

    garbageCollection_.Add(marker);
  }

  double metVal(0.);

  for(std::map<TString, int>::iterator tItr(metTags_.begin()); tItr != metTags_.end(); ++tItr){
    susy::METMap::const_iterator metItr(event_.metMap.find(tItr->first));
    if(metItr == event_.metMap.end()){
      std::cerr << tItr->first << " not found" << std::endl;
      canvas_.Clear();
      return;
    }

    susy::MET const& met(metItr->second);

    TLine* metLine(new TLine(-etaMax, TVector2::Phi_mpi_pi(met.mEt.Phi()), etaMax, TVector2::Phi_mpi_pi(met.mEt.Phi())));
    metLine->SetLineWidth(2);
    metLine->SetLineColor(tItr->second);
    metLine->SetLineStyle(kDashed);

    metLine->Draw();

    garbageCollection_.Add(metLine);

    if(tItr->first == mainMet_) metVal = met.met();
  }

  eventInfo_.Clear();

  eventInfo_.AddText(TString::Format("MET = %.2f GeV", metVal));
  eventInfo_.AddText(TString::Format("HT = %.2f GeV", ht));

  eventInfo_.Draw();
}

void
RA3EventDisplay::formLegend_(bool _clearOnly/* = false*/)
{
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

  colorLegend_.Clear();
  sizeLegend_.Clear();

  if(_clearOnly) return;

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
  TMarker* jet(new TMarker(0., 0., kOpenCircle));
  photon->SetMarkerColor(kGreen);
  electron->SetMarkerColor(kBlue);
  muon->SetMarkerColor(kRed);
  jetHadron->SetMarkerColor(kOrange);
  strayHadron->SetMarkerColor(kBlack);
  jet->SetMarkerColor(kOrange);

  colorLegend_.AddEntry(photon, "Photon", "P");
  colorLegend_.AddEntry(electron, "Electron", "P");
  colorLegend_.AddEntry(muon, "Muon", "P");
  colorLegend_.AddEntry(jetHadron, "Hadron in jet", "P");
  colorLegend_.AddEntry(strayHadron, "Hadron not in jet", "P");
  colorLegend_.AddEntry(jet, pfJetCollectionTag_ + "PFJet", "P");

  for(std::map<TString, int>::iterator mItr(metTags_.begin()); mItr != metTags_.end(); ++mItr){
    TLine* metLine(new TLine(-etaMax, 0., etaMax, 0.));
    metLine->SetLineWidth(2);
    metLine->SetLineColor(mItr->second);
    metLine->SetLineStyle(kDashed);

    colorLegend_.AddEntry(metLine, mItr->first, "L");
  }
}
