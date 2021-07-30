// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Example V0 tutorial task 1: a simple looper
//
// Step 1: add extra histogram for lambda and antilambda mass
// Step 2: add extra 3D histogram with masses, momenta and centralities
// Step 3: add configurable TPC dE/dx selection
// Step 4: vary selections
// Step 5: try with finder instead of builder
//
/// \brief A filter task for strangeness filter
//  usage:
/*  o2-analysis-timestamp --aod-file ~/Scaricati/AO2D_PbPb.root -b | \
    o2-analysis-multiplicity-table -b | \
    o2-analysis-centrality-table -b | \
    o2-analysis-event-selection -b | \
    o2-analysis-trackextension -b | \
    o2-analysis-pid-tpc | \
    o2-analysis-pid-tof | \
    o2-analysis-weak-decay-indices -b | \
    o2-analysis-lambdakzerobuilder --d_bz 5 -b | \
    o2-analysis-cascadebuilder --d_bz 5 -b | \
    o2-analysis-strangeness-filter -b 
*/
///
///
/// \author Chiara De Martin (chiara.de.martin@cern.ch) and Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since June 1, 2021

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "AnalysisCore/RecoDecay.h"
#include "AnalysisCore/trackUtilities.h"
#include "AnalysisDataModel/StrangenessTables.h"
#include "AnalysisCore/TrackSelection.h"
#include "AnalysisDataModel/TrackSelectionTables.h"
#include "AnalysisDataModel/EventSelection.h"
#include "AnalysisDataModel/Centrality.h"
#include "AnalysisDataModel/PID/PIDResponse.h"
#include "filterTables.h"

#include <TFile.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <Math/Vector4D.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct strangenessFilter {

  //Recall the output table
  Produces<aod::StrangenessFilters> strgtable;

  //Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry EventsvsMultiplicity{"EventsvsMultiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 4, 0., 4.)};

  //Selection criteria for cascades
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> v0cospa{"v0cospa", 0.97, "V0 CosPA"}; //is it with respect to Xi decay vertex?
  Configurable<float> casccospa{"casccospa", 0.995, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};       //is it in sigmas?
  Configurable<float> dcacascdau{"dcacascdau", 0.8, "DCA Casc Daughters"}; //is it in sigmas?
  Configurable<float> dcamesontopv{"dcamesontopv", 0.04, "DCA Meson To PV"};
  Configurable<float> dcabaryontopv{"dcabaryontopv", 0.03, "DCA Baryon To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.04, "DCA Bach To PV"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.02, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.02, "DCA Pos To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 1.2, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.2, "V0 Radius"};
  Configurable<float> v0radiusupperlimit{"v0radiusupperlimit", 34, "V0 Radius Upper Limit"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> cascradiusupperlimit{"cascradiusupperlimit", 34, "Casc Radius Upper Limit"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"};
  Configurable<float> masslambdalimit{"masslambdalimit", 0.01, "masslambdalimit"}; //0.006 Chiara
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"}; 
  Configurable<float> xirej{"xirej", 0.008, "xirej"};  //merge the two rejection variables into one?
  Configurable<float> ximasswindow{"ximasswindow", 0.075, "Xi Mass Window"};
  Configurable<float> omegamasswindow{"omegamasswindow", 0.075, "Omega Mass Window"}; //merge the two windows variables into one?   
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};
  //missing selections: OOB pileup?
  //eta and y selections: loose enough?
  //eta selections of daughters

  //Selections criteria for tacks                                                                                                                
  Configurable<float> hEta{"hEta", 0.8f, "Eta range for trigger particles"};
  Configurable<float> hMinPt{"hMinPt", 1.0f, "Min pt for trigger particles"};

  void init(o2::framework::InitContext&)
  {

    //  std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec ximassAxis = {100, 1.30f, 1.34f};
    AxisSpec omegamassAxis = {100, 1.5f, 1.8f};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}; 

    QAHistos.add("Centrality", "Centrality distribution (V0M)", HistType::kTH1F, {{100, 0, 100, "V0M (%)"}});
    QAHistos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", HistType::kTH1F, {{100, -20, 20}});
    QAHistos.add("hMassXiBefSel", "hMassXiBefSel", HistType::kTH1F, {ximassAxis});
    QAHistos.add("hMassXiAfterSel", "hMassXiAfterSel", HistType::kTH1F, {ximassAxis});
    QAHistos.add("hMassOmegaBefSel", "hMassOmegaBefSel", HistType::kTH1F, {omegamassAxis});
    QAHistos.add("hMassOmegaAfterSel", "hMassOmegaAfterSel", HistType::kTH1F, {omegamassAxis});
    QAHistos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaAfterSelvsPt", "hMassOmegaAfterSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hTriggeredParticles", "Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("PtTrigger", "PtTrigger", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particle"}});

    EventsvsMultiplicity.add("AllEventsvsMultiplicity", "Multiplicity distribution of all events", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("DXiEventsvsMultiplicity", "Multiplicity distribution of events with >1Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("hXiEventsvsMultiplicity", "Multiplicity distribution of events with h + Xi", HistType::kTH1F, {centAxis});
    EventsvsMultiplicity.add("OmegaEventsvsMultiplicity", "Multiplicity distribution of events with >=1 Omega", HistType::kTH1F, {centAxis});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "#Xi-#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "high-#it{p}_{T} hadron - #Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "#Omega");
  }

  //Filters
  Filter collisionFilter = (nabs(aod::collision::posZ) < cutzvertex);
  Filter trackFilter = (nabs(aod::track::eta) < hEta) && (aod::track::pt>hMinPt) && (aod::track::isGlobalTrack == static_cast<uint8_t>(1u)) ;   
  Filter preFilterCasc = nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& aod::cascdata::dcaV0daughters < dcav0dau&& aod::cascdata::dcacascdaughters < dcacascdau;

  //Tables
  using CollisionCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Cents>>::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;
  using DaughterTracks = soa::Join<aod::FullTracks, aod::TracksExtended, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;

  void process(CollisionCandidates const& collision, TrackCandidates const& tracks, aod::V0Datas const& V0s,soa::Join<aod::Cascades, aod::CascDataExt> const& fullCasc, DaughterTracks& dtracks)
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) { //what's that?
      return;
    }

    QAHistos.fill(HIST("VtxZAfterSel"), collision.posZ());
    QAHistos.fill(HIST("Centrality"), collision.centV0M());
    EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicity"), collision.centV0M());
    hProcessedEvents->Fill(0.5);

    //Is event good? [0] = DoubleXi, [1] = high-pT hadron + Xi
    bool keepEvent[3]{false};
    //
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int omegacounter = 0;
    const float ctauxi = 4.91;    //from PDG
    const float massxi = 1.32171; //from PDG
    const float ctauomega = 2.461;    //from PDG
    const float massomega = 1.67245; //from PDG

    for (auto& casc : fullCasc) { //loop over cascades

      bool isXi = false;
      bool isOmega = false;

      QAHistos.fill(HIST("hMassXiBefSel"), casc.mXi());
      QAHistos.fill(HIST("hMassOmegaBefSel"), casc.mOmega());

      //Position
      xipos = TMath::Sqrt(TMath::Power(casc.x() - collision.posX(), 2) + TMath::Power(casc.y() - collision.posY(), 2) + TMath::Power(casc.z() - collision.posZ(), 2));
      //Total momentum
      xiptotmom = TMath::Sqrt(casc.px() * casc.px() + casc.py() * casc.py() + casc.pz() * casc.pz());
      //Proper lifetime
      xiproperlifetime = massxi * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = massomega * xipos / (xiptotmom + 1e-13);

      //Selections
      auto v0 = casc.v0_as<aod::V0Datas>();

      if (casc.sign() == 1) {
         if (TMath::Abs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
          continue;
        if (TMath::Abs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc)
          continue;
      } else {
        if (TMath::Abs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPr()) > nsigmatpc)
          continue;
        if (TMath::Abs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
          continue;
      }
      //this selection differes for Xi and Omegas:
      // if (TMath::Abs(casc.bachelor_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc)
      //  continue; 

      //if (TMath::Abs(casc.bachelor_as<DaughterTracks>().eta()) > etadau)
      //  continue; //?
      if (TMath::Abs(v0.posTrack_as<DaughterTracks>().eta()) > etadau)
        continue; 
      if (TMath::Abs(v0.negTrack_as<DaughterTracks>().eta()) > etadau)
        continue; 
      
      if (casc.sign() == 1) {
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv)
          continue;
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv)
          continue;
      } else {
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv)
          continue;
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv)
          continue;
      }
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv)
        continue;
      if (casc.v0radius() > v0radiusupperlimit || casc.v0radius() < v0radius)
        continue;
      if (casc.cascradius() > cascradiusupperlimit || casc.cascradius() < cascradius)
        continue;
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa)
        continue; //is it calculated with respect to the Xi decay vertex?
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < casccospa)
        continue;
      if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < dcav0topv)
        continue;
      if (casc.dcaV0daughters() > dcav0dau)
        continue;
      if (casc.dcacascdaughters() > dcacascdau)
        continue;
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit)
        continue;
      if (TMath::Abs(casc.eta()) > eta)
        continue;

      isXi = (TMath::Abs(casc.mXi() - massxi) < ximasswindow) && (TMath::Abs(casc.mOmega() - massomega)>omegarej) && (xiproperlifetime < properlifetimefactor * ctauxi) && (TMath::Abs(casc.yXi()) < rapidity); //add PID on bachelor
      isOmega = (TMath::Abs(casc.mOmega() - massomega) < omegamasswindow) && (TMath::Abs(casc.mXi() - massxi)>xirej) && (omegaproperlifetime < properlifetimefactor * ctauomega) && (TMath::Abs(casc.yOmega()) < rapidity); //add PID on bachelor
      if (isXi){
	QAHistos.fill(HIST("hMassXiAfterSel"), casc.mXi());
	QAHistos.fill(HIST("hMassXiAfterSel2D"), casc.mXi(), casc.pt());
	//Count number of Xi candidates
	xicounter++;
      }
      if (isOmega){
	QAHistos.fill(HIST("hMassOmegaAfterSel"), casc.mOmega());
	QAHistos.fill(HIST("hMassOmegaAfterSel2D"), casc.mOmega(), casc.pt());
	//Count number of Omega candidates
	omegacounter++;
      }
    } //end loop over cascades

    //Double Xi trigger definition
    if (xicounter > 1) {
      keepEvent[0] = true;
    }

    //High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
      for (auto track : tracks) { // start loop over tracks
	//all needed selections applied via aod::track::isGlobalTrack == static_cast<uint8_t>(1u)
	//no we need track length selections?

        QAHistos.fill(HIST("hTriggeredParticles"), 1);
        QAHistos.fill(HIST("PtTrigger"), track.pt());

        keepEvent[1] = true;
      } // end loop over tracks
    }

    //omega trigger definition
    if (omegacounter > 0) {
      keepEvent[2] = true;
    }

    if (keepEvent[0]){
      hProcessedEvents->Fill(1.5);
      EventsvsMultiplicity.fill(HIST("DXiEventsvsMultiplicity"), collision.centV0M());
    }
    if (keepEvent[1]){
      hProcessedEvents->Fill(2.5);
      EventsvsMultiplicity.fill(HIST("hXiEventsvsMultiplicity"), collision.centV0M());
    }
    if (keepEvent[2]){
      hProcessedEvents->Fill(3.5);
      EventsvsMultiplicity.fill(HIST("OmegaEventsvsMultiplicity"), collision.centV0M());
    }

    //Filling the table
    strgtable(keepEvent[0],keepEvent[1], keepEvent[2]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"strangeness-filter"})};
}
