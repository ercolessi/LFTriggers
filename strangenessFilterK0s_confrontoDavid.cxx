// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief A filter task for strangeness filter
//  usage:
/*  o2-analysis-timestamp --aod-file /data/dataalice/cdemart/O2/AO2D_15o_t180.root -b | \
    o2-analysis-multiplicity-table -b | \
    o2-analysis-centrality-table -b | \
    o2-analysis-event-selection -b | \
    o2-analysis-trackextension -b | \
    o2-analysis-pid-tpc | \
    o2-analysis-pid-tof | \
    o2-analysis-weak-decay-indices -b | \
    o2-analysis-lambdakzerobuilder --d_bz 5 -b | \
    o2-analysis-strangeness-filter-K0s -b 
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
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/PID/PIDResponse.h"

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

#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct strangenessFilterK0s {

  //Recall the output table
  Produces<aod::StrangenessFiltersK0s> strgtableK0s;

  //Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 3, 0., 3.)};

  //Selection criteria for V0s (David task like)
  Configurable<float> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //is it with respect to Xi decay vertex?
  Configurable<float> dcav0dau{"dcav0dau", 1., "DCA V0 Daughters"};       //is it in sigmas?
  Configurable<float> dcanegtopv{"dcanegtopv", 0.1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.1, "DCA Pos To PV"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<float> v0radius{"v0radius", 5., "V0 Radius"};
  //David does not have these so skip
  /*
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> dcav0topv{"dcav0topv", 0.5, "DCA V0 To PV"};
  Configurable<float> v0radiusupperlimit{"v0radiusupperlimit", 34, "V0 Radius Upper Limit"};
  Configurable<float> LRej{"LRej", 0.005, "LRej"}; 
  Configurable<float> eta{"eta", 0.8, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.8, "EtaDaughters"};
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N Sigmas TPC"};
  Configurable<float> k0smasswindow{"k0smasswindow", 0.03, "K0s Mass Window"};
  */
  //missing selections: OOB pileup?

  //Selections criteria for tracks
  Configurable<float> hEta{"hEta", 0.8f, "Eta range for trigger particles"};
  Configurable<float> hMinPt{"hMinPt", 1.0f, "Min pt for trigger particles"};

  void init(o2::framework::InitContext&)
  {

    //  std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    //  std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    //  AxisSpec centAxis = {centBinning, "V0M (%)"};

    QAHistos.add("VtxZAfterSel", "Vertex distribution in Z;Z (cm)", HistType::kTH1F, {{100, -20, 20}});
    QAHistos.add("hMassK0sBefSel", "hMassK0sBefSel", HistType::kTH1F, {{100, 0.46f, 0.54f}});
    QAHistos.add("hMassK0sAfterSel", "hMassK0sAfterSel", HistType::kTH1F, {{200, 0.450f, 0.550f}});
    QAHistos.add("hMassK0sAfterSelTrig", "hMassK0sAfterSelTrig", HistType::kTH1F, {{100, 0.46f, 0.54f}});
    QAHistos.add("hMassK0sAfterSelEta", "hMassK0sAfterSelEta", HistType::kTH1F, {{100, 0.46f, 0.54f}});
    QAHistos.add("hTriggeredParticles_hasK0S", "Total Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("hTriggeredParticles", "Total Selected triggered particles", HistType::kTH1F, {{10, 0.5, 10.5, "Trigger counter"}});
    QAHistos.add("hHasTriggerParticle", "Events with trigger particles", HistType::kTH1F, {{1, 0.5, 1.5, "Trigger counter"}});
    QAHistos.add("hTriggeredParticles_perevent", "Selected triggered particles per event", HistType::kTH1F, {{10, 0.5, 100.5, "Number of trg particles"}});
    QAHistos.add("PtTrg_hasK0S", "PtTrg_hasK0S", HistType::kTH1F, {{300, 0, 30, "#it{p}_{T} of trg particle"}});
    QAHistos.add("PtTrg", "PtTrg", HistType::kTH1F, {{300, 0, 30, "#it{p}_{T} of trg particle"}});
    QAHistos.add("hDCAdauK0s", "hDCAdauK0s", HistType::kTH1F, {{100, 0.f, 0.54f, "DCA V0 Daughters"}});
    QAHistos.add("hDCAdauK0sSel", "hDCAdauK0sSel", HistType::kTH1F, {{100, 0.f, 0.54f,"DCA V0 Daughters"}});
    QAHistos.add("PtK0s", "PtK0s", HistType::kTH1F, {{300, 0, 30, "#it{p}_{T}"}});
    QAHistos.add("PtK0sSel", "PtK0sSel", HistType::kTH1F, {{300, 0, 30, "#it{p}_{T}"}});
    QAHistos.add("EtaK0s", "EtaK0s", HistType::kTH1F, {{300, -3.0, 3.0, "#it{eta}"}});
    QAHistos.add("EtaK0sSel", "EtaK0sSel", HistType::kTH1F, {{100, -1.0, 1.0, "#it{eta}"}}); 
    QAHistos.add("V0RadiusK0s", "V0RadiusK0s", HistType::kTH1F, {{200, 0.0, 200.0, "V0 Radius"}});
    QAHistos.add("V0RadiusK0sSel", "V0RadiusK0sSel", HistType::kTH1F, {{200, 0.0, 200.0, "V0 Radius"}});
    QAHistos.add("hTOFnsigmaPiPlus", "hTOFnsigmaPiPlus", HistType::kTH1F, {{10, 0, +10, "TOFnsigmaPiPlus"}});
    QAHistos.add("DCANegPVK0s", "DCANegPVK0s", HistType::kTH1F, {{30, 0.0, 30.0, "DCA Neg to PV"}});
    QAHistos.add("DCANegPVK0sSel", "DCANegPVK0sSel", HistType::kTH1F, {{30, 0.0, 30.0, "DCA Neg to PV"}});
    QAHistos.add("DCAPosPVK0s", "DCAPosPVK0s", HistType::kTH1F, {{30, 0.0, 30.0, "DCA Pos to PV"}});
    QAHistos.add("DCAPosPVK0sSel", "DCAPosPVK0sSel", HistType::kTH1F, {{30, 0.0, 30.0, "DCA Pos to PV"}});    
    QAHistos.add("PropLifetimeK0s", "PropLifetimeK0s", HistType::kTH1F, {{200, 0.0, 200.0, "mL/p (cm)"}});    
    QAHistos.add("PropLifetimeK0sSel", "PropLifetimeK0sSel", HistType::kTH1F, {{400, 0.0, 40.0, "mL/p (cm)"}});    
    QAHistos.add("TPCNSigmaPos", "TPCNSigmaPos", HistType::kTH1F, {{10, 0.0, 10.0, "#pi^{+} TPC PID dE/dx (#sigma)"}});    
    QAHistos.add("TPCNSigmaPosSel", "TPCNSigmaPosSel", HistType::kTH1F, {{10, 0.0, 10.0, "#pi^{+} TPC PID dE/dx (#sigma)"}});    
    QAHistos.add("TPCNSigmaNeg", "TPCNSigmaNeg", HistType::kTH1F, {{10, 0.0, 10.0, "#pi^{-} TPC PID dE/dx (#sigma)"}});    
    QAHistos.add("TPCNSigmaNegSel", "TPCNSigmaNegSel", HistType::kTH1F, {{10, 0.0, 10.0, "#pi^{-} TPC PID dE/dx (#sigma)"}});    
    QAHistos.add("V0cosPA", "V0cosPA", HistType::kTH1F, {{100, 0.99, 1., "V0 cos(PA)"}});    
    QAHistos.add("V0cosPASel", "V0cosPASel", HistType::kTH1F, {{100, 0.99, 1., "V0 cos(PA)"}});    
    //TH2
    QAHistos.add("h2MassK0sPtAllSel", "h2MassK0sPtAllSel", HistType::kTH2F, {{65, 0, 6.5,"#it{p}_{T} (GeV/#it{c})"},{100, 0.46f, 0.54f,"M (#pi^{+}#pi^{-})"}});    

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "K0s-K0s");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, Form("#it{p}_{T} hadron - K0s",hMinPt));
  }

  //Filters
  Filter trackFilter = (nabs(aod::track::eta) < hEta) && (aod::track::isGlobalTrack == static_cast<uint8_t>(1u)) && (aod::track::pt>hMinPt); 
 
  //Tables
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Cents>::iterator;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;

  void process(CollisionCandidates const& collision,  TrackCandidates const& tracks, aod::V0Datas const& fullV0, DaughterTracks& dtracks)
 
  {
    
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) { //what's that?
      return;
    }

    hProcessedEvents->Fill(0.5);

    //David does not have this so skip
    //if (TMath::Abs(collision.posZ()) > cutzvertex) return; 
    QAHistos.fill(HIST("VtxZAfterSel"), collision.posZ());

    //Is event good? [0] = DoubleK0s, [1] = high-pT hadron + K0s
    bool keepEvent[2]{false};
    //
    float v0pos = -1.;
    float K0sproperlifetime = -1.;
    float v0ptotmom = -1.;
    int K0scounter = 0;
    const float ctauK0s = 2.6844;   //from PDG
    const float massK0s = 0.497611; //from PDG

    int howmany_highpttracks = 0;
    for (auto track : tracks) { // start loop over tracks
      howmany_highpttracks++;    
      QAHistos.fill(HIST("PtTrg"), track.pt());
      QAHistos.fill(HIST("hTriggeredParticles"), 1);  
    }
    //
    QAHistos.fill(HIST("hTriggeredParticles_perevent"), howmany_highpttracks); 
    if (howmany_highpttracks>0) QAHistos.fill(HIST("hHasTriggerParticle"), 1); 
    
    for (auto& v0 : fullV0) { //loop over V0s
     
      //Position
      v0pos = TMath::Sqrt(TMath::Power(v0.x() - collision.posX(), 2) + TMath::Power(v0.y() - collision.posY(), 2) + TMath::Power(v0.z() - collision.posZ(), 2));
      //Total momentum
      v0ptotmom = TMath::Sqrt(v0.px() * v0.px() + v0.py() * v0.py() + v0.pz() * v0.pz());
      //Proper lifetime
      K0sproperlifetime = massK0s * v0pos / (v0ptotmom + 1e-13);

      //Bef selections
      QAHistos.fill(HIST("hMassK0sBefSel"), v0.mK0Short());
      QAHistos.fill(HIST("hDCAdauK0s"), v0.dcaV0daughters());
      QAHistos.fill(HIST("PtK0s"), v0.pt());
      QAHistos.fill(HIST("EtaK0s"), v0.eta());
      QAHistos.fill(HIST("V0RadiusK0s"), v0.v0radius());
      QAHistos.fill(HIST("DCAPosPVK0s"), v0.dcapostopv());
      QAHistos.fill(HIST("DCANegPVK0s"), v0.dcanegtopv());
      QAHistos.fill(HIST("PropLifetimeK0s"), K0sproperlifetime);
      QAHistos.fill(HIST("TPCNSigmaNeg"), v0.negTrack_as<DaughterTracks>().tpcNSigmaPi());
      QAHistos.fill(HIST("TPCNSigmaPos"), v0.posTrack_as<DaughterTracks>().tpcNSigmaPi());
      QAHistos.fill(HIST("V0cosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      QAHistos.fill(HIST("hTOFnsigmaPiPlus"), TMath::Abs(v0.posTrack_as<DaughterTracks>().tofNSigmaPi()));
      if (TMath::Abs(v0.eta()) < 0.8 && TMath::Abs(v0.posTrack_as<DaughterTracks>().eta()) < 0.8 && TMath::Abs(v0.negTrack_as<DaughterTracks>().eta()) < 0.8 )
      {
        QAHistos.fill(HIST("hMassK0sAfterSelEta"), v0.mK0Short());
      }

      //-----------------------                                                                                                                  
      // TOPOLOGICAL - KINEMATIC SELECTIONS                                                                                                      
      //-----------------------   

      //SELECTION DAVID LIKE
      if (TMath::Abs(v0.yK0Short()) > rapidity) continue;
      if (v0.v0radius() < v0radius) continue;
      if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) continue;
      if (TMath::Abs(v0.dcapostopv()) < dcapostopv) continue;
      if (TMath::Abs(v0.dcanegtopv()) < dcanegtopv) continue;
      if (v0.dcaV0daughters() > dcav0dau) continue;
            
      //David does not have these so skip
      /*                                                                               
      if (TMath::Abs(v0.eta()) > 0.8) continue;
      if (TMath::Abs(v0.posTrack_as<DaughterTracks>().eta()) > 0.8) continue;
      if (TMath::Abs(v0.negTrack_as<DaughterTracks>().eta()) > 0.8) continue;
      if (K0sproperlifetime > 30.) continue;
      if (TMath::Abs(v0.mK0Short() - massK0s) > k0smasswindow) continue;    
      if (v0.pt() < 0.1) continue; 
      if (TMath::Abs(v0.posTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue;
      if (TMath::Abs(v0.negTrack_as<DaughterTracks>().tpcNSigmaPi()) > nsigmatpc) continue;
      */

      //Fill QA 
      QAHistos.fill(HIST("hMassK0sAfterSel"), v0.mK0Short());
      QAHistos.fill(HIST("PtK0sSel"), v0.pt());
      QAHistos.fill(HIST("hDCAdauK0sSel"), v0.dcaV0daughters());
      QAHistos.fill(HIST("EtaK0sSel"), v0.eta());
      QAHistos.fill(HIST("V0RadiusK0sSel"), v0.v0radius());
      QAHistos.fill(HIST("DCAPosPVK0sSel"), v0.dcapostopv());
      QAHistos.fill(HIST("DCANegPVK0sSel"), v0.dcanegtopv());
      QAHistos.fill(HIST("PropLifetimeK0sSel"), K0sproperlifetime);
      QAHistos.fill(HIST("TPCNSigmaNegSel"), v0.negTrack_as<DaughterTracks>().tpcNSigmaPi());
      QAHistos.fill(HIST("TPCNSigmaPosSel"), v0.posTrack_as<DaughterTracks>().tpcNSigmaPi());
      QAHistos.fill(HIST("V0cosPASel"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));


      //Count number of K0s candidates
      K0scounter++;

      if (howmany_highpttracks>0) {
        QAHistos.fill(HIST("hMassK0sAfterSelTrig"), v0.mK0Short());
        QAHistos.fill(HIST("h2MassK0sPtAllSel"),v0.pt(), v0.mK0Short());
      }

    } //end loop over V0s
    
    //Double K0s trigger definition
    if (K0scounter > 1) {
      keepEvent[0] = true;
    }
    
    //High-pT hadron + K0s trigger definition
    if (K0scounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        QAHistos.fill(HIST("hTriggeredParticles_hasK0S"), 1);
        QAHistos.fill(HIST("PtTrg_hasK0S"), track.pt());

        keepEvent[1] = true;       
        break;
      } // end loop over tracks
    }
    
    if (keepEvent[0])
      hProcessedEvents->Fill(1.5);
    if (keepEvent[1])
      hProcessedEvents->Fill(2.5);
    
    //Filling the table
    strgtableK0s(keepEvent[0],keepEvent[1]);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilterK0s>(cfgc, TaskName{"strangeness-filter-K0s"})};
}
