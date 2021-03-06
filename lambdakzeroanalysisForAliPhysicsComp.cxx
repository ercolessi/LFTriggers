// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Example V0 analysis task
// ========================
//
// This code loops over a V0Data table and produces some
// standard analysis output. It requires either
// the lambdakzerofinder or the lambdakzeroproducer tasks
// to have been executed in the workflow (before).
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//

//o2-analysis-timestamp --aod-file AO2D_17p_pass1_FAST_AOD234_001To5.root -b |
//o2-analysis-event-selection -b |
//o2-analysis-multiplicity-table -b |
//o2-analysis-centrality-table -b |
//o2-analysis-pid-tof -b |
//o2-analysis-pid-tpc -b |
//o2-analysis-trackextension -b |
//o2-analysis-trackselection -b |
//o2-analysis-trackqa -b |
//o2-analysis-weak-decay-indices -b |
//o2-analysis-lf-lambdakzerobuilder --d_bz 5 --v0cospa 0.95 --dcav0dau 2 --v0radius 0.9 --dcanegtopv 0.03 --dcapostopv 0.03 -b |
//o2-analysis-lf-lambdakzeroanalysis -b

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/StrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
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

struct lambdakzeroQA {
  //Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{100, 0.46f, 0.54f}}}},
      {"hMassLambda", "hMassLambda", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}}},
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{1000, 0.0f, 100.0f}}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{5000, 0.5f, 1.0f}}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{1000, 0.0f, 10.0f}}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {{2000, 0.0f, 20.0f}}}},
    },
  };

  OutputObj<TH1F> henum{TH1F("henum", "Event counter", 4, 0., 4.)};
  OutputObj<TH1F> hvertexZ{TH1F("hvertexZ", "Z vertex", 400, -20., 20.)};
  OutputObj<TH1F> hMult{TH1F("hMult", "Multiplicity distribution", 100, 0., 100.)};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};

  void init(InitContext const&)
  {
    henum->GetXaxis()->SetBinLabel(1, "All events");
    henum->GetXaxis()->SetBinLabel(2, "kINT7");
    henum->GetXaxis()->SetBinLabel(3, "|Zvtx|<10 cm");
    henum->GetXaxis()->SetBinLabel(4, "sel7");
  }

  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator;
  
  void process(CollisionCandidates const& collision, aod::V0Datas const& fullV0s)
  {

    henum->Fill(0.5);
    if (!collision.alias()[kINT7]) {
      return;
    }
    henum->Fill(1.5);
    if (TMath::Abs(collision.posZ()) >= cutzvertex){
      return;
    }
    henum->Fill(2.5);
    if (!collision.sel7()) {
      //return; 
    }
    //    henum->Fill(3.5);
    hvertexZ->Fill(collision.posZ());
    hMult->Fill(collision.centV0M());

    for (auto& v0 : fullV0s) {
      registry.fill(HIST("hMassK0Short"), v0.mK0Short());
      registry.fill(HIST("hMassLambda"), v0.mLambda());
      registry.fill(HIST("hMassAntiLambda"), v0.mAntiLambda());

      registry.fill(HIST("hV0Radius"), v0.v0radius());
      registry.fill(HIST("hV0CosPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), v0.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), v0.dcanegtopv());
      registry.fill(HIST("hDCAV0Dau"), v0.dcaV0daughters());
    }
  }
};

struct lambdakzeroanalysis {

  HistogramRegistry registry{
    "registry",
    {
      {"henumBis", "henumBis", {HistType::kTH1F, {{3,0.0f,3.0f}}}},
      {"h3dMassK0Short", "h3dMassK0Short", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 0.450f, 0.550f}}}},
      {"h3dMassLambda", "h3dMassLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
      {"h3dMassAntiLambda", "h3dMassAntiLambda", {HistType::kTH3F, {{20, 0.0f, 100.0f}, {200, 0.0f, 10.0f}, {200, 1.015f, 1.215f}}}},
    },
  };

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisK0Short = {200, 0.450f, 0.550f, "Inv. Mass (GeV)"};
    AxisSpec massAxisLambda = {200, 1.015f, 1.215f, "Inv. Mass (GeV)"};

    registry.add("h3dMassK0ShortDca", "h3dMassK0ShortDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisK0Short}});
    registry.add("h3dMassLambdaDca", "h3dMassLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
    registry.add("h3dMassAntiLambdaDca", "h3dMassAntiLambdaDca", {HistType::kTH3F, {dcaAxis, ptAxis, massAxisLambda}});
  }

  //Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.995, "V0 CosPA"}; //double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", .1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", .1, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 5.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "rapidity"};
  Configurable<int> saveDcaHist{"saveDcaHist", 0, "saveDcaHist"};

  static constexpr float defaultLifetimeCuts[1][2] = {{25., 20.}};
  Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {defaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0S"}}, "lifetimecut"};

  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range"};
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;


  void process(soa::Join<aod::Collisions, aod::EvSels, aod::CentV0Ms>::iterator const& collision, soa::Filtered<aod::V0Datas> const& fullV0s)
    {
    registry.fill(HIST("henumBis"),0.5);
    if (TMath::Abs(collision.posZ()) >= cutzvertex){
      return;
    }
    registry.fill(HIST("henumBis"),1.5);
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      //      return;
    }
    registry.fill(HIST("henumBis"),2.5);
    for (auto& v0 : fullV0s) {
      //FIXME: could not find out how to filter cosPA and radius variables (dynamic columns)
      if (v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa) {
        if (TMath::Abs(v0.yLambda()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kLambda0) < lifetimecut->get("lifetimecutLambda")) {
            registry.fill(HIST("h3dMassLambda"), collision.centV0M(), v0.pt(), v0.mLambda());
            registry.fill(HIST("h3dMassAntiLambda"), collision.centV0M(), v0.pt(), v0.mAntiLambda());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mLambda());
              registry.fill(HIST("h3dMassAntiLambdaDca"), v0.dcaV0daughters(), v0.pt(), v0.mAntiLambda());
            }
          }
        }
        if (TMath::Abs(v0.yK0Short()) < rapidity) {
          if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(kK0Short) < lifetimecut->get("lifetimecutK0S")) {
            registry.fill(HIST("h3dMassK0Short"), collision.centV0M(), v0.pt(), v0.mK0Short());
            if (saveDcaHist == 1) {
              registry.fill(HIST("h3dMassK0ShortDca"), v0.dcaV0daughters(), v0.pt(), v0.mK0Short());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeroanalysis>(cfgc, TaskName{"lf-lambdakzeroanalysis"}),
    adaptAnalysisTask<lambdakzeroQA>(cfgc, TaskName{"lf-lambdakzeroQA"})};
}
