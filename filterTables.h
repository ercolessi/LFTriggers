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
#ifndef O2_ANALYSIS_TRIGGER_H_
#define O2_ANALYSIS_TRIGGER_H_

#include <array>
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace filtering
{
DECLARE_SOA_COLUMN(H2, hasH2, bool);   //!
DECLARE_SOA_COLUMN(H3, hasH3, bool);   //!
DECLARE_SOA_COLUMN(He3, hasHe3, bool); //!
DECLARE_SOA_COLUMN(He4, hasHe4, bool); //!

// diffraction
DECLARE_SOA_COLUMN(DG, hasDG, bool); //! Double Gap events, DG

//Strangeness
 DECLARE_SOA_COLUMN(DXi, hasDXi, bool);
 DECLARE_SOA_COLUMN(hadronXi, hashadronXi, bool);

} // namespace filtering

DECLARE_SOA_TABLE(NucleiFilters, "AOD", "Nuclei Filters", //!
                  filtering::H2, filtering::H3, filtering::He3, filtering::He4);

 constexpr std::array<char[32], 3> AvailableFilters{"NucleiFilters", "DiffractionFilters", "StrangenessFilters"};
 constexpr std::array<char[16], 3> FilterDescriptions{"Nuclei Filters", "DiffFilters", "Strg Filters"};

using NucleiFilter = NucleiFilters::iterator;

// diffraction
DECLARE_SOA_TABLE(DiffractionFilters, "AOD", "DiffFilters", //! Diffraction filters
                  filtering::DG);
using DiffractionFilter = DiffractionFilters::iterator;

//Strangeness
 DECLARE_SOA_TABLE(StrangenessFilters, "AOD", "Strg Filters", filtering::DXi, filtering::hadronXi);
 using StrangenessFilter = StrangenessFilters::iterator;

} // namespace o2::aod

#endif // O2_ANALYSIS_TRIGGER_H_
