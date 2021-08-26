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

// heavy flavours
DECLARE_SOA_COLUMN(HfHighPt, hasHfHighPt, bool); //! high-pT charm hadron
DECLARE_SOA_COLUMN(HfBeauty, hasHfBeauty, bool); //! beauty hadron
DECLARE_SOA_COLUMN(HfFemto, hasHfFemto, bool);   //! charm-hadron - N pair

// strangeness
DECLARE_SOA_COLUMN(Omega, hasOmega, bool); //! at leat 1 Omega
DECLARE_SOA_COLUMN(hadronXi, hashadronXi, bool); //! at least 1 Xi + high-pt hadron
DECLARE_SOA_COLUMN(DoubleXi, hasDoubleXi, bool); //! at least 2 Xi
DECLARE_SOA_COLUMN(TripleXi, hasTripleXi, bool); //! at least 3 Xi
DECLARE_SOA_COLUMN(QuadrupleXi, hasQuadrupleXi, bool); //! at least 4 Xi

} // namespace filtering

// nuclei
DECLARE_SOA_TABLE(NucleiFilters, "AOD", "NucleiFilters", //!
                  filtering::H2, filtering::H3, filtering::He3, filtering::He4);
using NucleiFilter = NucleiFilters::iterator;

// diffraction
DECLARE_SOA_TABLE(DiffractionFilters, "AOD", "DiffFilters", //! Diffraction filters
                  filtering::DG);
using DiffractionFilter = DiffractionFilters::iterator;

// heavy flavours
DECLARE_SOA_TABLE(HfFilters, "AOD", "HF Filters", //!
                  filtering::HfHighPt, filtering::HfBeauty, filtering::HfFemto);

using HfFilter = HfFilters::iterator;

// strangeness
DECLARE_SOA_TABLE(StrangenessFilters, "AOD", "StrgFilters", //!
                   filtering::Omega, filtering::hadronXi, filtering::DoubleXi, filtering::TripleXi, filtering::QuadrupleXi);

using StrangenessFilter = StrangenessFilters::iterator;

/// List of the available filters, the description of their tables and the name of the tasks
constexpr int NumberOfFilters{4};
constexpr std::array<char[32], NumberOfFilters> AvailableFilters{"NucleiFilters", "DiffractionFilters", "HeavyFlavourFilters", "StrangenessFilters"};
constexpr std::array<char[16], NumberOfFilters> FilterDescriptions{"NucleiFilters", "DiffFilters", "HFFilters", "StrgFilters"};
constexpr std::array<char[128], NumberOfFilters> FilteringTaskNames{"o2-analysis-nuclei-filter", "o2-analysis-diffraction-filter", "o2-analysis-hf-filter", "o2-analysis-strangeness-filter"};
constexpr o2::framework::pack<NucleiFilters, DiffractionFilters, HfFilters, StrangenessFilters> FiltersPack;
static_assert(o2::framework::pack_size(FiltersPack) == NumberOfFilters);

template <typename T, typename C>
void addColumnToMap(std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  map[MetadataTrait<T>::metadata::tableLabel()][C::columnLabel()] = 1.f;
}

template <typename T, typename... C>
void addColumnsToMap(o2::framework::pack<C...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  (addColumnToMap<T, C>(map), ...);
}

template <typename... T>
void FillFiltersMap(o2::framework::pack<T...>, std::unordered_map<std::string, std::unordered_map<std::string, float>>& map)
{
  (addColumnsToMap<T>(typename T::iterator::persistent_columns_t{}, map), ...);
}

template <typename... C>
static std::vector<std::string> ColumnsNames(o2::framework::pack<C...>)
{
  return {C::columnLabel()...};
}

template <typename T>
unsigned int NumberOfColumns()
{
  return o2::framework::pack_size(typename T::iterator::persistent_columns_t{});
}

} // namespace o2::aod

#endif // O2_ANALYSIS_TRIGGER_H_
