// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

MapGeoToTrsAndIdx::MapGeoToTrsAndIdx() :
  m_trs(),
  m_idxInTrs(),
  m_isBGeo()
{
}

//////////////////////////////////////////////////////////////////////////////

MapGeoToTrsAndIdx::~MapGeoToTrsAndIdx()
{
}

//////////////////////////////////////////////////////////////////////////////

void MapGeoToTrsAndIdx::resize(const CFuint totalNbGeos)
{
  deallocate();
  m_trs.resize(totalNbGeos);
  m_idxInTrs.resize(totalNbGeos);
  m_isBGeo.resize(totalNbGeos);
}

//////////////////////////////////////////////////////////////////////////////

void MapGeoToTrsAndIdx::setMappingData(CFuint geoID,
		      Common::SafePtr<TopologicalRegionSet> trs,
		      CFuint idxInTrs,
		      bool isBGeo)
{
  cf_assert(geoID < m_trs.size());
  cf_assert(geoID < m_idxInTrs.size());
  cf_assert(geoID < m_isBGeo.size());
  m_trs[geoID] = trs;
  m_idxInTrs[geoID] = idxInTrs;
  m_isBGeo[geoID] = isBGeo;
}

//////////////////////////////////////////////////////////////////////////////

void MapGeoToTrsAndIdx::deallocate()
{
  std::vector< Common::SafePtr<TopologicalRegionSet> >().swap(m_trs);
  std::vector<CFuint>().swap(m_idxInTrs);
  std::vector<bool>().swap(m_isBGeo);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

