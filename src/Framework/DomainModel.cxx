// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BadFormatException.hh"
#include "Framework/DomainModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void DomainModel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("MapTopologicalRegions","The names of the TRSs and the local index of the TRs, for each TR.");
  options.addConfigOption< std::vector<CFuint> >("MapTopologicalRegionsCADids","Coresponding CAD ids for TRs from MapTopologicalRegions");
}

//////////////////////////////////////////////////////////////////////////////

DomainModel::DomainModel(const std::string& name) :
  OwnedObject(),
  ConfigObject(name),
  NonCopyable<DomainModel>(),
  m_trsNamesAndTRIdxs(),
  m_mapTRSName2TRIdx()
{
  addConfigOptionsTo(this);

  m_trsNamesAndTRIdxs = vector<std::string>();
  setParameter("MapTopologicalRegions",&m_trsNamesAndTRIdxs);
  setParameter("MapTopologicalRegionsCADids",&m_tabCADid);

}

//////////////////////////////////////////////////////////////////////////////

DomainModel::~DomainModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void DomainModel::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);

  // check if the size of m_trNames
  if ( m_trsNamesAndTRIdxs.size() % 2 != 0)
  {
    throw BadFormatException (FromHere(),"DomainModel::configure(): The size of MapTopologicalRegions is odd!!!");
  }

  // number of TRs
  const CFuint nbrTRs = m_trsNamesAndTRIdxs.size()/2;


  // if there are no CAD ids in CFcase then use consecutive numbers
  if ( m_tabCADid.size() != nbrTRs)
  {
	m_tabCADid.resize( nbrTRs);
	for ( CFuint i=0; i<m_tabCADid.size(); ++i)
		m_tabCADid[i] = i;
  }

  // create map from TRS name to TR ID
  for (CFuint iTR = 0; iTR < nbrTRs; ++iTR)
  {
    // create a pair < TRSname, TR idx >
    const std::string trKey = m_trsNamesAndTRIdxs[2*iTR] + m_trsNamesAndTRIdxs[2*iTR+1];

    // insert in the map
    m_mapTRSName2TRIdx.insert(trKey,iTR);
  }

  // sort the map
  m_mapTRSName2TRIdx.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
