// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BaseCFMeshFileSource.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

BaseCFMeshFileSource::BaseCFMeshFileSource() :
  _isWithSolution(false),
  _dimension(0),
  _nbEquations(0),
  _nbUpdatableNodes(0),
  _nbNonUpdatableNodes(0),
  _nbUpdatableStates(0),
  _nbNonUpdatableStates(0),
  _nbElements(0),
  _nbElementTypes(0),
  _nbTRSs(0),
  _geometricPolyOrder(CFPolyOrder::ORDER0),
  _solutionPolyOrder(CFPolyOrder::ORDER0),
  _geometricPolyType(CFPolyForm::INVALID),
  _solutionPolyType(CFPolyForm::INVALID),
  _elementTypeData(),
  _nameTRS(),
  _nbTRs(),
  _nbGeomEntsPerTR(),
  _geomType(),
  _geoConn(),
  _nbExtraVars(0),
  _extraVarNames(),
  _nbExtraNodalVars(0),
  _nbExtraStateVars(0),
  _extraNodalVarNames(),
  _extraStateVarNames(),
  _nbGroups(0),
  _groupNames(),
  _groupsNbElem(),
  _groupElemList(),
  _storePastStates(false),
  _storePastNodes(false),
  _storeInterStates(false),
  _storeInterNodes(false),
  dynamicSockets(CFNULL),
  _extraNVarMap(),
  _extraSVarMap(),
  _extraVarMap(),
  _extraStateVector(),
  _extraNodeVector(),
  _extraVector()

{
  // by default we will put the polynomial type to be CFPolyForm::LAGRANGE
  // not to have to change all the CFmesh files created

  /// @todo eventually CFmesh should evolve to a CFmesh2 file format
  /// where this information about elements is mandatory

  _geometricPolyType = CFPolyForm::LAGRANGE;
  _solutionPolyType  = CFPolyForm::LAGRANGE;

}

//////////////////////////////////////////////////////////////////////////////

void BaseCFMeshFileSource::releaseMemory()
{
  vector<TRGeoConn>().swap(_geoConn);
  vector<CFGeoEnt::Type>().swap(_geomType);
  vector<vector<CFuint> >().swap(_nbGeomEntsPerTR);
  vector<CFuint>().swap(_nbTRs);
  vector<std::string>().swap(_nameTRS);
  vector<ElementTypeData>().swap(_elementTypeData);

  _isWithSolution = 0;
  _dimension = 0;
  _nbEquations = 0;
  _nbUpdatableNodes = 0;
  _nbNonUpdatableNodes = 0;
  _nbUpdatableStates = 0;
  _nbNonUpdatableStates = 0;
  _nbElements = 0;
  _nbElementTypes = 0;
  _nbTRSs = 0;
  _nbExtraNodalVars = 0;
  _nbExtraStateVars = 0;
  _nbExtraVars = 0;
  _nbGroups = 0;
}

//////////////////////////////////////////////////////////////////////////////

bool BaseCFMeshFileSource::consistencyCheck() const
{
  cf_assert(_nbElementTypes == _elementTypeData.size());

  CFuint sum = 0;
  for(CFuint i = 0; i< _nbElementTypes; ++i) {
    sum += _elementTypeData[i].getNbElems();
  }
  cf_assert(sum == _nbElements);
//CFout << _nameTRS.size() << " " << _nbTRSs << "\n";
  cf_assert(_nameTRS.size() == _nbTRSs);
  cf_assert(_geoConn.size() == _nbTRSs);
  cf_assert(_nbTRs.size()   == _nbTRSs);
  cf_assert(_nbGeomEntsPerTR.size() == _nbTRSs);

  for(CFuint iTRS = 0; iTRS < _nbTRSs; ++iTRS ) {
    cf_assert(_nbTRs[iTRS] == _nbGeomEntsPerTR[iTRS].size());
    cf_assert(_nbTRs[iTRS] == _geoConn[iTRS].size());
    for(CFuint iTR = 0; iTR < _nbTRs[iTRS]; ++iTR ) {
      cf_assert(_geoConn[iTRS][iTR].size() == _nbGeomEntsPerTR[iTRS][iTR]);
    }
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
