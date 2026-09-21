// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionTurb/DiffBndCorrectionsRHSJacobFluxReconstructionTurb.hh"
#include "FluxReconstructionTurb/TurbWallDistance.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
diffBndCorrectionsRHSJacobTurbFluxReconstructionProvider("DiffBndCorrectionsRHSJacobTurb");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionTurb::DiffBndCorrectionsRHSJacobFluxReconstructionTurb(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL),
  m_navierStokesVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionTurb::~DiffBndCorrectionsRHSJacobFluxReconstructionTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;

  // setup parent class
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);

  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();

  m_navierStokesVarSet = m_diffusiveVarSet.d_castTo< NavierStokesVarSet >();
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::unsetup()
{
  CFAUTOTRACE;

  // unsetup parent class
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    DiffBndCorrectionsRHSJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionTurb::prepareFlxPntFluxComputation(const CFuint iFlx)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  // the interior cell, at its solution point closest to the flux point
  m_navierStokesVarSet->setWallDistance(wallDistanceAtFlxPnt(wallDist,*m_cellStates,*m_closestSolToFlxIdx,
                                                             (*m_faceFlxPntConn)[m_orient][iFlx]));

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
