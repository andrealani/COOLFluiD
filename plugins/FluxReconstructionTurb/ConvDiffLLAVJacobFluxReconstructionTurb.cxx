// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"

#include "FluxReconstructionTurb/ConvDiffLLAVJacobFluxReconstructionTurb.hh"
#include "FluxReconstructionTurb/TurbWallDistance.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvDiffLLAVJacobFluxReconstructionTurb,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
convDiffLLAVRHSJacobTurbFluxReconstructionProvider("ConvDiffLLAVRHSJacobTurb");

//////////////////////////////////////////////////////////////////////////////

ConvDiffLLAVJacobFluxReconstructionTurb::ConvDiffLLAVJacobFluxReconstructionTurb(const std::string& name) :
  ConvDiffLLAVJacobFluxReconstructionNS(name),
  socket_wallDistance("wallDistance"),
  m_closestSolToFlxIdx(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::configure ( Config::ConfigArgs& args )
{
  ConvDiffLLAVJacobFluxReconstructionNS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    ConvDiffLLAVJacobFluxReconstructionTurb::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = ConvDiffLLAVJacobFluxReconstructionNS::needsSockets();

  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::prepareSolPntFluxComputation(const CFuint stateID)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  m_diffusiveVarSetNS->setWallDistance(wallDist[stateID]);

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::prepareFlxPntFluxComputation(const CFuint iFlx)
{
  DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

  // average of the two cells, each taken at its solution point closest to the flux point
  const CFreal wallDistL = wallDistanceAtFlxPnt(wallDist,*m_states[LEFT],*m_closestSolToFlxIdx,
                                                (*m_faceFlxPntConnPerOrient)[m_orient][LEFT][iFlx]);
  const CFreal wallDistR = wallDistanceAtFlxPnt(wallDist,*m_states[RIGHT],*m_closestSolToFlxIdx,
                                                (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlx]);

  m_diffusiveVarSetNS->setWallDistance(0.5*(wallDistL + wallDistR));

  prepareFluxComputation();
}

//////////////////////////////////////////////////////////////////////////////

void ConvDiffLLAVJacobFluxReconstructionTurb::setup()
{
  CFAUTOTRACE;

  // setup parent class
  ConvDiffLLAVJacobFluxReconstructionNS::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);

  m_closestSolToFlxIdx = frLocalData[0]->getClosestSolToFlxIdx();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
