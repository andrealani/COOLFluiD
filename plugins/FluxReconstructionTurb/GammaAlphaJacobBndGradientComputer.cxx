// Copyright (C) 2019 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionTurb/GammaAlphaJacobBndGradientComputer.hh"
#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< GammaAlphaJacobBndGradientComputer,
		       FluxReconstructionSolverData,
		       FluxReconstructionNavierStokesModule >
GammaAlphaJacobBndGradientComputerProvider("ConvBndCorrectionsRHSJacobGammaAlpha");

//////////////////////////////////////////////////////////////////////////////

GammaAlphaJacobBndGradientComputer::GammaAlphaJacobBndGradientComputer(const std::string& name) :
  NSJacobBndGradientComputer(name),
  socket_volumes("volumes"),
  socket_solPntNormals("solPntNormals"),
  m_transitionCriterion()
{
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaJacobBndGradientComputer::setup()
{
  CFAUTOTRACE;

  // setup parent class
  NSJacobBndGradientComputer::setup();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);

  m_transitionCriterion.setup(*frLocalData[0],m_nbrEqs,m_dim,m_updateVarSet,m_diffusiveVarSet);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
GammaAlphaJacobBndGradientComputer::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = NSJacobBndGradientComputer::needsSockets();

  result.push_back(&socket_solPntNormals);
  result.push_back(&socket_volumes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaJacobBndGradientComputer::prepareGhostStates()
{
  m_transitionCriterion.setTransitionFlags(*m_bcStateComputer,*m_cellStates,m_cellStatesFlxPnt,m_flxPntGhostSol,
                                           m_unitNormalFlxPnts,m_flxPntCoords,m_faceJacobVecSizeFlxPnts,
                                           (*m_faceFlxPntConn)[m_orient],m_nbrFaceFlxPnts,m_corrFctDiv,
                                           socket_solPntNormals.getDataHandle(),socket_volumes.getDataHandle());
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
