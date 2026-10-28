// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/DiffRHSJacobFluxReconstructionNSSubcellBlending.hh"
#include "FluxReconstructionNavierStokes/DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending,
                       FluxReconstructionSolverData,
                       FluxReconstructionNavierStokesModule >
diffBndCorrectionsRHSJacobNSSubcellBlendingProvider("DiffBndCorrectionsRHSJacobNSSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending::DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending(const std::string& name) :
  DiffBndCorrectionsRHSJacobFluxReconstructionNS(name),
  socket_alpha("alpha"),
  m_closestSolToFlx(CFNULL),
  m_faceFactors()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = DiffBndCorrectionsRHSJacobFluxReconstructionNS::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending::computeCorrection(vector< RealVector >& corrections)
{
  DiffBndCorrectionsRHSJacobFluxReconstructionNS::computeCorrection(corrections);

  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  const CFreal a = alpha[(*m_cellStates)[0]->getLocalID()];
  if (a <= 0.0) return;

  // alpha F_f (c_f/w_f at the closest point - divh(., f)) for each face flux point f
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iFlx];
    const RealVector& faceFlux = m_cellFlx[iFlx];

    const CFuint closest = (*m_closestSolToFlx)[flxIdx];
    const CFreal factor = a*m_faceFactors[flxIdx];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[closest][iEq] += factor*faceFlux[iEq];
    }

    const std::vector< CFuint >& solDep = (*m_flxSolDep)[flxIdx];
    for (CFuint k = 0; k < solDep.size(); ++k)
    {
      const CFuint solIdx = solDep[k];
      const CFreal divh = a*m_corrFctDiv[solIdx][flxIdx];
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        corrections[solIdx][iEq] -= divh*faceFlux[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffBndCorrectionsRHSJacobFluxReconstructionNSSubcellBlending::setup()
{
  CFAUTOTRACE;

  DiffBndCorrectionsRHSJacobFluxReconstructionNS::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_closestSolToFlx = frLocalData[0]->getClosestSolToFlxIdx();
  DiffRHSJacobFluxReconstructionNSSubcellBlending::computeFaceSubcellFactors(frLocalData[0],m_corrFctDiv,m_faceFactors);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
