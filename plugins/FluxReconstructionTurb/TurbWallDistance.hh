// Copyright (C) 2026 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_TurbWallDistance_hh
#define COOLFluiD_FluxReconstructionMethod_TurbWallDistance_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataHandle.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Wall distance rules of the turbulence commands. The turbulent diffusive
 * variable set reads the wall distance (NavierStokesVarSet::setWallDistance)
 * when it evaluates a diffusive flux, so every command sets it for the point
 * the flux is evaluated at: at a solution point the wall distance of its state;
 * at a flux point of a face the wall distance of the solution point of the cell
 * closest to that flux point (FluxReconstructionElementData::getClosestSolToFlxIdx),
 * averaged over the two cells at an interior face, the interior cell's alone at
 * a boundary face.
 *
 * @author Rayan Dhib
 */

/**
 * Wall distance of one cell at a flux point: the wall distance of the solution
 * point of the cell closest to the flux point.
 * @param wallDist           wall distance of every state
 * @param states             states of the cell
 * @param closestSolToFlxIdx index of the closest solution point of every flux point of the cell
 * @param flxPntIdx          index of the flux point among the flux points of the cell
 */
inline CFreal wallDistanceAtFlxPnt(const Framework::DataHandle< CFreal >& wallDist,
                                   const std::vector< Framework::State* >& states,
                                   const std::vector< CFuint >& closestSolToFlxIdx,
                                   const CFuint flxPntIdx)
{
  return wallDist[states[closestSolToFlxIdx[flxPntIdx]]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_TurbWallDistance_hh
