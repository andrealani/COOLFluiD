// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_GradientVariables_hh
#define COOLFluiD_FluxReconstructionMethod_GradientVariables_hh

//////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

#include "Common/COOLFluiD.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/State.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Kernels shared by the gradient computations of the FR diffusive terms:
 * the gradient variables of a cell extrapolated to a flux point, the volume
 * term of their gradient, the lifting of a face jump, the jumps at interior
 * and boundary flux points, the extrapolation of a gradient to the flux
 * points of a face, and the variables the artificial viscosity takes the
 * gradient of. The residual commands and the Jacobian commands call the same
 * functions, so both compute the same gradient.
 *
 * @author Rayan Dhib
 */

//////////////////////////////////////////////////////////////////////////////

/**
 * Extrapolates the gradient variables of a cell from its solution points to
 * one flux point with the solution polynomial:
 *
 *   g^D_f = sum_{j in solDep} phi_j(f) * g(:,j)
 *
 * The volume term of the gradient differentiates the same polynomial, so this
 * is also the value g^D_f the face correction starts from.
 *
 * @param gradVarsSolPnts   gradient variables at the solution points, nbrEqs x nbrSolPnts
 * @param solDep            solution points this flux point depends on
 * @param solPolyValsAtFlx  basis values at this flux point, indexed by solution point
 * @param nbrEqs            number of equations
 * @param gradVarsFlxPnt    the gradient variables at the flux point, size nbrEqs
 */
inline void extrapolateGradVarsToFlxPnt(const RealMatrix& gradVarsSolPnts,
                                        const std::vector< CFuint >& solDep,
                                        const std::vector< CFreal >& solPolyValsAtFlx,
                                        const CFuint nbrEqs,
                                        RealVector& gradVarsFlxPnt)
{
  gradVarsFlxPnt = 0.0;

  const CFuint nbrSolDep = solDep.size();

  for (CFuint iSol = 0; iSol < nbrSolDep; ++iSol)
  {
    const CFuint solIdx = solDep[iSol];
    const CFreal solPolyVal = solPolyValsAtFlx[solIdx];

    for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
    {
      gradVarsFlxPnt[iEq] += solPolyVal*gradVarsSolPnts(iEq,solIdx);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Adds the volume term of the gradient of the gradient variables of one cell
 * to gradUpdates, before the division by the Jacobian determinant:
 *
 *   gradUpdates[j][iEq] += sum_i sum_dir dphi_j/dxi_dir(x_i) * g(iEq,i) * S_dir(x_i)
 *
 * with S_dir the mapped coordinate plane normals at the solution points.
 *
 * @param gradVarsSolPnts        gradient variables at the solution points, nbrEqs x nbrSolPnts
 * @param fluxProjVects          mapped coordinate plane normals [dir][solPnt]
 * @param solSolDep              solution point dependencies of each solution point
 * @param solPolyDerivAtSolPnts  basis derivatives [solPnt][dir][solPnt]
 * @param nbrSolSolDep           number of dependencies per solution point
 * @param nbrSolPnts             number of solution points
 * @param nbrEqs                 number of equations
 * @param dim                    number of dimensions
 * @param projectedCorr          scratch vector of size dim
 * @param gradUpdates            updates [solPnt][eq], added to
 */
inline void addGradVarsVolumeTerm(const RealMatrix& gradVarsSolPnts,
                                  const std::vector< std::vector< RealVector > >& fluxProjVects,
                                  const std::vector< std::vector< CFuint > >& solSolDep,
                                  const std::vector< std::vector< std::vector< CFreal > > >& solPolyDerivAtSolPnts,
                                  const CFuint nbrSolSolDep,
                                  const CFuint nbrSolPnts,
                                  const CFuint nbrEqs,
                                  const CFuint dim,
                                  RealVector& projectedCorr,
                                  std::vector< std::vector< RealVector > >& gradUpdates)
{
  for (CFuint iSolPnt = 0; iSolPnt < nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
    {
      for (CFuint iDir = 0; iDir < dim; ++iDir)
      {
        projectedCorr = gradVarsSolPnts(iEq,iSolPnt) * fluxProjVects[iDir][iSolPnt];

        for (CFuint jSolPnt = 0; jSolPnt < nbrSolSolDep; ++jSolPnt)
        {
          const CFuint jSolIdx = solSolDep[iSolPnt][jSolPnt];

          gradUpdates[jSolIdx][iEq] += solPolyDerivAtSolPnts[jSolIdx][iDir][iSolPnt]*projectedCorr;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Adds the correction of one equation at one flux point to gradUpdates:
 *
 *   gradUpdates[s][iEq] += gradVarsJump * faceJacobVecSize * n * divh_s(f)
 *
 * for the solution points s the flux point depends on: the term
 * (g^I_f - g^D_f) grad h_f of the reconstructed gradient. gradVarsJump is the
 * jump g^I_f - g^D_f at this flux point (interface value minus extrapolated
 * value, times the BR2 multiplier eta for a compact face gradient),
 * faceJacobVecSize the signed face Jacobian vector size of this side, n the
 * unit normal and divh_s(f) the divergence of the correction function of the
 * flux point at solution point s.
 *
 * @param nbrSolDep      number of solution points to loop over in solDep
 * @param projectedCorr  scratch vector of size dim
 */
inline void addGradVarsLifting(const CFreal gradVarsJump,
                               const CFreal faceJacobVecSize,
                               const RealVector& unitNormal,
                               const std::vector< CFuint >& solDep,
                               const CFuint nbrSolDep,
                               const std::vector< std::vector< CFreal > >& corrFctDiv,
                               const CFuint flxIdx,
                               const CFuint iEq,
                               RealVector& projectedCorr,
                               std::vector< std::vector< RealVector > >& gradUpdates)
{
  projectedCorr = gradVarsJump*faceJacobVecSize*unitNormal;

  for (CFuint iSolPnt = 0; iSolPnt < nbrSolDep; ++iSolPnt)
  {
    const CFuint iSolIdx = solDep[iSolPnt];

    gradUpdates[iSolIdx][iEq] += projectedCorr*corrFctDiv[iSolIdx][flxIdx];
  }
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Jump g^I_f - g^D_f of the gradient variables at an interior flux point: the
 * interface value is the average of the values extrapolated from both sides,
 * and this side jumps from its own value to it.
 *
 * @param ownGradVars    value extrapolated from this side
 * @param otherGradVars  value extrapolated from the other side
 */
inline CFreal interiorGradVarsJump(const CFreal ownGradVars, const CFreal otherGradVars)
{
  const CFreal avgGradVars = (ownGradVars+otherGradVars)/2.0;
  return avgGradVars-ownGradVars;
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Jump g_b - g^D_f of the gradient variables at a boundary flux point: from the
 * extrapolated value to the boundary value of the boundary condition.
 *
 * @param gradVars     value extrapolated from the cell
 * @param bndGradVars  boundary value
 */
inline CFreal bndGradVarsJump(const CFreal gradVars, const CFreal bndGradVars)
{
  const CFreal ghostGradVars = 2.0*bndGradVars - gradVars;
  const CFreal avgGradVars = (gradVars+ghostGradVars)/2.0;
  return avgGradVars-gradVars;
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Extrapolates a gradient stored at the solution points of a cell to the flux
 * points of one face:
 *
 *   gradsFlxPnt[iFlx][iEq] = sum_s phi_s(f) * gradsSolPnts[s][iEq]
 *
 * @param flxPntConn      cell flux point index of each face flux point
 * @param nbrFaceFlxPnts  number of flux points of the face
 */
inline void extrapolateGradsToFlxPnts(const std::vector< std::vector< RealVector > >& gradsSolPnts,
                                      const std::vector< CFuint >& flxPntConn,
                                      const CFuint nbrFaceFlxPnts,
                                      const std::vector< std::vector< CFuint > >& flxSolDep,
                                      const std::vector< std::vector< CFreal > >& solPolyValsAtFlxPnts,
                                      const CFuint nbrEqs,
                                      std::vector< std::vector< RealVector* > >& gradsFlxPnt)
{
  for (CFuint iFlxPnt = 0; iFlxPnt < nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx = flxPntConn[iFlxPnt];

    for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
    {
      *(gradsFlxPnt[iFlxPnt][iEq]) = 0.0;
    }

    const CFuint nbrSolDep = flxSolDep[flxIdx].size();

    for (CFuint iSolPnt = 0; iSolPnt < nbrSolDep; ++iSolPnt)
    {
      const CFuint iSolIdx = flxSolDep[flxIdx][iSolPnt];

      for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
      {
        *(gradsFlxPnt[iFlxPnt][iEq]) += solPolyValsAtFlxPnts[flxIdx][iSolIdx]*gradsSolPnts[iSolIdx][iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Sets the variables whose gradients the artificial viscosity uses, one column
 * per state: the states themselves for Cons update variables, the solution
 * (conservative) variables for Puvt and RhoivtTv. Returns false, and leaves
 * values untouched, when no artificial viscosity is added or the update
 * variables are none of these.
 *
 * @param nbrStates                 number of states to set
 * @param hasArtificialViscosity    true when an artificial viscosity is added
 * @param updateVarStr              name of the update variables
 * @param updateToSolutionVecTrans  transformer from update to solution variables
 * @param nbrEqs                    number of equations
 * @param values                    output, nbrEqs x nbrStates
 */
inline bool setAVGradientVars(const std::vector< Framework::State* >& states,
                              const CFuint nbrStates,
                              const bool hasArtificialViscosity,
                              const std::string& updateVarStr,
                              Common::SafePtr< Framework::VarSetTransformer > updateToSolutionVecTrans,
                              const CFuint nbrEqs,
                              RealMatrix& values)
{
  if (!hasArtificialViscosity)
  {
    return false;
  }

  if (updateVarStr == "Cons")
  {
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
      {
        values(iEq,iState) = (*(states[iState]))[iEq];
      }
    }
    return true;
  }

  if (updateVarStr == "Puvt" || updateVarStr == "RhoivtTv")
  {
    for (CFuint iState = 0; iState < nbrStates; ++iState)
    {
      const RealVector transformedState = static_cast<RealVector&>(*updateToSolutionVecTrans->transform(states[iState]));

      for (CFuint iEq = 0; iEq < nbrEqs; ++iEq)
      {
        values(iEq,iState) = transformedState[iEq];
      }
    }
    return true;
  }

  return false;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GradientVariables_hh
