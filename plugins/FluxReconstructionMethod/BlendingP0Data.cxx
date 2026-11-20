// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/BaseTerm.hh"
#include "Framework/PhysicalModel.hh"

#include "FluxReconstructionMethod/BlendingP0Data.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

BlendingP0Data::BlendingP0Data() :
  m_nbrSolPnts(0),
  m_nbrEqs(0),
  m_c(),
  m_w(),
  m_updateToSol(CFNULL),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_pData(),
  m_solAvg(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

BlendingP0Data::~BlendingP0Data()
{
  unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BlendingP0Data::setup(FluxReconstructionElementData* frLocalData,
                           const vector< vector< CFreal > >& corrFctDiv,
                           SafePtr< vector< vector< CFuint > > > flxSolDep,
                           const CFuint nbrEqs,
                           SafePtr< VarSetTransformer > updateToSol,
                           SafePtr< ConvectiveVarSet > solutionVar,
                           SafePtr< ConvectiveVarSet > updateVar)
{
  m_nbrSolPnts = frLocalData->getNbrOfSolPnts();
  m_nbrEqs = nbrEqs;

  // c_j from the element data
  SafePtr< RealVector > cTable = frLocalData->getCellAvgSolCoefs();
  cf_assert(cTable->size() == m_nbrSolPnts);
  m_c.resize(m_nbrSolPnts);
  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    m_c[j] = (*cTable)[j];
  }

  // W_fp = sum_j c_j divh_j,fp over the solution points fp depends on
  const CFuint nbrFlxPnts = flxSolDep->size();
  m_w.assign(nbrFlxPnts, 0.0);
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const vector< CFuint >& dep = (*flxSolDep)[iFlx];
    for (CFuint iSol = 0; iSol < dep.size(); ++iSol)
    {
      m_w[iFlx] += m_c[dep[iSol]]*corrFctDiv[dep[iSol]][iFlx];
    }
  }

  // update to solution variables by the transformer of the method data, back through the physical data
  m_updateToSol = updateToSol;
  m_solutionVar = solutionVar;
  m_updateVar = updateVar;
  SafePtr< PhysicalModel > physModel = PhysicalModelStack::getActive();
  physModel->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_pData);

  RealVector dummyCoord(physModel->getDim());
  dummyCoord = 0.0;
  m_solAvg = new State();
  m_solAvg->setSpaceCoordinates(new Node(dummyCoord, false));
}

//////////////////////////////////////////////////////////////////////////////

void BlendingP0Data::unsetup()
{
  // each state owns its coordinates (Node built with isOnMesh = false)
  deletePtr(m_solAvg);
}

//////////////////////////////////////////////////////////////////////////////

void BlendingP0Data::computeCellAverage(const vector< State* >& states,
                                        DataHandle< CFreal >& volumes,
                                        State& avg)
{
  cf_assert(states.size() == m_nbrSolPnts);

  *m_solAvg = 0.0;
  CFreal wSum = 0.0;
  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    // transform() returns a buffer overwritten by the next call, so use it straight away
    const State* const solState = m_updateToSol->transform(states[j]);
    const CFreal w = m_c[j]*volumes[states[j]->getLocalID()];
    wSum += w;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_solAvg)[iEq] += w*(*solState)[iEq];
    }
  }

  cf_assert(wSum > 0.0);
  *m_solAvg /= wSum;

  m_solutionVar->computePhysicalData(*m_solAvg, m_pData);
  m_updateVar->computeStateFromPhysicalData(m_pData, avg);
}

//////////////////////////////////////////////////////////////////////////////

CFreal BlendingP0Data::computeModeRelaxation(const vector< State* >& states,
                                             DataHandle< CFreal >& volumes,
                                             const CFreal alpha,
                                             const CFreal factor,
                                             const CFuint order,
                                             const CFreal h,
                                             const vector< RealVector >& normals,
                                             vector< RealVector >& res)
{
  cf_assert(states.size() == m_nbrSolPnts);
  cf_assert(res.size() == m_nbrSolPnts);

  // conserved average U_bar of the solution variables
  *m_solAvg = 0.0;
  CFreal wSum = 0.0;
  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    // transform() returns a buffer overwritten by the next call, so use it straight away
    const State* const solState = m_updateToSol->transform(states[j]);
    const CFreal w = m_c[j]*volumes[states[j]->getLocalID()];
    wSum += w;
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_solAvg)[iEq] += w*(*solState)[iEq];
    }
  }

  cf_assert(wSum > 0.0);
  *m_solAvg /= wSum;

  const CFreal sigma = computeModeRelaxationSigma(*m_solAvg, factor, order, h, normals);

  // -alpha sigma |J_j| (U_j - U_bar)
  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    const State* const solState = m_updateToSol->transform(states[j]);
    const CFreal f = alpha*sigma*volumes[states[j]->getLocalID()];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      res[j][iEq] = -f*((*solState)[iEq] - (*m_solAvg)[iEq]);
    }
  }

  return alpha*sigma;
}

//////////////////////////////////////////////////////////////////////////////

CFreal BlendingP0Data::computeModeRelaxationSigma(const RealVector& solAvg,
                                                  const CFreal factor,
                                                  const CFuint order,
                                                  const CFreal h,
                                                  const vector< RealVector >& normals)
{
  // largest |eigenvalue| lambda at U_bar over the normals; m_solAvg holds U_bar as a state for the physical data
  if (&solAvg != m_solAvg)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_solAvg)[iEq] = solAvg[iEq];
    }
  }
  m_solutionVar->computePhysicalData(*m_solAvg, m_pData);

  CFreal lambda = 0.0;
  for (CFuint i = 0; i < normals.size(); ++i)
  {
    lambda = std::max(lambda, m_updateVar->getMaxAbsEigenValue(m_pData, normals[i]));
  }

  return factor*(2.0*order + 1.0)*lambda/h;
}

//////////////////////////////////////////////////////////////////////////////

void BlendingP0Data::setModeRelaxationNormals(const CFuint dim, vector< RealVector >& normals)
{
  normals.clear();

  for (CFuint iDim = 0; iDim < dim; ++iDim)
  {
    RealVector n(dim);
    n = 0.0;
    n[iDim] = 1.0;
    normals.push_back(n);
  }

  if (dim == 2)
  {
    const CFreal s = 1.0/std::sqrt(2.0);
    RealVector n(2);
    n[0] = s; n[1] = s;
    normals.push_back(n);
    n[0] = s; n[1] = -s;
    normals.push_back(n);
  }
  else if (dim == 3)
  {
    const CFreal s = 1.0/std::sqrt(3.0);
    const CFreal signs[4][2] = {{1.0, 1.0}, {1.0, -1.0}, {-1.0, 1.0}, {-1.0, -1.0}};
    for (CFuint i = 0; i < 4; ++i)
    {
      RealVector n(3);
      n[0] = s; n[1] = signs[i][0]*s; n[2] = signs[i][1]*s;
      normals.push_back(n);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
