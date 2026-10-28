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

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
