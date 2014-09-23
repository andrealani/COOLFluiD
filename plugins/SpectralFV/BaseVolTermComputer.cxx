#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/BaseVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

BaseVolTermComputer::BaseVolTermComputer(const std::string& name) :
  SpectralFVMethodStrategy(name),
  socket_faceNormTransfMatrices("faceNormTransfMatrices"),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_flxPntsRecCoefs(CFNULL),
  m_cellFaceNormTransfM(),
  m_solInFlxPnts(),
  m_extraVarsInFlxPnts(),
  m_solRVInFlxPnts(),
  m_gradVarsInFlxPnts(),
  m_gradInFlxPnts(),
  m_backupPhysVar(),
  m_nbrFlxPnts(),
  m_nbrEqs(),
  m_nbrExtraVars()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseVolTermComputer::~BaseVolTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellData()
{
  // get faceNormTransfMatrices
  DataHandle< RealMatrix > faceNormTransfMatrices = socket_faceNormTransfMatrices.getDataHandle();

  // set the cell face normal transformation matrix
  m_cellFaceNormTransfM = faceNormTransfMatrices[m_cell->getID()];
}


//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::reconstructStates(const vector< State* >& cellStates)
{
  m_statesReconstr->reconstructStates(cellStates,m_solInFlxPnts,
                                      *m_flxPntsRecCoefs,
                                      cellStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::reconstructGradients(const vector< vector< RealVector >* >& cellGradients,
                                               const CFuint nbrCellGrads)
{
  m_statesReconstr->reconstructGradients(cellGradients,
                                         m_gradInFlxPnts,
                                         *m_flxPntsRecCoefs,
                                         nbrCellGrads);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::backupAndReconstructPhysVar(const CFuint iVar, const vector< State* >& cellStates)
{
  // backup
  backupPhysVar(iVar);

  // reconstruct
  m_statesReconstr->reconstructPhysVar(iVar,cellStates,m_solInFlxPnts,
                                       *m_flxPntsRecCoefs,
                                       cellStates.size());
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::backupPhysVar(const CFuint iVar)
{
  const CFuint nbrFlxPnts = m_solInFlxPnts.size();
  cf_assert(nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_solInFlxPnts[iFlx])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::restorePhysVar(const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    (*m_solInFlxPnts[iFlx])[iVar] = m_backupPhysVar[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellConvVolumeTerm(RealVector& resUpdates)
{
  // compute the states data in the flux points
 /// @todo broken after release 2009.3
//   m_updateVarSet->computeStatesData(m_solInFlxPnts,m_extraVarsInFlxPnts,m_nbrFlxPnts);

  // compute the actual volume term
  computeConvVolTermFromFlxPntSol(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeGradientVolumeTerm(vector< vector< RealVector > >& gradUpdates)
{
  // set the gradient variables in the flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
  {
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      m_gradVarsInFlxPnts(iGrad,iFlx) = (*m_solRVInFlxPnts[iFlx])[iGrad];
    }
  }

  // compute the actual volume term contribution to the gradient
  computeGradVolTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::computeCellDiffVolumeTerm(RealVector& resUpdates)
{
  // compute the actual volume term
  computeDiffVolTermFromFlxPntSolAndGrad(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseVolTermComputer::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of variables
//   const CFuint dim     = PhysicalModelStack::getActive()->getDim ();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  // get number of extra variables from the update variable set
  m_nbrExtraVars = m_updateVarSet->getExtraPhysicalVarsSize();

  // get the diffusive varset
  if (getMethodData().hasDiffTerm())
  {
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
    BaseVolTermComputer::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_faceNormTransfMatrices);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

