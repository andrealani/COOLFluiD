#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "SpectralFD/BR2FaceTermComputer.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
  BR2FaceTermComputer,SpectralFDMethodData,BaseFaceTermComputer,SpectralFDModule >
BR2FaceTermComputerProvider("BR2FaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Alpha","Parameter of Bassi-Rebay 2 approach determining the amount of damping.");
}

//////////////////////////////////////////////////////////////////////////////

BR2FaceTermComputer::BR2FaceTermComputer(const std::string& name) :
  CompactFaceTermComputer(name),
  m_solPntsLocalCoords(CFNULL),
  m_solJacobDet(),
  m_backupPhysVarGrad(),
  m_alpha(),
  m_liftingOperators(),
  m_liftingOperatorsPtrs(),
  m_liftOperatorTerms(),
  m_liftingOperator(),
  m_backupLiftingOperatorTerms()
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_alpha = 0.0;
  setParameter("Alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

BR2FaceTermComputer::~BR2FaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::configure ( Config::ConfigArgs& args )
{
  CompactFaceTermComputer::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::setFaceTermData()
{
  // call setFaceTermData() of parent class
  CompactFaceTermComputer::setFaceTermData();

  // get the SpectralFDElementData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get solution point mapped coordinates for Jacobian determinant computation
  m_solPntsLocalCoords = sdLocalData[0]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::computeNeighbourCellData()
{
  // call computeNeighbourCellData() of parent class
  CompactFaceTermComputer::computeNeighbourCellData();

  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbouring cell
    GeometricEntity* cell = m_face->getNeighborGeo(iSide);

    // compute solution point Jacobian determinants
    m_solJacobDet[iSide] = cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::reconstructFluxPntsSolPolyGrads(const vector< vector< State* >* >& cellStates)
{
  // compute lifting operators
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx] = m_alpha*(*m_flxPntRVSol[LEFT][iFlx] - *m_flxPntRVSol[RIGHT][iFlx]);
  }
  computeGradFaceTermFromFlxPntSol(m_liftingOperators);

  cf_assert(cellStates.size() == 2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // compute polynomial gradients
    m_statesReconstr->computePolyGradients(*cellStates[iSide],m_flxPntGrads[iSide],
                                           *m_coefSolPolyDerivInFlxPnts,
                                           (*m_faceFlxPntConn)[m_orient][iSide],
                                           *m_solPntIdxsSolPolyDerivInFlxPnts,
                                           m_cellInvJacobMatr[iSide]);

    // subtract lifting operator terms
    const CFuint nbrSol = m_liftingOperators[iSide].size();
    if (iSide == LEFT)
    {
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_liftingOperators[iSide][iSol][iEq] /= m_solJacobDet[iSide][iSol];
        }
      }
    }
    else
    {
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          m_liftingOperators[iSide][iSol][iEq] /= -m_solJacobDet[iSide][iSol];
        }
      }
    }
    m_statesReconstr->reconstructGradients(m_liftingOperatorsPtrs[iSide],m_liftOperatorTerms[iSide],*m_flxPntsRecCoefs,
                                           (*m_faceFlxPntConn)[m_orient][iSide],
                                           *m_flxPntMatrixIdxForReconstruction,
                                           *m_solPntIdxsForReconstruction);
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        *m_flxPntGrads[iSide][iFlx][iEq] -= *m_liftOperatorTerms[iSide][iFlx][iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::backupAndReconstrFluxPntsSolPolyGrad(const CFuint side,
                                                               const CFuint iVar,
                                                               const vector< State* >& cellStates)
{
  // other side
  const CFuint otherSide = side == LEFT ? RIGHT : LEFT;

  // backup
  backupPhysVarGrad(side,iVar);

  // reconstruct polynomial gradient
  m_statesReconstr->computePolyGradient(iVar,cellStates,m_flxPntGrads[side],
                                        *m_coefSolPolyDerivInFlxPnts,
                                        (*m_faceFlxPntConn)[m_orient][side],
                                        *m_solPntIdxsSolPolyDerivInFlxPnts,
                                        m_cellInvJacobMatr[side]);

  // remove previous lifting operator term from other side gradient
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    *m_flxPntGrads[otherSide][iFlx][iVar] += *m_liftOperatorTerms[otherSide][iFlx][iVar];
  }

  // compute new lifting operator
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx][iVar] = m_alpha*((*m_flxPntRVSol[LEFT][iFlx])[iVar] - (*m_flxPntRVSol[RIGHT][iFlx])[iVar]);
  }
  computePhysVarGradFaceTermFromFlxPntSol(iVar,m_liftingOperator);

  // add new lifting operator terms
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    const CFuint nbrSol = m_liftingOperators[iSide].size();
    if (iSide == LEFT)
    {
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
      {
        m_liftingOperators[iSide][iSol][iVar] = m_liftingOperator[iSide][iSol]/m_solJacobDet[iSide][iSol];
      }
    }
    else
    {
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
      {
        m_liftingOperators[iSide][iSol][iVar] = m_liftingOperator[iSide][iSol]/(-m_solJacobDet[iSide][iSol]);
      }
    }
    m_statesReconstr->reconstructPhysVarGrad(iVar,m_liftingOperatorsPtrs[iSide],
                                             m_liftOperatorTerms[iSide],*m_flxPntsRecCoefs,
                                             (*m_faceFlxPntConn)[m_orient][iSide],
                                             *m_flxPntMatrixIdxForReconstruction,
                                             *m_solPntIdxsForReconstruction);
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      *m_flxPntGrads[iSide][iFlx][iVar] -= *m_liftOperatorTerms[iSide][iFlx][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::backupPhysVarGrad(const CFuint side, const CFuint iVar)
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad[iSide].size());
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      m_backupPhysVarGrad         [iSide][iFlx] = *m_flxPntGrads      [iSide][iFlx][iVar];
      m_backupLiftingOperatorTerms[iSide][iFlx] = *m_liftOperatorTerms[iSide][iFlx][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::restorePhysVarGrad(const CFuint side, const CFuint iVar)
{
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVarGrad[iSide].size());
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      *m_flxPntGrads      [iSide][iFlx][iVar] = m_backupPhysVarGrad         [iSide][iFlx];
      *m_liftOperatorTerms[iSide][iFlx][iVar] = m_backupLiftingOperatorTerms[iSide][iFlx];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  CompactFaceTermComputer::setup();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  // resize m_backupPhysVarGrad
  m_backupPhysVarGrad.resize(2);
  m_backupPhysVarGrad[LEFT ].resize(nbrFlxPnts,RealVector(m_dim));
  m_backupPhysVarGrad[RIGHT].resize(nbrFlxPnts,RealVector(m_dim));

  // check if damping value is sufficiently large
//   const CFuint nbrCellFaces = sdLocalData[0]->getNbrCellFaces();
//   m_alpha = nbrCellFaces > m_alpha ? nbrCellFaces : m_alpha;
//   if (sdLocalData[0]->getPolyOrder() == 0)
//   {
//     m_alpha = 1.0;
//   }
  m_alpha = 1.0;

  // resize m_liftOperatorTerms, m_liftOperatorTermPtrs and m_solJacobDet
  m_solJacobDet               .resize(2);
  m_liftingOperator           .resize(2);
  m_liftingOperators          .resize(2);
  m_liftingOperatorsPtrs      .resize(2);
  m_liftOperatorTerms         .resize(2);
  m_backupLiftingOperatorTerms.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_solJacobDet         [iSide].resize(nbrSolPnts);
    m_liftingOperator     [iSide].resize(nbrSolPnts,RealVector(m_dim));
    m_liftingOperators    [iSide].resize(nbrSolPnts);
    m_liftingOperatorsPtrs[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_liftingOperators    [iSide][iSol].resize(m_nbrEqs,RealVector(m_dim));
      m_liftingOperatorsPtrs[iSide][iSol] = &m_liftingOperators[iSide][iSol];
    }
    m_liftOperatorTerms         [iSide].resize(nbrFlxPnts);
    m_backupLiftingOperatorTerms[iSide].resize(nbrFlxPnts,RealVector(m_dim));
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_liftOperatorTerms[iSide][iFlx].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_liftOperatorTerms[iSide][iFlx][iEq] = new RealVector(m_dim);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BR2FaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_matrInverter);

  for (CFuint iSide = 0; iSide < m_liftOperatorTerms.size(); ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < m_liftOperatorTerms[iSide].size(); ++iFlx)
    {
      for (CFuint iEq = 0; iEq < m_liftOperatorTerms[iSide][iFlx].size(); ++iEq)
      {
        deletePtr(m_liftOperatorTerms[iSide][iFlx][iEq]);
      }
    }
  }

  // call unsetup of parent class
  CompactFaceTermComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
