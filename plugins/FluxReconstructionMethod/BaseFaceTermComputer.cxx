#include "Common/BadValueException.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/ConsistencyException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/BaseFaceTermComputer.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BaseFaceTermComputer,FluxReconstructionSolverData,BaseFaceTermComputer,FluxReconstructionModule >
    BaseFaceTermComputerProvider("BaseFaceTermComputer");

//////////////////////////////////////////////////////////////////////////////

BaseFaceTermComputer::BaseFaceTermComputer(const std::string& name) :
  FluxReconstructionSolverStrategy(name),
  socket_extraVars("meanflow",false),/// @todo get the name of this socket in another way, not hardcoded, to make it more general
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_statesReconstr(CFNULL),
  m_riemannFluxComputer(CFNULL),
  //m_faceDiffFluxComputer(CFNULL),
  m_cflConvDiffRatio(),
  m_flxPntsRecCoefs(CFNULL),
  m_solPntsDerivCoefs(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_faceMappedCoordDir(CFNULL),
  m_flxPntMatrixIdxForReconstruction(CFNULL),
  m_solPntIdxsForReconstruction(CFNULL),
  m_flxPntMatrixIdxForDerivation(CFNULL),
  m_solPntIdxsForDerivation(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_faceFlxPntsFaceLocalCoords(CFNULL),
  m_faceFlxPntCellMappedCoords(CFNULL),
  m_face(CFNULL),
  m_orient(),
  m_cellVolumes(),
  m_unitNormalFlxPnts(),
  m_faceJacobVecSizeFlxPnts(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceInvCharLengths(),
  m_cellExtraVars(),
  m_flxPntSol(),
  m_allSol(),
  m_flxPntExtraVars(),
  m_allExtraVars(),
  m_flxPntRVSol(),
  m_flxPntRiemannFlux(),
  m_flxPntGrads(),
  m_flxPntGradPtrs(),
  m_backupPhysVar(),
  m_nbrFlxPnts(),
  m_gradTerm(),
  m_nbrEqs(),
  m_dim(),
  m_nbrExtraVars()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BaseFaceTermComputer::~BaseFaceTermComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverStrategy::configure(args);

  // initialize socket_extraVars
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::setFaceTermData()
{
  // get the FluxReconstructionElementData
  vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();

  // get derivation coefficients for the solution points
  m_solPntsDerivCoefs = sdLocalData[0]->getDerivCoefsSolPnts1D();

  // get indexes of internal flux points
  m_faceFlxPntConn = sdLocalData[0]->getFaceFlxPntConnPerOrient();

  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = sdLocalData[0]->getFaceMappedCoordDirPerOrient();

  // get flux point index (in the matrix flxPntRecCoefs) for reconstruction
  m_flxPntMatrixIdxForReconstruction = sdLocalData[0]->getFlxPntMatrixIdxForReconstruction();

  // get flux point index (in the matrix m_solPntsDerivCoefs) for derivation
  m_flxPntMatrixIdxForDerivation = sdLocalData[0]->getFlxPntMatrixIdxForDerivation();

  // get solution point index (in the cell) for derivation
  m_solPntIdxsForDerivation = sdLocalData[0]->getSolPntIdxsForDerivation();

  // get coefficients for integration over a face
  m_faceIntegrationCoefs = sdLocalData[0]->getFaceIntegrationCoefs();

  // get face flux points face local coordinates
  m_faceFlxPntsFaceLocalCoords = sdLocalData[0]->getFaceFlxPntsFaceLocalCoords();

  // get convective/diffusive CFL ratio
  m_cflConvDiffRatio = sdLocalData[0]->getCFLConvDiffRatio();

  // set number of flux points
  const CFuint nbrOrients = m_faceFlxPntConn->size();
  m_nbrFlxPnts.resize(nbrOrients);
  for (CFuint iOrient = 0; iOrient < nbrOrients; ++iOrient)
  {
    m_nbrFlxPnts[iOrient] = (*m_faceFlxPntConn)[iOrient][LEFT].size();
  }

//   // set interpolation type dependent data
//   const std::string interpolationType = getMethodData().getInterpolationType();
//   if (interpolationType == "standard")
//   {
//     // get reconstruction coefficients for the flux points
//     m_flxPntsRecCoefs = sdLocalData[0]->getRecCoefsFlxPnts1D();
// 
//     // get solution point index (in the cell) for reconstruction
//     m_solPntIdxsForReconstruction = sdLocalData[0]->getSolPntIdxsForReconstruction();
//   }
//   else if (interpolationType == "optimized")
//   {
//     // get reconstruction coefficients for the flux points
//     m_flxPntsRecCoefs = sdLocalData[0]->getRecCoefsFlxPnts1DOptim();
// 
//     // get solution point index (in the cell) for reconstruction
//     m_solPntIdxsForReconstruction = sdLocalData[0]->getSolPntIdxsForRecOptim();
//   }
//   else
//   {
//     throw BadValueException (FromHere(),"BaseVolTermComputer::setFaceTermData --> Interpolation type should be standard or optimized");
//   }

  // get face flux point cell mapped coordinates
  m_faceFlxPntCellMappedCoords = sdLocalData[0]->getFaceFlxPntCellMappedCoordsPerOrient();
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeFaceData()
{
  // face ID
  const CFuint faceID = m_face->getID();

  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_face->computeFaceJacobDetVectorAtMappedCoords(*m_faceFlxPntsFaceLocalCoords);

  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > >
      faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();

  // set face Jacobian vector sizes and unit normals in the flux points
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlx] = faceJacobVecSizeFaceFlxPnts[faceID][iFlx];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlx][LEFT ] =
        m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][LEFT ];
    m_faceJacobVecSizeFlxPnts[iFlx][RIGHT] =
        m_faceJacobVecAbsSizeFlxPnts[iFlx]*(*m_faceMappedCoordDir)[m_orient][RIGHT];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlx] = faceJacobVecs[iFlx]/m_faceJacobVecAbsSizeFlxPnts[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeNeighbourCellData()
{
  vector< std::valarray<CFreal> > jacobDets(2,std::valarray<CFreal>(m_nbrFlxPnts[m_orient]));
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    // get neighbouring cell
    GeometricEntity* cell = m_face->getNeighborGeo(iSide);

    // compute volume
    m_cellVolumes[iSide] = cell->computeVolume();

    // compute Jacobian determinants
    jacobDets[iSide] = cell->computeGeometricShapeFunctionJacobianDeterminant((*m_faceFlxPntCellMappedCoords)[m_orient][iSide]);
  }

  // compute inverse characteristic lengths
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_faceInvCharLengths[iFlx] = m_faceJacobVecAbsSizeFlxPnts[iFlx]/(jacobDets[LEFT][iFlx] + jacobDets[RIGHT][iFlx]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsStates(const vector< vector< State* >* >& cellStates, bool onlyExtraVars)
{
  cf_assert(cellStates.size() == 2);
  if (!onlyExtraVars)
  {
    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      m_statesReconstr->reconstructStates(*cellStates[iSide],m_flxPntSol[iSide],*m_flxPntsRecCoefs,
                                          (*m_faceFlxPntConn)[m_orient][iSide],
                                          *m_flxPntMatrixIdxForReconstruction,
                                          *m_solPntIdxsForReconstruction);
    }
  }

  // if needed, reconstruct the extra variables
  if (m_nbrExtraVars > 0)
  {
    cf_assert(socket_extraVars.isConnected());
    DataHandle<RealVector> extraVars = socket_extraVars.getDataHandle();

    for (CFuint iSide = 0; iSide < 2; ++iSide)
    {
      // get extra vars at the solution points in this cell
      const CFuint nbrSolPnts = cellStates[iSide]->size();
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        const CFuint stateID = (*cellStates[iSide])[iSol]->getLocalID();
        m_cellExtraVars[iSide][iSol] = &extraVars[stateID];
      }

      // reconstruct extra variables
      m_statesReconstr->reconstructExtraVars(m_cellExtraVars[iSide],m_flxPntExtraVars[iSide],*m_flxPntsRecCoefs,
                                            (*m_faceFlxPntConn)[m_orient][iSide],
                                              *m_flxPntMatrixIdxForReconstruction,
                                              *m_solPntIdxsForReconstruction);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsGradients(const vector< vector< vector< RealVector >* > >& cellGrads)
{
  cf_assert(cellGrads.size() == 2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_statesReconstr->reconstructGradients(cellGrads[iSide],m_flxPntGrads[iSide],*m_flxPntsRecCoefs,
                                           (*m_faceFlxPntConn)[m_orient][iSide],
                                           *m_flxPntMatrixIdxForReconstruction,
                                           *m_solPntIdxsForReconstruction);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::reconstructFluxPntsGradients(const CFuint side,
                                                        const vector< vector< RealVector >* >& cellGrads)
{
  m_statesReconstr->reconstructGradients(cellGrads,m_flxPntGrads[side],*m_flxPntsRecCoefs,
                                         (*m_faceFlxPntConn)[m_orient][side],
                                         *m_flxPntMatrixIdxForReconstruction,
                                         *m_solPntIdxsForReconstruction);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::backupAndReconstructPhysVar(const CFuint side,
                                                       const CFuint iVar,
                                                       const vector< State* >& cellStates)
{
  // backup
  backupPhysVar(side,iVar);

  // reconstruct
  m_statesReconstr
      ->reconstructPhysVar(iVar,cellStates,m_flxPntSol[side],*m_flxPntsRecCoefs,
                           (*m_faceFlxPntConn)[m_orient][side],
                           *m_flxPntMatrixIdxForReconstruction,
                           *m_solPntIdxsForReconstruction);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::backupPhysVar(const CFuint side, const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_backupPhysVar[iFlx] = (*m_flxPntSol[side][iFlx])[iVar];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::restorePhysVar(const CFuint side, const CFuint iVar)
{
  cf_assert(m_nbrFlxPnts[m_orient] <= m_backupPhysVar.size());
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    (*m_flxPntSol[side][iFlx])[iVar] = m_backupPhysVar[iFlx];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeConvFaceTerm(vector< RealVector>& resUpdates)
{
  // evaluate the riemann fluxes in all the flux points
  m_flxPntRiemannFlux = m_riemannFluxComputer->computeFlux(m_flxPntSol[LEFT],m_flxPntSol[RIGHT],
                                                           m_unitNormalFlxPnts,m_nbrFlxPnts[m_orient]);

  // compute the actual face term
  computeFaceTermFromFlxPntFluxes(resUpdates);

  // change the sign of the updates to the residuals
  resUpdates[LEFT ] *= -1.0;
  resUpdates[RIGHT] *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeConvFaceTermAndWaveSpeedUpdates(vector< RealVector>& resUpdates,
                                                                  vector< CFreal >& waveSpeedUpd)
{
  // compute the face term
  computeConvFaceTerm(resUpdates);

// cout << m_flxPntSol.size() << "\n" << flush;
// cout << m_flxPntSol[0].size() << "\n" << flush;


  // compute the wave speed updates for the neighbouring cells
  cf_assert(waveSpeedUpd.size() == 2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    waveSpeedUpd[iSide] = 0.0;
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      const CFreal jacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
                                   (*m_faceIntegrationCoefs)[iFlx];
      // transform update states to physical data to calculate eigenvalues

// cout << m_flxPntSol[iSide][iFlx]->getLocalID() << " " << flush;

      m_updateVarSet->computePhysicalData(*m_flxPntSol[iSide][iFlx], m_pData);
      waveSpeedUpd[iSide] += jacobXIntCoef*
          m_updateVarSet->getMaxAbsEigenValue(m_pData,m_unitNormalFlxPnts[iFlx]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeGradientFaceTerm(vector< vector< vector< RealVector > > >& gradUpdates)
{
//   // compute the averaged gradient variables in the flux points (stored in m_flxPntRiemannFlux)
//   // (m_leftRVSol and m_rightRVSol are RealVector pointers to the left and right states)
//   m_flxPntRiemannFlux = m_faceDiffFluxComputer->computeAvgGradVars(m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
//                                                                    m_nbrFlxPnts[m_orient]);
// 
//   // compute the actual volume term contribution to the gradient
//   computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeGradientExtraVarsFaceTerm
                                                (vector< vector< vector< RealVector > > >& gradUpdates)
{
  // compute the averaged extra variables in the flux points (stored in m_flxPntRiemannFlux)
  for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
  {
    m_flxPntRiemannFlux[iFlx]
      = 0.5*(*m_flxPntExtraVars[LEFT][iFlx] + *m_flxPntExtraVars[RIGHT][iFlx]);
  }

  // compute the actual volume term contribution to the gradient
  computeGradFaceTermFromFlxPntSol(gradUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeDiffFaceTerm(vector< RealVector>& resUpdates)
{
//   // compute the diffusive fluxes in the face flux points
//   m_flxPntRiemannFlux = m_faceDiffFluxComputer->computeDiffFlux(m_flxPntGradPtrs[LEFT],m_flxPntGradPtrs[RIGHT],
//                                                                 m_flxPntRVSol[LEFT],m_flxPntRVSol[RIGHT],
//                                                                 m_faceInvCharLengths,m_unitNormalFlxPnts,
//                                                                 m_nbrFlxPnts[m_orient]);
// 
//   // compute the actual face term
//   computeFaceTermFromFlxPntFluxes(resUpdates);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeDiffFaceTermAndUpdateCoefContributions(vector< RealVector>& resUpdates,
                                                                         vector< CFreal >& updateCoefContr)
{
//   // compute the face term
//   computeDiffFaceTerm(resUpdates);
// 
//   // compute the update coefficient contributions for the neighbouring cells
//   const CFreal visc = 1.0;
//   for (CFuint iSide = 0; iSide < 2; ++iSide)
//   {
//     updateCoefContr[iSide] = 0.0;
//     for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
//     {
//       const CFreal jacobXJacobXIntCoef = m_faceJacobVecAbsSizeFlxPnts[iFlx]*
//                                          m_faceJacobVecAbsSizeFlxPnts[iFlx]*
//                                          (*m_faceIntegrationCoefs)[iFlx]*m_cflConvDiffRatio;
//       updateCoefContr[iSide] += visc*jacobXJacobXIntCoef/m_cellVolumes[iSide];
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeFaceTermFromFlxPntFluxes(vector< RealVector>& resUpdates)
{
  cf_assert(resUpdates.size() == 2);

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over sides
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
  // set updates to zero
    resUpdates[iSide] = 0.0;

    // loop over flux points
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      // flux point index
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iSide][iFlx];

      // flux point index in the matrix m_solPntsDerivCoefs
      const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

      // evaluate flux projected on the projection vector
      const RealVector fluxXProjVect = m_flxPntRiemannFlux[iFlx]*
                                       m_faceJacobVecSizeFlxPnts[iFlx][iSide];

      // add contribution of this flux point to the solution points
      for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
      {
        // first residual in solution point index
        CFuint resIdx = m_nbrEqs*(*m_solPntIdxsForDerivation)[flxIdx][iSol];
        for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar, ++resIdx)
        {
          resUpdates[iSide][resIdx] +=
              (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*fluxXProjVect[iVar];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::computeGradFaceTermFromFlxPntSol(vector< vector< vector< RealVector > > >& gradUpdates)
{
  cf_assert(gradUpdates.size() == 2);

  // number of solution points one flux point contributes to
  cf_assert(m_solPntIdxsForDerivation->size() > 0);
  const CFuint nbrSolsPerFlx = (*m_solPntIdxsForDerivation)[0].size();

  // loop over sides
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
//     CF_DEBUG_OBJ(iSide);
    // set updates to zero
    const CFuint nbrSolPnts = gradUpdates[iSide].size();
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        gradUpdates[iSide][iSol][iGrad] = 0.0;
      }
    }

    // loop over flux points
    for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts[m_orient]; ++iFlx)
    {
      // flux point index
      const CFuint flxIdx = (*m_faceFlxPntConn)[m_orient][iSide][iFlx];

      // flux point index in the matrix m_solPntsDerivCoefs
      const CFuint flxPntMatrixIdx = (*m_flxPntMatrixIdxForDerivation)[flxIdx];

      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        // compute gradient term
        m_gradTerm = m_flxPntRiemannFlux[iFlx][iGrad]*
                     m_faceJacobVecSizeFlxPnts[iFlx][iSide]*
                     m_unitNormalFlxPnts[iFlx];
//         CF_DEBUG_OBJ(m_faceJacobVecSizeFlxPnts[iFlx][iSide]);
//         CF_DEBUG_OBJ(m_unitNormalFlxPnts[iFlx]);

        // add contribution of this flux point to the solution points
        for (CFuint iSol = 0; iSol < nbrSolsPerFlx; ++iSol)
        {
          // solution point index
          const CFuint solIdx = (*m_solPntIdxsForDerivation)[flxIdx][iSol];
          gradUpdates[iSide][solIdx][iGrad] +=
              (*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx)*m_gradTerm;
//           CF_DEBUG_OBJ((*m_solPntsDerivCoefs)(iSol,flxPntMatrixIdx));
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  FluxReconstructionSolverStrategy::setup();

  // m_dimensionality and number of variables
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVar();

  // get number of extra variables from the update variable set
  m_nbrExtraVars = m_updateVarSet->getExtraPhysicalVarsSize();
  if (m_nbrExtraVars > 0 && !socket_extraVars.isConnected())
  {
    throw ConsistencyException(FromHere(), "Extra variables socket not connected but required...");
  }

  // get the states reconstructor
  m_statesReconstr = getMethodData().getStatesReconstructor();

//   // get the Riemann flux
//   m_riemannFluxComputer = getMethodData().getRiemannFlux();
// 
//   // get the face diffusive flux
//   m_faceDiffFluxComputer = getMethodData().getFaceDiffusiveFlux();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();
  cf_assert(sdLocalData.size() > 0);

  // number of solution points
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // resize variables
  m_cellVolumes  .resize(2);

  // resize m_cellExtraVars
  m_cellExtraVars.resize(2);
  m_cellExtraVars[LEFT ].resize(nbrSolPnts);
  m_cellExtraVars[RIGHT].resize(nbrSolPnts);

  // number of flux points
  const CFuint nbrFlxPnts = sdLocalData[0]->getNbrOfFaceFlxPnts();

  // resize m_faceJacobVecAbsSizeFlxPnts
  m_faceJacobVecAbsSizeFlxPnts.resize(nbrFlxPnts);

  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(nbrFlxPnts,vector<CFreal>(2));

  // resize m_faceInvCharLengths
  m_faceInvCharLengths.resize(nbrFlxPnts);

  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(nbrFlxPnts,RealVector(m_dim));

  // resize m_flxPntRiemannFlux
  m_flxPntRiemannFlux.resize(nbrFlxPnts,RealVector(m_nbrEqs));

  // create left and right states and extra variables
  m_flxPntSol.resize(2);
  m_flxPntExtraVars.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_flxPntSol[iSide]      .push_back(new State()                   );
      m_flxPntExtraVars[iSide].push_back(new RealVector(m_nbrExtraVars));
    }
  }

  // resize m_backupPhysVar
  m_backupPhysVar.resize(nbrFlxPnts);

  // set m_allSol pointers to left and right states
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    m_allSol      .push_back(m_flxPntSol      [LEFT ][iFlx]);
    m_allSol      .push_back(m_flxPntSol      [RIGHT][iFlx]);
    m_allExtraVars.push_back(m_flxPntExtraVars[LEFT ][iFlx]);
    m_allExtraVars.push_back(m_flxPntExtraVars[RIGHT][iFlx]);
  }

  // setup variables for gradient computation
  if (false) //getMethodData().hasDiffTerm())
  {
    // get the diffusive varset
    m_diffusiveVarSet = getMethodData().getDiffusiveVar();
  }

  // set RealVector pointers to states
  m_flxPntRVSol.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_flxPntRVSol[iSide].push_back(m_flxPntSol[iSide][iFlx]);
    }
  }

  // create gradients for flux points
  m_flxPntGrads.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_flxPntGrads[iSide].resize(nbrFlxPnts);
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_flxPntGrads[iSide][iFlx].resize(m_nbrEqs);
      for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
      {
        m_flxPntGrads[iSide][iFlx][iGrad] = new RealVector(m_dim);
      }
    }
  }

  // create pointers to gradients for flux points
  m_flxPntGradPtrs.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
    {
      m_flxPntGradPtrs[iSide].push_back(&m_flxPntGrads[iSide][iFlx]);
    }
  }

  // resize m_gradTerm
  m_gradTerm.resize(m_dim);

  // set maximum number of flux points that will be passed at one time to the physical model
  const CFuint maxNbrFlxPnts = 2*nbrFlxPnts;
  const CFuint prevMaxNbrStatesData = getMethodData().getMaxNbrStatesData();
  getMethodData().
      setMaxNbrStatesData(prevMaxNbrStatesData > maxNbrFlxPnts ? prevMaxNbrStatesData : maxNbrFlxPnts);

  // set maximum number of points in which the Riemann flux has to be evaluated at the same time
  const CFuint prevMaxNbrRFluxPnts = getMethodData().getMaxNbrRFluxPnts();
  getMethodData().
      setMaxNbrRFluxPnts
      (
       prevMaxNbrRFluxPnts > nbrFlxPnts ? prevMaxNbrRFluxPnts : nbrFlxPnts
      );
  
  // resize the physical data temporary vector
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  convTerm->resizePhysicalData(m_pData);
}

//////////////////////////////////////////////////////////////////////////////

void BaseFaceTermComputer::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iSide = 0; iSide < m_flxPntSol.size(); ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < m_flxPntSol[iSide].size(); ++iFlx)
    {
      deletePtr(m_flxPntSol[iSide][iFlx]);
      deletePtr(m_flxPntExtraVars[iSide][iFlx]);
    }
    m_flxPntSol[iSide].resize(0);
    m_flxPntExtraVars[iSide].resize(0);
  }
  m_flxPntSol.resize(0);
  m_flxPntExtraVars.resize(0);

  for (CFuint iSide = 0; iSide < m_flxPntGrads.size(); ++iSide)
  {
    for (CFuint iFlx = 0; iFlx < m_flxPntGrads[iSide].size(); ++iFlx)
    {
      for (CFuint iGrad = 0; iGrad < m_flxPntGrads[iSide][iFlx].size(); ++iGrad)
      {
        deletePtr(m_flxPntGrads[iSide][iFlx][iGrad]);
      }
      m_flxPntGrads[iSide][iFlx].resize(0);
    }
    m_flxPntGrads[iSide].resize(0);
  }
  m_flxPntGrads.resize(0);

  // call unsetup of parent class
  FluxReconstructionSolverStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
    BaseFaceTermComputer::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_extraVars                  );
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD
