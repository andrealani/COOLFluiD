#include "AeroCoef/AeroCoefSpectralFD.hh"
#include "AeroCoef/NavierStokesConsComputeAeroSpectralFD.hh"

#include "Common/PE.hh"

#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/ParserException.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/BaseVolTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/CellToFaceGEBuilder.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/SpectralFDMethod.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/TensorProductGaussIntegrator.hh"

#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::SpectralFD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesConsComputeAeroSpectralFD, DataProcessingData, AeroCoefSpectralFDModule>
    NavierStokesConsComputeAeroSpectralFDProvider("NavierStokesConsComputeAeroSpectralFD");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coeficients.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFreal >("uInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("rhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("pInf","Pressure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
   options.addConfigOption< vector< CFreal > >("FlowDir","Flow direction, used to determine orientation of skin friction.");
   options.addConfigOption< CFreal >("ProjSurf","Projected surface, needed for computation of drag coefficient");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesConsComputeAeroSpectralFD::NavierStokesConsComputeAeroSpectralFD(const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  m_spectralFDData(CFNULL),
  m_faceBuilder(CFNULL),
  m_cellBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_volTermComputer(CFNULL),
  m_faceTermComputers(),
  m_bndFaceTermComputers(),
  m_bcStateComputer(CFNULL),
  m_bcStateComputers(CFNULL),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_face(CFNULL),
  m_intCell(CFNULL),
  m_faces(CFNULL),
  m_faceNghbrStates(),
  m_solPntsLocalCoords(CFNULL),
  m_cellGrads(),
  m_gradUpdates(),
  m_solJacobDet(),
  m_otherFaceLocalIdxs(),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_nbrEqs(),
  m_lift(),
  m_drag(),
  m_fricDrag(),
  m_forceCoef(),
  m_fricForceCoef(),
  m_projSurf(),
  m_grads()
{
  m_fileWall = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFileWall = "Wall";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_saveRateWall = 1;
  setParameter("SaveRateWall",&m_saveRateWall);

  m_nameOutputFileAero = "AeroCoef.plt";
  setParameter("OutputFileAero",&m_nameOutputFileAero);

  m_saveRateAero = 1;
  setParameter("SaveRateAero",&m_saveRateAero);

  m_uInf = 0.;
  setParameter("uInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("rhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("pInf",&m_pInf);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

  m_flowDir = vector<CFreal>();
  setParameter("FlowDir",&m_flowDir);

  m_projSurf = 1.0;
  setParameter("ProjSurf",&m_projSurf);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesConsComputeAeroSpectralFD::~NavierStokesConsComputeAeroSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesConsComputeAeroSpectralFD::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NavierStokesConsComputeAeroSpectralFD::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % m_saveRateWall))
  {
      computeWall();
  }

  if(!(iter % m_saveRateAero))
  {
    computeAero();
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::setFaceAndVolumeTermData()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = m_spectralFDData->getSDLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[0]->getSolPntsLocalCoords();

  // set the BC state computer in the current boundary face term computer
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // set the data in the current boundary face term computer
  m_bndFaceTermComputer->setFaceTermData();

  // set the data in the face term computers
  for (CFuint iFace = 0; iFace < m_faceTermComputers.size(); ++iFace)
  {
    m_faceTermComputers[iFace]->setFaceTermData();
  }

  // set the data in the boundary face term computers
  for (CFuint iFace = 0; iFace < m_bndFaceTermComputers.size(); ++iFace)
  {
    m_bndFaceTermComputers[iFace]->setFaceTermData();
  }

  // set the data in the volume term computer
  m_volTermComputer->setVolumeTermData(0);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::setOtherFacesLocalIdxs()
{
  // get face ID of current face
  const CFuint currFaceID = m_face->getID();

  // loop over the faces of the cells
  const CFuint nbrFaces = m_intCell->nbNeighborGeos();
  cf_assert(m_faces->size() == nbrFaces);
  CFuint iFace = 0;
  for (CFuint faceIdx = 0; faceIdx < nbrFaces; ++faceIdx)
  {
    if ((*m_faces)[faceIdx]->getID() != currFaceID)
    {
      m_otherFaceLocalIdxs[iFace] = faceIdx;
      ++iFace;
    }
  }
  cf_assert(iFace == m_otherFaceLocalIdxs.size());
  cf_assert(iFace == nbrFaces-1);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::setFaceNeighbourStates()
{
  // get neighbour states of other faces
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];

    // get neigbouring states
    const CFuint nbrFaceNeighbours = (*m_isFaceOnBoundary)[faceIdx] ? 1 : 2;
    for (CFuint iSide = 0; iSide < nbrFaceNeighbours; ++iSide)
    {
      GeometricEntity* cell = (*m_faces)[faceIdx]->getNeighborGeo(iSide);
      m_faceNghbrStates[iFace][iSide] = cell->getStates();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::computeFaceData()
{
  // set the face in the boundary face term computer
  m_bndFaceTermComputer->setCurrentFace(m_face);

  // compute the face data in the boundary face term computer
  m_bndFaceTermComputer->computeFaceData();

  // compute the neighbouring cell data in the boundary face term computer
  m_bndFaceTermComputer->computeNeighbourCellData();// probably not necessary

  // reconstruct the states in the face flux points
  m_bndFaceTermComputer->reconstructFluxPntsStates(*m_cellStates);

  // compute the face data in the boundary face term computer
  m_bndFaceTermComputer->computeFacePntSetData();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::setCellData()
{
  // compute solution points Jacobian determinants
  m_solJacobDet = m_intCell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  // set the neighbouring cell and compute neighbouring cell data in the volume term computer
  m_volTermComputer->setCurrentCell(m_intCell);
  m_volTermComputer->computeCellData();
  m_volTermComputer->reconstructStates(*m_cellStates);

  // set orientation or boundary condition state computers of the other faces
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];

    // set orientation or boundary condition
    if ((*m_isFaceOnBoundary)[faceIdx])
    {
      // set BC state computer
      m_bndFaceTermComputers[iFace]
          ->setBcStateComputer((*m_bcStateComputers)[(*m_faceBCIdx)[faceIdx]]);

      // set the orientation of the face
      m_bndFaceTermComputers[iFace]->setFaceOrientation(faceIdx);

      // set the face in the boundary face term computer
      m_bndFaceTermComputers[iFace]->setCurrentFace((*m_faces)[faceIdx]);

      // compute the face data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeFaceData();

      // compute the neighbouring cell data in the boundary face term computer
      m_bndFaceTermComputers[iFace]->computeNeighbourCellData();// probably not necessary

      // reconstruct the states in the face flux points
      m_bndFaceTermComputers[iFace]->reconstructFluxPntsStates(*m_cellStates);
    }
    else
    {
      // set the orientation of the face
      m_faceTermComputers[iFace]->setFaceOrientation((*m_faceOrients)[faceIdx]);

      // set the face in the face term computer
      m_faceTermComputers[iFace]->setCurrentFace((*m_faces)[faceIdx]);

      // compute the face data in the face term computer
      m_faceTermComputers[iFace]->computeFaceData();

      // compute the neighbouring cell data in the face term computer
      m_faceTermComputers[iFace]->computeNeighbourCellData();// probably not necessary

      // reconstruct the states in the face flux points
      m_faceTermComputers[iFace]->
          reconstructFluxPntsStates(m_faceNghbrStates[iFace]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::computeCellGradients()
{
  // compute volume contribution to current cell gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates[0]);
  const CFuint nbrSolPnts = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq] = m_gradUpdates[0][iSol][iEq];
    }
  }

  // compute face contributions to the gradients
  const CFuint nbrOtherFaces = m_otherFaceLocalIdxs.size();
  for (CFuint iFace = 0; iFace < nbrOtherFaces; ++iFace)
  {
    // get local face index
    const CFuint faceIdx = m_otherFaceLocalIdxs[iFace];
    if ((*m_isFaceOnBoundary)[faceIdx])
    {
      // compute the boundary face contribution to the gradients
      m_bndFaceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates[0]);

      // add the contribution to the gradients
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
        }
      }
    }
    else
    {
      // compute the face contribution to the gradients
      m_faceTermComputers[iFace]->computeGradientFaceTerm(m_gradUpdates);

      // add the contribution to the gradients
      const CFuint side = (*m_currCellSide)[faceIdx];
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*m_cellGrads[iSol])[iEq] += m_gradUpdates[side][iSol][iEq];
        }
      }
    }
  }

  // compute the contribution of the current boundary face to the gradients
  m_bndFaceTermComputer->computeGradientFaceTerm(m_gradUpdates[0]);

  // add the contribution to the gradients
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq] += m_gradUpdates[0][iSol][iEq];
    }
  }

  // divide by solution point Jacobian determinant
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    const CFreal invJacobDet = 1.0/m_solJacobDet[iSol];
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file handle is opened here

  // get some data needed further
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsM1 = nbEqs - 1;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get SpectralFDMethodData
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get bndFacesStartIdxs from SpectralFDMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = spectralFDData->getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;

  // SET THE VOLUME AND FACE TERM COMPUTER DATA
  setFaceAndVolumeTermData();

  // GET FACE MAPPED COORDINATES OF OUTPUT POINTS
  vector< SpectralFDElementData* >& sdLocalData = spectralFDData->getSDLocalData();
  SafePtr< vector< RealVector > > outputPntsMappedCoord = sdLocalData[0]->getFaceOutputPntFaceMappedCoords();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
      m_bndFaceTermComputer->setFaceOrientation(m_orient);
      m_bndFaceTermComputer->setPointSet(*outputPntsMappedCoord);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // GET THE NEIGHBOURING CELL
        m_intCell = m_face->getNeighborGeo(0);

        // GET THE STATES IN THE NEIGHBOURING CELL
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // BUILD THE CELL WITH CONNECTIVITY TO ITS FACES
          geoDataCB.idx = m_intCell->getID();
          m_intCell = m_cellBuilder->buildGE();

          // GET ALL THE FACES NEIGHBOURING THE CELLS
          m_faces = m_intCell->getNeighborGeos();

          // SET THE LOCAL INDEXES OF THE OTHER FACES THAN THE CURRENT FACES
          setOtherFacesLocalIdxs();

          // GET THE NEIGBOURING STATES OF THE OTHER FACES
          setFaceNeighbourStates();

          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // SET THE INTERNAL CELL DATA
          setCellData();

          // COMPUTE INTERNAL CELL GRADIENTS
          computeCellGradients();

          // COMPUTE STATES AND GRADIENTS IN OUTPUT POINTS
          const vector< RealVector > outStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);
          const vector< vector< RealVector > > outGrads = m_bndFaceTermComputer->reconstructGivenPntsGrads(m_cellGrads);

          // GET UNIT NORMALS
          const vector< RealVector >& outUnitNormals = *m_bndFaceTermComputer->getUnitNormalPntSet();

          // LOOP OVER OUTPUT POINTS
          const CFuint nbrOutPnts = outStates.size();
          cf_assert(nbrOutPnts == outputPntsMappedCoord->size());
          for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
          {
            // compute coordinates of output point
            const RealVector coord = m_face->computeCoordFromMappedCoord((*outputPntsMappedCoord)[iPnt]);

            // dereference current state
            const RealVector& state = outStates[iPnt];

            // dereference current unit normal
            const RealVector& normal = outUnitNormals[iPnt];

            // get current gradients
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
            {
              *m_grads[iEq] = outGrads[iPnt][iEq];
            }

            // pressure
            const CFreal rho = state[0];
            const CFreal invRho = 1./rho;
            CFreal rhoK2 = 0.0;
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              rhoK2 += state[iDim+1]*state[iDim+1];
            }
            rhoK2 *= 0.5*invRho;
            const CFreal p = gammaMinus1*(state[nbEqsM1] - rhoK2);

            // temperature
            const CFreal T = p*invRho/R;

            // dimensional values
            const CFreal TDim   = T   * m_updateVarSet->getModel()->getTempRef();
            const CFreal pDim   = p   * m_updateVarSet->getModel()->getPressRef();
            const CFreal rhoDim = rho * (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];

            // pressure coefficient
            const CFreal Cp = (pDim - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

            // compute dynamic viscosity
            m_diffusiveVarSet->setWallDistance(0.); // we are at the wall
            const CFreal mu  = m_diffusiveVarSet->getDynViscosity(state,m_grads);

            // compute the friction at the wall
            const CFreal refU = m_updateVarSet->getModel()->getVelRef();
            const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
            const CFreal refMu = (m_diffusiveVarSet->getModel().getReferencePhysicalData())[NSTerm::MU];
            const CFreal refTau = refMu*refU/refLength;
            const bool adim = Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional();

            // compute shear stress and skin friction tensors
            RealMatrix tauTensor(dim,dim);
            RealMatrix skinFrictionTensor(dim,dim);
            RealVector frictionForce(dim);
            CFreal skinFriction = 0.0;
            if (dim == DIM_2D)
            {
              const RealVector& gradU = *m_grads[1];
              const RealVector& gradV = *m_grads[2];
              const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY]);

              tauTensor(XX,XX) = mu*(2.0*gradU[XX] - divU);
              tauTensor(YY,YY) = mu*(2.0*gradV[YY] - divU);

              tauTensor(XX,YY) = tauTensor(YY,XX) = mu*(gradU[YY] + gradV[XX]);

              if(adim)
              {
                const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
                const CFreal Re = m_diffusiveVarSet->getModel().getReynolds();

                /// @todo check this
                skinFrictionTensor = tauTensor/(0.5*machInf*sqrt(gamma)/Re);
              }
              else
              {
                skinFrictionTensor = tauTensor/(0.5*m_rhoInf*m_uInf*m_uInf);
              }

              // compute the viscous force on this wall face.
              frictionForce[XX] = refTau*(tauTensor(XX,XX)*normal[XX] + tauTensor(XX,YY)*normal[YY]);
              frictionForce[YY] = refTau*(tauTensor(YY,XX)*normal[XX] + tauTensor(YY,YY)*normal[YY]);

              // compute adimensional skin friction
              RealVector skinFrictionVector(2);
              skinFrictionVector[XX] = skinFrictionTensor(XX,XX)*normal[XX] + skinFrictionTensor(XX,YY)*normal[YY];
              skinFrictionVector[YY] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY];
              if (m_flowDir[XX]*normal[YY] - m_flowDir[YY]*normal[XX] > 0.0)
              {
                // sign is reversed to obtain the force on the body
                skinFriction = - skinFrictionVector[XX]*normal[YY] + skinFrictionVector[YY]*normal[XX];
              }
              else
              {
                // sign is reversed to obtain the force on the body
                skinFriction = skinFrictionVector[XX]*normal[YY] - skinFrictionVector[YY]*normal[XX];
              }
            }
            else
            {
              cf_assert(dim == DIM_3D);

              const RealVector& gradU = *m_grads[1];
              const RealVector& gradV = *m_grads[2];
              const RealVector& gradW = *m_grads[3];
              const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY] + gradW[ZZ]);

              tauTensor(XX,XX) = mu*(2.0*gradU[XX] - divU);
              tauTensor(YY,YY) = mu*(2.0*gradV[YY] - divU);
              tauTensor(ZZ,ZZ) = mu*(2.0*gradW[ZZ] - divU);

              tauTensor(XX,YY) = tauTensor(YY,XX) = mu*(gradU[YY] + gradV[XX]);
              tauTensor(XX,ZZ) = tauTensor(ZZ,XX) = mu*(gradU[ZZ] + gradW[XX]);
              tauTensor(YY,ZZ) = tauTensor(ZZ,YY) = mu*(gradV[ZZ] + gradW[YY]);

              if(adim)
              {
                const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
                const CFreal Re = m_diffusiveVarSet->getModel().getReynolds();

                skinFrictionTensor = tauTensor/(0.5*machInf*sqrt(gamma)/Re);
              }
              else
              {
                skinFrictionTensor = tauTensor/(0.5*m_rhoInf*m_uInf*m_uInf);
              }

              //Compute the viscous force on this wall face.
              frictionForce[XX] = refTau*(tauTensor(XX,XX)*normal[XX] + tauTensor(XX,YY)*normal[YY] + tauTensor(XX,ZZ)*normal[ZZ]);
              frictionForce[YY] = refTau*(tauTensor(YY,XX)*normal[XX] + tauTensor(YY,YY)*normal[YY] + tauTensor(YY,ZZ)*normal[ZZ]);
              frictionForce[ZZ] = refTau*(tauTensor(ZZ,XX)*normal[XX] + tauTensor(ZZ,YY)*normal[YY] + tauTensor(ZZ,ZZ)*normal[ZZ]);

              // compute adimensional skin friction
              RealVector skinFrictionVector(3);
              skinFrictionVector[XX] = skinFrictionTensor(XX,XX)*normal[XX] + skinFrictionTensor(XX,YY)*normal[YY] + skinFrictionTensor(XX,ZZ)*normal[ZZ];
              skinFrictionVector[YY] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY] + skinFrictionTensor(YY,ZZ)*normal[ZZ];
              skinFrictionVector[ZZ] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY] + skinFrictionTensor(ZZ,ZZ)*normal[ZZ];
              skinFriction = skinFrictionVector.norm2();// like this, the sign of this coefficient is always positive. How to keep the sign?
            }

            // heat flux
            const CFreal heatFluxDim = m_diffusiveVarSet->getHeatFlux(state,m_grads,normal); // check dimensionality

            // stanton number
            const CFreal stanton = 0.0;
//             const CFreal stanton = heatFluxDim/(m_rhoInf*pow(m_uInf,3.0));

            // compute y+ value
            // not implemented
            const CFreal yPlus = 0.0;

            // Output to File
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              fout << coord[iDim]  << " ";
            }
            fout << pDim         << " "
                 << TDim         << " "
                 << rhoDim       << " "
                 << Cp           << " "
                 << heatFluxDim  << " "
                 << skinFriction << " "
                 << stanton      << " "
                 << yPlus        << " "
                 << mu           << "\n";
          }
        }

        // release the face and the cell
        m_faceBuilder->releaseGE();
        m_cellBuilder->releaseGE();
      }
    }
  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::computeAero()
{
  CFAUTOTRACE;

  // get some data needed further
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsM1 = nbEqs - 1;
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // initialize Cp, skin friction and surface
  RealVector Cp(dim);
  Cp = 0.;
  RealVector fric(dim);
  fric = 0;
  CFreal surface = 0.0;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get SpectralFDMethodData
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get bndFacesStartIdxs from SpectralFDMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = spectralFDData->getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // get the geodata of the cell builders and set the TRS
  CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
  geoDataCB.trs = cellTrs;

  // SET THE VOLUME AND FACE TERM COMPUTER DATA
  setFaceAndVolumeTermData();

  // GET FACE MAPPED COORDINATES OF OUTPUT POINTS
  vector< SpectralFDElementData* >& sdLocalData = spectralFDData->getSDLocalData();
  SafePtr< vector< vector< RealVector > > > faceNodeMappedCoords = sdLocalData[0]->getFaceNodeCoords();

  // create integrator object (assume P2 geometrical polynomials) for integration over faces
  // and get quadrature point face mapped coordinates and wheights
  const CFuint solPolyOrder = sdLocalData[0]->getPolyOrder();
  CFuint integratorOrder = 0;
  if (dim == 2)
  {
    integratorOrder = solPolyOrder+1;
  }
  else if (dim == 3)
  {
    integratorOrder = solPolyOrder+3;
  }
  TensorProductGaussIntegrator tpIntegrator(static_cast<CFDim>(dim-1),static_cast<CFPolyOrder::Type>(integratorOrder));
  vector< RealVector > qPntsFaceMapCoords = tpIntegrator.getQuadPntsMappedCoords();
  vector< CFreal >     qPntsWheights  = tpIntegrator.getQuadPntsWheights();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
      m_bndFaceTermComputer->setFaceOrientation(m_orient);
      m_bndFaceTermComputer->setPointSet(qPntsFaceMapCoords);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // BUILD THE CELL WITH CONNECTIVITY TO ITS FACES
          geoDataCB.idx = m_intCell->getID();
          m_intCell = m_cellBuilder->buildGE();

          // GET ALL THE FACES NEIGHBOURING THE CELLS
          m_faces = m_intCell->getNeighborGeos();

          // SET THE LOCAL INDEXES OF THE OTHER FACES THAN THE CURRENT FACES
          setOtherFacesLocalIdxs();

          // GET THE NEIGBOURING STATES OF THE OTHER FACES
          setFaceNeighbourStates();

          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // SET THE INTERNAL CELL DATA
          setCellData();

          // COMPUTE INTERNAL CELL GRADIENTS
          computeCellGradients();

          // COMPUTE STATES AND GRADIENTS IN OUTPUT POINTS
          const vector< RealVector > pntStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);
          const vector< vector< RealVector > > pntGrads = m_bndFaceTermComputer->reconstructGivenPntsGrads(m_cellGrads);

          // GET FACE JACOBIAN VECTORS IN THE QUADRATURE POINTS
          const vector<RealVector>& jacobVecPntSet = *m_bndFaceTermComputer->getFaceJacobPntSet();

          // COMPUTE INTEGRAL OF PRESSURE, SKIN FRICTION AND FACE SURFACE
          const CFuint nbrPnts = pntStates.size();
          RealVector cpFace(dim);
          cpFace = 0;
          RealVector fricFace(dim);
          fricFace = 0;
          CFreal surfFace = 0.0;
          cf_assert(nbrPnts == jacobVecPntSet.size());
          cf_assert(nbrPnts == qPntsWheights.size());
          for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
          {
            // compute coordinates of output point
            const RealVector coord = m_face->computeCoordFromMappedCoord(qPntsFaceMapCoords[iPnt]);

            // dereference current state
            const RealVector& state = pntStates[iPnt];

            // get current gradients
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
            {
              *m_grads[iEq] = pntGrads[iPnt][iEq];
            }

            // pressure
            const CFreal rho = state[0];
            const CFreal invRho = 1./rho;
            CFreal rhoK2 = 0.0;
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              rhoK2 += state[iDim+1]*state[iDim+1];
            }
            rhoK2 *= 0.5*invRho;
            const CFreal p = gammaMinus1*(state[nbEqsM1] - rhoK2);

            // compute shear stress
            m_diffusiveVarSet->setWallDistance(0.); // we are at the wall
            const CFreal mu  = m_diffusiveVarSet->getDynViscosity(state,m_grads);
            RealMatrix tauTensor(dim,dim);
            RealVector tauTensorXJacobVec(dim);
            if (dim == DIM_2D)
            {
              const RealVector& gradU = *m_grads[1];
              const RealVector& gradV = *m_grads[2];
              const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY]);

              tauTensor(XX,XX) = mu*(2.0*gradU[XX] - divU);
              tauTensor(YY,YY) = mu*(2.0*gradV[YY] - divU);

              tauTensor(XX,YY) = tauTensor(YY,XX) = mu*(gradU[YY] + gradV[XX]);

              tauTensorXJacobVec[XX] = tauTensor(XX,XX)*jacobVecPntSet[iPnt][XX] + tauTensor(XX,YY)*jacobVecPntSet[iPnt][YY];
              tauTensorXJacobVec[YY] = tauTensor(YY,XX)*jacobVecPntSet[iPnt][XX] + tauTensor(YY,YY)*jacobVecPntSet[iPnt][YY];
            }
            else
            {
              cf_assert(dim == DIM_3D);

              const RealVector& gradU = *m_grads[1];
              const RealVector& gradV = *m_grads[2];
              const RealVector& gradW = *m_grads[3];
              const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY] + gradW[ZZ]);

              tauTensor(XX,XX) = mu*(2.0*gradU[XX] - divU);
              tauTensor(YY,YY) = mu*(2.0*gradV[YY] - divU);
              tauTensor(ZZ,ZZ) = mu*(2.0*gradW[ZZ] - divU);

              tauTensor(XX,YY) = tauTensor(YY,XX) = mu*(gradU[YY] + gradV[XX]);
              tauTensor(XX,ZZ) = tauTensor(ZZ,XX) = mu*(gradU[ZZ] + gradW[XX]);
              tauTensor(YY,ZZ) = tauTensor(ZZ,YY) = mu*(gradV[ZZ] + gradW[YY]);

              tauTensorXJacobVec[XX] =
                  tauTensor(XX,XX)*jacobVecPntSet[iPnt][XX] + tauTensor(XX,YY)*jacobVecPntSet[iPnt][YY] + tauTensor(XX,ZZ)*jacobVecPntSet[iPnt][ZZ];
              tauTensorXJacobVec[YY] =
                  tauTensor(YY,XX)*jacobVecPntSet[iPnt][XX] + tauTensor(YY,YY)*jacobVecPntSet[iPnt][YY] + tauTensor(YY,ZZ)*jacobVecPntSet[iPnt][ZZ];
              tauTensorXJacobVec[ZZ] =
                  tauTensor(ZZ,XX)*jacobVecPntSet[iPnt][XX] + tauTensor(ZZ,YY)*jacobVecPntSet[iPnt][YY] + tauTensor(ZZ,ZZ)*jacobVecPntSet[iPnt][ZZ];
            }

            // add contribution to face pressure integral
            cpFace += qPntsWheights[iPnt]*jacobVecPntSet[iPnt]*p;

            // add contribution to face skin friction integral
            fricFace -= qPntsWheights[iPnt]*tauTensorXJacobVec;

            // add contribution to face surface
            surfFace += qPntsWheights[iPnt]*jacobVecPntSet[iPnt].norm2();
          }

          // add to global Cp, skin friction and surface
          Cp += cpFace;
          fric += fricFace;
          surface += surfFace;
        }

        // release the face and the cell
        m_faceBuilder->releaseGE();
        m_cellBuilder->releaseGE();
      }
    }
  }

  // adimensionalize Cp
  const CFreal pressRef = m_updateVarSet->getModel()->getPressRef();
  Cp *= pressRef/(0.5*m_rhoInf*m_uInf*m_uInf*m_projSurf);

  // adimensionalize skin friction
  const CFreal refU = m_updateVarSet->getModel()->getVelRef();
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  const CFreal refMu = (m_diffusiveVarSet->getModel().getReferencePhysicalData())[NSTerm::MU];
  const CFreal refTau = refMu*refU/refLength;
  fric *= refTau/(0.5*m_rhoInf*m_uInf*m_uInf*m_projSurf);

  // set variables
  m_fricForceCoef = fric;
  m_forceCoef = fric + Cp;

  // project m_forceCoef
  /// @todo implement the following in a rigorous way. The following is only for 2D for the lift
  m_drag     = 0.;
  m_fricDrag = 0.;
  for (CFuint iDim = 0; iDim < m_flowDir.size(); ++iDim)
  {
    m_drag     += m_flowDir[iDim]*m_forceCoef    [iDim];
    m_fricDrag += m_flowDir[iDim]*m_fricForceCoef[iDim];
  }
  m_lift = m_forceCoef[YY]*m_flowDir[XX] - m_forceCoef[XX]*m_flowDir[YY];

  //Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::setup()
{
  CFAUTOTRACE;

  // set lift and drag to zero
  m_lift = 0.;
  m_drag = 0.;

  // get SpectralFDMethodData
  // assume that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  m_spectralFDData = spectralFDMethod->getData();

  // get Euler var set
  m_updateVarSet = m_spectralFDData->getUpdateVar().d_castTo<EulerVarSet>();

  // get Navier-Stokes varset
  m_diffusiveVarSet = m_spectralFDData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();

  // get face builder
  m_faceBuilder = m_spectralFDData->getFaceBuilder();

  // get cell builder
  m_cellBuilder = m_spectralFDData->getCellBuilder();

  // get some additional data for cell building
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  // get boundary face term computer
  m_bndFaceTermComputer = m_spectralFDData->getBndFaceTermComputer();

  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = m_spectralFDData->getSDLocalData();
  cf_assert(sdLocalData.size() > 0);

  // get the number of faces in a cell
  const CFuint nbrFaces = sdLocalData[0]->getNbrCellFaces();
  const CFuint nbrFacesM1 = nbrFaces - 1;

  // get face term computers and additional boundary face term computers
  CFuint computerIdx = 0;
  m_faceTermComputers   .resize(nbrFacesM1);
  m_bndFaceTermComputers.resize(nbrFacesM1);
  for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace, ++computerIdx)
  {
    m_faceTermComputers   [iFace] = m_spectralFDData->getAdditionalFaceTermComputer   (computerIdx);
    m_bndFaceTermComputers[iFace] = m_spectralFDData->getAdditionalBndFaceTermComputer(computerIdx);
  }

  // get the volume term computer
  m_volTermComputer = m_spectralFDData->getVolTermComputer();

  // get current TRS name
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  const std::string currTRSName = faceTrs->getName();

  // get BCStateComputers and BCStateComputer corresponding to current TRS
  m_bcStateComputers = m_spectralFDData->getBCStateComputers();
  const CFuint nbrBCs = m_bcStateComputers->size();
  bool bcFound = false;
  for (CFuint iBC = 0; iBC < nbrBCs && !bcFound; ++iBC)
  {
    SafePtr< vector< std::string > > bcTRSNames = (*m_bcStateComputers)[iBC]->getTRSNames();
    const CFuint nbrTRSs = bcTRSNames->size();
    for (CFuint iTRS = 0; iTRS < nbrTRSs; ++iTRS)
    {
      if (currTRSName == (*bcTRSNames)[iTRS])
      {
        m_bcStateComputer = (*m_bcStateComputers)[iBC];
        bcFound = true;
        break;
      }
    }
  }

  // get the number of solution points in a cell
  const CFuint nbrSolPnts = sdLocalData[0]->getNbrOfSolPnts();

  // dimensionality and number of equations
  const CFuint dim    = PhysicalModelStack::getActive()->getDim ();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // resize m_otherFaceLocalIdxs
  m_otherFaceLocalIdxs.resize(nbrFacesM1);

  // resize m_faceNghbrStates
  m_faceNghbrStates.resize(nbrFacesM1);
  for (CFuint iFace = 0; iFace < nbrFacesM1; ++iFace)
  {
    m_faceNghbrStates[iFace].resize(2);
  }

  // resize gradient updates
  m_gradUpdates.resize(2);
  for (CFuint iSide = 0; iSide < 2; ++iSide)
  {
    m_gradUpdates[iSide].resize(nbrSolPnts);
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_gradUpdates[iSide][iSol].resize(m_nbrEqs);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_gradUpdates[iSide][iSol][iEq].resize(dim);
      }
    }
  }

  // resize m_solJacobDet
  m_solJacobDet.resize(nbrSolPnts);

  // allocate memory for m_cellGrads
  m_cellGrads.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_cellGrads[iSol] = new vector<RealVector>(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*m_cellGrads[iSol])[iEq].resize(dim);
    }
  }

  // allocate m_grads
  m_grads.resize(0);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_grads.push_back(new RealVector(dim));
  }

  // resize variables
  m_forceCoef    .resize(dim);
  m_fricForceCoef.resize(dim);

  // check if values at infinity are set
  if (!(m_uInf > 0. && m_rhoInf > 0. && m_pInf > 0.))
  {
    throw Common::BadValueException(FromHere(),"Values at infinity are not set!!!");
  }

  if (m_flowDir.size() != dim)
  {
    throw Common::BadValueException(FromHere(),"Flow direction vector has wrong dimensionality");
  }
  else
  {
    CFreal size = 0.;
    for (CFuint iDim = 0; iDim < dim; ++iDim)
    {
      size += m_flowDir[iDim]*m_flowDir[iDim];
    }
    size = sqrt(size);
    for (CFuint iDim = 0; iDim < dim; ++iDim)
    {
      m_flowDir[iDim] /= size;
    }
  }

  prepareOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::unsetup()
{
  for (CFuint iSol = 0; iSol < m_cellGrads.size(); ++iSol)
  {
    deletePtr(m_cellGrads[iSol]);
  }
  m_cellGrads.resize(0);

  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    deletePtr(m_grads[iEq]);
  }
  m_grads.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::prepareOutputFileWall()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_assert (!m_fileWall->isopen());

  using namespace boost::filesystem;
  path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
  file = PathAppender::getInstance().appendAllInfo( file );

  ofstream& fout = m_fileWall->open(file);

    fout << "TITLE = Wall_Values" << "\n";
  if (dim == 2)
  {
    fout << "VARIABLES = x y Pressure Temperature Density Cp heatF skinF Stanton yplus muwall" << "\n";
  }
  else
  {
    fout << "VARIABLES = x y z Pressure Temperature Density Cp heatF skinF Stanton yplus muwall" << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::prepareOutputFileAero()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

    convergenceFile << "TITLE = Aerodynamic_Coeficients"  << "\n";
  if (dim == 2)
  {
    convergenceFile << "VARIABLES = Iter PhysTime CL CD CDf Fx Fy Fx_f Fy_f" << "\n";
  }
  else
  {
    convergenceFile << "VARIABLES = Iter PhysTime CL CD CDf Fx Fy Fz Fx_f Fy_f Fz_f" << "\n";
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFD::updateOutputFileAero()
{
  CFAUTOTRACE;

  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  convergenceFile
  << subSysStatus->getNbIter()        << " "
  << subSysStatus->getCurrentTime()   << " "
  << m_lift                           << " "
  << m_drag                           << " "
  << m_fricDrag                       << " "
  << m_forceCoef                      << " "
  << m_fricForceCoef                  << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoeff

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

