#include "AeroCoef/AeroCoefSpectralFD.hh"
#include "AeroCoef/NavierStokesConsComputeAeroSpectralFDCompact.hh"

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

#include "SpectralFD/CompactBndFaceTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
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

MethodCommandProvider<NavierStokesConsComputeAeroSpectralFDCompact, DataProcessingData, AeroCoefSpectralFDModule>
    NavierStokesConsComputeAeroSpectralFDCompactProvider("NavierStokesConsComputeAeroSpectralFDCompact");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::defineConfigOptions(Config::OptionList& options)
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

NavierStokesConsComputeAeroSpectralFDCompact::NavierStokesConsComputeAeroSpectralFDCompact(const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  m_faceBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_bcStateComputer(CFNULL),
  m_updateVarSet(CFNULL),
  m_diffusiveVarSet(CFNULL),
  m_lift(),
  m_drag(),
  m_fricDrag(),
  m_forceCoef(),
  m_fricForceCoef(),
  m_projSurf(),
  m_nbrEqs(),
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

NavierStokesConsComputeAeroSpectralFDCompact::~NavierStokesConsComputeAeroSpectralFDCompact()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesConsComputeAeroSpectralFDCompact::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NavierStokesConsComputeAeroSpectralFDCompact::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::executeOnTrs()
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

void NavierStokesConsComputeAeroSpectralFDCompact::computeFaceData()
{
  // set the face in the boundary face term computer
  m_bndFaceTermComputer->setCurrentFace(m_face);

  // compute the face data in the boundary face term computer
  m_bndFaceTermComputer->computeFacePntSetData();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::computeWall()
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

  // SET THE BOUNDARY CONDITION STATE COMPUTER
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

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

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // COMPUTE STATES AND GRADIENTS IN OUTPUT POINTS
          const vector< RealVector > outStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);
          const vector< vector< RealVector > > outGrads = m_bndFaceTermComputer->reconstructGivenPntsGrads(*m_cellStates);

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

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::computeAero()
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

  // SET THE BOUNDARY CONDITION STATE COMPUTER
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

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
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // GET FACE JACOBIAN VECTORS IN THE QUADRATURE POINTS
          const vector<RealVector>& jacobVecPntSet = *m_bndFaceTermComputer->getFaceJacobPntSet();

          // COMPUTE STATES AND GRADIENTS IN QUADRATURE POINTS
          const vector< RealVector > pntStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);
          const vector< vector< RealVector > > pntGrads = m_bndFaceTermComputer->reconstructGivenPntsGrads(*m_cellStates);

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

        // release the face
        m_faceBuilder->releaseGE();
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

void NavierStokesConsComputeAeroSpectralFDCompact::setup()
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
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get Euler var set
  m_updateVarSet = spectralFDData->getUpdateVar().d_castTo<EulerVarSet>();

  // get Navier-Stokes varset
  m_diffusiveVarSet = spectralFDData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();

  // get face builder
  m_faceBuilder = spectralFDData->getFaceBuilder();

  // get boundary face term computer
  m_bndFaceTermComputer = spectralFDData->getBndFaceTermComputer().d_castTo<CompactBndFaceTermComputer>();

  // get current TRS name
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  const std::string currTRSName = faceTrs->getName();

  // get corresponding BCStateComputer
  SafePtr< vector< SafePtr< BCStateComputer > > >
      bcStateComputers = spectralFDData->getBCStateComputers();
  const CFuint nbrBCs = bcStateComputers->size();
  bool bcFound = false;
  for (CFuint iBC = 0; iBC < nbrBCs && !bcFound; ++iBC)
  {
    SafePtr< vector< std::string > > bcTRSNames = (*bcStateComputers)[iBC]->getTRSNames();
    const CFuint nbrTRSs = bcTRSNames->size();
    for (CFuint iTRS = 0; iTRS < nbrTRSs; ++iTRS)
    {
      if (currTRSName == (*bcTRSNames)[iTRS])
      {
        m_bcStateComputer = (*bcStateComputers)[iBC];
        bcFound = true;
        break;
      }
    }
  }

  // dimensionality
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

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

  // allocate m_grads
  m_grads.resize(0);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_grads.push_back(new RealVector(dim));
  }

  prepareOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::unsetup()
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    deletePtr(m_grads[iEq]);
  }
  m_grads.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAeroSpectralFDCompact::prepareOutputFileWall()
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

void NavierStokesConsComputeAeroSpectralFDCompact::prepareOutputFileAero()
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

void NavierStokesConsComputeAeroSpectralFDCompact::updateOutputFileAero()
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

