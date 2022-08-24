#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BCSuperInletProjMHD.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Common/NotImplementedException.hh"

#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperInletProjMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCSuperInletProjMHDProvider("SuperInletProjMHD");

//////////////////////////////////////////////////////////////////////////////

BCSuperInletProjMHD::BCSuperInletProjMHD(const std::string& name) :
  BCStateComputer(name),
  m_initialSolutionMap(),
  m_faceBuilder(CFNULL),
  m_intCell(CFNULL),
  m_cellStates(CFNULL),
  m_currFace(CFNULL),
  m_orient(),
  m_cellStatesFlxPnt(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_faceFlxPntConn(CFNULL),
  m_thisTRS(),
  m_globalToLocalTRSFaceID()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_rhoBC = 1.0;
  setParameter("RhoBC",&m_rhoBC);
  
  m_pBC = 0.108;
  setParameter("pBC",&m_pBC);
  
  m_VrBC = 1935.07; //848.15;
  setParameter("VrBC",&m_VrBC);
  
  m_rotation = false; 
  setParameter("Rotation",&m_rotation);
  
  m_enforceFullB = false;
  setParameter("EnforceFullB",&m_enforceFullB);
  
  m_enforeDipoleB = false;
  setParameter("EnforceDipoleB",&m_enforeDipoleB);
  
  m_initialSolutionIDs = std::vector<CFuint>();
  setParameter("InitialSolutionIDs",&m_initialSolutionIDs);
}

//////////////////////////////////////////////////////////////////////////////

BCSuperInletProjMHD::~BCSuperInletProjMHD()
{
  CFAUTOTRACE;
  
  if (m_initialSolutionMap.size() > 0) 
  {
    for (CFuint i = 0; i < m_initialSolutionMap.size(); ++i) 
    {
      deletePtr(m_initialSolutionMap[i]);
    }
  }
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
  }
  m_cellStatesFlxPnt.clear();
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("RhoBC","Boundary rho value.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("pBC","Boundary p value.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("VrBC","radial velocity at the boundary");
  options.addConfigOption< bool >("Rotation","rotation");
  options.addConfigOption< bool >("EnforceFullB","Enforce full B_in field.");
  options.addConfigOption< bool >("EnforceDipoleB","Enforce analytical dipole B_in field.");
  options.addConfigOption< std::vector<CFuint> > ("InitialSolutionIDs", "IDs of initial solution components that will be used as BC value.");
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::preProcess()
{  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  if (iter == 1)
  {
  
  const CFuint initialSSize = m_initialSolutionIDs.size();
  
  if (initialSSize > 0)// && this->getMethodData().isPreProcessedSolution()) 
  {    
    if (!m_initialSolutionMap.exists(m_thisTRS->getName()) && m_thisTRS->getLocalNbGeoEnts() > 0) 
    {
      CFLog(INFO, "BCSuperInletProjMHD::preProcess() => TRS[ " << m_thisTRS->getName() << " ] => START\n");	
      
      // get InnerCells TopologicalRegionSet
      SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

      // get bndFacesStartIdxs from SpectralFRMethodData
      map< std::string , vector< vector< CFuint > > >& bndFacesStartIdxsPerTRS = getMethodData().getBndFacesStartIdxs();
      vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[m_thisTRS->getName()];

      // number of face orientations (should be the same for all TRs)
      cf_assert(bndFacesStartIdxs.size() != 0);
      const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

      // number of TRs
      const CFuint nbTRs = m_thisTRS->getNbTRs();
      cf_assert(bndFacesStartIdxs.size() == nbTRs);

      // get the geodata of the face builder and set the TRSs
      FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
      geoData.cellsTRS = cellTrs;
      geoData.facesTRS = m_thisTRS;
      geoData.isBoundary = true;
      
      const CFuint nbTrsFaces = m_thisTRS->getLocalNbGeoEnts();
          
      RealVector* initialState = new RealVector(nbTrsFaces*m_nbrFaceFlxPnts*initialSSize);
      cf_assert(initialState != CFNULL);

      // get the geodata of the cell builders and set the TRS
      //CellToFaceGEBuilder::GeoData& geoDataCB = m_cellBuilder->getDataGE();
      //geoDataCB.m_thisTRS = cellTrs;
      
      //CFuint counter = 0;
      CFuint localFaceID = 0;
      
      const CFuint nbVars = m_initialSolutionIDs.size();

      // loop over TRs
      for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
      {
        // loop over different orientations
        for (m_orient = 0; m_orient < nbOrients; ++m_orient)
        {
          // start and stop index of the faces with this orientation
          const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
          const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

          // loop over faces with this orientation
          for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
          {
	    // build the face GeometricEntity
            geoData.idx = faceID;
            m_currFace = m_faceBuilder->buildGE();
            
            const CFuint faceGlobalID = m_currFace->getID();
            m_globalToLocalTRSFaceID.insert(faceGlobalID,localFaceID);
            
            //const CFuint faceID = m_currFace->getID();

            // GET THE NEIGHBOURING CELL
            m_intCell = m_currFace->getNeighborGeo(0);

            // GET THE STATES IN THE NEIGHBOURING CELL
            m_cellStates = m_intCell->getStates();
	
            // BUILD THE CELL WITH CONNECTIVITY TO ITS FACES
            //geoDataCB.idx = m_intCell->getID();
            //m_intCell = m_cellBuilder->buildGE();

	    const CFuint nbrStates = m_cellStates->size();
  
            // Loop over flux points to extrapolate the states and gradients to the flux points
            for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
            {
              // reset the states in flx pnts
              *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;
              
              const CFuint startID = localFaceID*m_nbrFaceFlxPnts*nbVars + iFlxPnt*nbVars;

              // index of current flx pnt
              const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
              // Loop over sol points to add the contributions to each sol pnt
              for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
              {
                *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_cellStates)[iSol]));
              }
              
              const State& innerState = *(m_cellStatesFlxPnt[iFlxPnt]);
                            
              for (CFuint i = 0; i < nbVars; ++i)//, ++counter) 
              {
	        const CFuint varID = m_initialSolutionIDs[i];
	        //cf_assert(counter < initialState->size());
	  
	        cf_assert(varID < innerState.size());
	        (*initialState)[startID + i] = innerState[varID];
	      }
            }	

	    // release the face and the cell
            m_faceBuilder->releaseGE();
            //m_cellBuilder->releaseGE();
	
            localFaceID++;
          } 
        }
      }
      
      m_initialSolutionMap.insert(m_thisTRS->getName(), initialState);
            
      CFLog(INFO, "BCSuperInletProjMHD::preProcess() => TRS[ " << m_thisTRS->getName() << " ] => END\n");
    }
  }
  else
  {
    CFLog(INFO,"BCSuperInletProjMHD::preProcess() => No variable IDs specified for potential field tranfer!\n");
    cf_assert(false);
  } 
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  //bool is3D = (normals[0]).size()==3;

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 9 || 7);
    cf_assert(ghostState.size() == 9 || 7);

    // during initialization phase we store the initial solution values to be used as BC
    // This BC (the PFSS solution fixed on the inlet) was suggested by Jon Linker and
    // is evaluated by default. When using Jens' or Dana's BCs for the magnetic field
    // the B-field values in the ghost cells are overwritten.

    // initialize PFSS B solution vector
    RealVector B_PFSS_dimless(3);

    // get PFSS B solution from m_initialSolutionMap
    if (m_initialSolutionMap.size() > 0) 
    {
      /// map faces to corresponding TRS and index inside that TRS
      //SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs = MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");

      //const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(m_face->getID());
      
      // Get the localFaceID from the map, knowing the faceGlobalID
      const CFuint faceLocalID = m_globalToLocalTRSFaceID.find(m_face->getID());

      const string name = m_thisTRS->getName();

      SafePtr<RealVector> initialValues = m_initialSolutionMap.find(name);

      const CFuint nbVars = m_initialSolutionIDs.size();//is3D ? 3 : 2;//
      const CFuint startID = faceLocalID*nbVars*m_nbrFaceFlxPnts + iState*nbVars;

      for (CFuint i = 0; i < nbVars; ++i) 
      {
        const CFuint varID = m_initialSolutionIDs[i];//is3D ? i+4 : i+3;//
        const CFuint idx = startID+i;

        cf_assert(idx < initialValues->size());

        // save the Poisson PFSS solution in a vector
        B_PFSS_dimless[i] = (*initialValues)[idx]; 

        //(*ghostState)[varID] = 2.*(*initialValues)[idx] - (*innerState)[varID];
        //CFLog(DEBUG_MIN, "SuperInletProjectionParallel5::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
      }
    }
    
    const RealVector& normal = normals[iState];
      
    const CFreal xI_dimless = coords[iState][XX];
    const CFreal yI_dimless = coords[iState][YY];
    const CFreal zI_dimless = coords[iState][ZZ];
    const CFreal rI_dimless = (coords[iState]).norm2();

    // to be checked if imposed minimum is needed (suspected problems at poles otherwise)
    const CFreal rhoI_dimless = std::max(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),1.0e-4);
    
    // density
    ghostState[0] = 2.0*m_rhoBC - intState[0];
    
    // p
    ghostState[7] = 2.0*m_pBC - intState[7];

    // phi
    ghostState[8] = intState[8];
    
    // B
    // Br bnd from PFSS solution
    CFreal BrBoundary_dimless = xI_dimless/rI_dimless*B_PFSS_dimless[0] + yI_dimless/rI_dimless*B_PFSS_dimless[1] + zI_dimless/rI_dimless*B_PFSS_dimless[2];

    if (m_enforeDipoleB) BrBoundary_dimless = 2.0/3.0*0.666*zI_dimless;
    
    const CFreal BxI_dimless = intState[4];
    const CFreal ByI_dimless = intState[5];
    const CFreal BzI_dimless = intState[6];
    const CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
    
    CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
    
    if (m_enforceFullB) BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*B_PFSS_dimless[0] + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*B_PFSS_dimless[1] - rhoI_dimless/rI_dimless*B_PFSS_dimless[2];
    if (m_enforeDipoleB) BthetaI_dimless = 0.666/3.0*rhoI_dimless;
    
    CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;
    
    if (m_enforceFullB) BphiI_dimless = -yI_dimless/rhoI_dimless*B_PFSS_dimless[0] + xI_dimless/rhoI_dimless*B_PFSS_dimless[1];
    if (m_enforeDipoleB) BphiI_dimless = 0.0;

    // B bnd, Btheta and Bphi bnd are set equal to inner values
    const CFreal BxBoundary_dimless = xI_dimless/rI_dimless*BrBoundary_dimless - yI_dimless/rhoI_dimless*BphiI_dimless + xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BthetaI_dimless;
    const CFreal ByBoundary_dimless = yI_dimless/rI_dimless*BrBoundary_dimless + xI_dimless/rhoI_dimless*BphiI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BthetaI_dimless;
    const CFreal BzBoundary_dimless = zI_dimless/rI_dimless*BrBoundary_dimless - rhoI_dimless/rI_dimless*BthetaI_dimless;

    ghostState[4] = 2.0*BxBoundary_dimless - BxI_dimless;
    ghostState[5] = 2.0*ByBoundary_dimless - ByI_dimless;
    ghostState[6] = 2.0*BzBoundary_dimless - BzI_dimless;
    
    // V
    //ghostState[1] = 0.0;//-intState[1];
    //ghostState[2] = 0.0;//-intState[2];
    //ghostState[3] = 0.0;//-intState[3]; 
    
    const CFreal BxB = BxBoundary_dimless;
    const CFreal ByB = ByBoundary_dimless;
    const CFreal BzB = BzBoundary_dimless; 
    const CFreal BBmag = std::sqrt(BxB*BxB + ByB*ByB + BzB*BzB);
    const CFreal VBmag = m_VrBC/(2.2e-4/sqrt(1.2566e-6*1.67e-13)); //1935.07/(2.2e-4/sqrt(1.2566e-6*1.67e-13));
    const CFreal BxuB = BxB/BBmag;
    const CFreal ByuB = ByB/BBmag;
    const CFreal BzuB = BzB/BBmag;
    CFreal VxBB = BxuB * VBmag;
    CFreal VyBB = ByuB * VBmag;
    CFreal VzBB = BzuB * VBmag;
    
    CFreal VrBB = xI_dimless/rI_dimless*VxBB + yI_dimless/rI_dimless*VyBB + zI_dimless/rI_dimless*VzBB;

    //If the radial component would indicate inflow, switch the orientation
    if (VrBB < 0.0)
    {
      VxBB = -VxBB;
      VyBB = -VyBB;
      VzBB = -VzBB;
    }

    //With the fixed-orientation components, recalculate the spherical components
    VrBB = xI_dimless/rI_dimless*VxBB + yI_dimless/rI_dimless*VyBB + zI_dimless/rI_dimless*VzBB;
    CFreal VthetaBB = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxBB + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyBB - rhoI_dimless/rI_dimless*VzBB;
    CFreal VphiBB = -yI_dimless/rhoI_dimless*VxBB + xI_dimless/rhoI_dimless*VyBB;

    //Add the rotation motion
    if (m_rotation)
    {
      VphiBB += rI_dimless*std::sin(BthetaI_dimless)*3.86e-3;
    } 

    //Recompute Cartesian components
    const CFreal VxBBR = xI_dimless/rI_dimless*VrBB - yI_dimless/rhoI_dimless*VphiBB + xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VthetaBB;
    const CFreal VyBBR = yI_dimless/rI_dimless*VrBB + xI_dimless/rhoI_dimless*VphiBB + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VthetaBB;
    const CFreal VzBBR = zI_dimless/rI_dimless*VrBB - rhoI_dimless/rI_dimless*VthetaBB;  

    ghostState[1] = 2.0*VxBBR - intState[1];
    ghostState[2] = 2.0*VyBBR - intState[2];
    ghostState[3] = 2.0*VzBBR - intState[3];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  const CFuint nbrGradVars = intGrads[0].size();

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletProjMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = true;

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  m_nbrFaceFlxPnts = frLocalData[0]->getFaceFlxPntsFaceLocalCoords()->size();

  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();

  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt.push_back(new State());
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
  }

  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    if (m_trsNames[0]==trsList[iTRS]->getName())
    {
      m_thisTRS = trsList[iTRS];
      CFLog(INFO, "Inlet: Matching BC "<<m_trsNames[0]<<" with "<<m_thisTRS->getName() << "\n");
    }
  }

//  // get MHD 3D varset
//  m_varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
//  if (m_varSet.isNull())
//  {
//    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in BCSuperInletProjMHD!");
//  }
//
//  // resize the physical data for internal and ghost solution points
//  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
//  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

