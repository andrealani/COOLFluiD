#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"

#include "AeroCoef/AeroCoef.hh"
#include "AeroCoef/AeroForcesFR.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::FluxReconstructionMethod;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////
      
MethodCommandProvider<AeroForcesFR, DataProcessingData, AeroCoefModule>
aeroForcesFRProvider("AeroForcesFR");

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Alpha","Function defining the angle of attack.");
  options.addConfigOption< std::string >("Beta","Function defining the sideslip angle.");
  options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
  options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
  options.addConfigOption< CFreal >("uInf","Velocity at infinity.");
  options.addConfigOption< CFreal >("rhoInf","Density at infinity.");
  options.addConfigOption< CFreal >("pInf","Pressure at infinity.");
  options.addConfigOption< CFreal >("TInf","Value of the fresstream temperature (in the dimensional case)");
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name."); 
  options.addConfigOption< bool >("ReorderWallData","Reorder the wall data to make the file structured.");
  options.addConfigOption< CFuint >("TID","Position of T in the state vector");
  options.addConfigOption< CFuint >("UID","Position of u in the state vector");
  options.addConfigOption< CFuint >("VID","Position of v in the state vector");
  options.addConfigOption< CFuint >("WID","Position of w in the state vector");
  options.addConfigOption< CFreal >("RefLength2D","Reference length (e.g. chord) for scaling 2D aerodynamic coefficients");
  options.addConfigOption< CFreal >("RefArea","Reference area (e.g. chord) for scaling aerodynamic coefficients");
  options.addConfigOption< vector<CFreal> >("GravityCenter","Gravity center coordinates for computing aerodynamic moment");
  options.addConfigOption< std::string >("OutputFileConv","Name of convergence file for surface quantities.");
}

//////////////////////////////////////////////////////////////////////////////

AeroForcesFR::AeroForcesFR(const std::string& name) :
  DataProcessingCom(name),
  socket_gradients("gradients"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  m_sockets(),
  socket_states("states"), 
  m_updateVarSet(CFNULL),
  m_frData(CFNULL),
  m_mapTrsFaceToID(),
  m_dataState(),
  m_coord(),
  m_frictionForces(),
  m_xyzMoment(),
  m_v01(),
  m_v02(),
  m_vCross0102(),
  m_v12(),
  m_xg(),
  m_midFaceNode(),
  m_rotMat(),
  m_xyzForce(),
  m_force(),
  m_aeroForce(),
  m_valuesMatL2(),
  m_l2Norm(),
  m_valuesMat(),
  m_valuesMatRes(),
  m_varNames(),
  m_functionAlphaParser(),
  m_functionBetaParser(),
  m_vars("t"),
  m_eval(0.0,1),
  m_currFace(CFNULL),
  m_iFace(0),
  m_p(),
  m_Cp(0.),
  m_Cf(0.),
  m_Mach(),
  m_lift (),
  m_lateral(),
  m_drag(),
  m_wetSurface(),
  m_alphadeg(0.),
  m_alpha(0.),
  m_betadeg(0.),
  m_beta(0.),
  m_faceBuilder(CFNULL),
  m_cellBuilder(CFNULL),
  m_isFaceOnBoundary(CFNULL),
  m_nghbrCellSide(CFNULL),
  m_currCellSide(CFNULL),
  m_faceOrients(CFNULL),
  m_faceBCIdx(CFNULL),
  m_orient(),
  m_intCell(CFNULL),
  m_cellStates(CFNULL),
  m_nbrFaceFlxPnts(),
  m_flxLocalCoords(CFNULL),
  m_cellGrads(),
  m_cellGradFlxPnt(),
  m_faceFlxPntConn(CFNULL),
  m_faceConnPerOrient(CFNULL),
  m_faceIntegrationCoefs(CFNULL),
  m_cellStatesFlxPnt(),
  m_unitNormalFlxPnts(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_dim(),
  m_faceJacobVecAbsSizeFlxPnts(),
  m_faceMappedCoordDir(CFNULL),
  m_faceJacobVecSizeFlxPnts()
{
  addConfigOptionsTo(this);

  m_nameOutputFileWall = "wall.plt";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_nameOutputFileAero = "aeroCoef.plt";
  setParameter("OutputFileAero",&m_nameOutputFileAero);
  
  m_functionAlpha = std::string("0.");
  setParameter("Alpha",&m_functionAlpha);

  m_functionBeta = std::string("0.");
  setParameter("Beta",&m_functionBeta);

  m_uInf = 0.;
  setParameter("uInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("rhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("pInf",&m_pInf);

  m_TInf = 0.;
  setParameter("TInf",&m_TInf);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_reorderWallData = true;
  setParameter("ReorderWallData",&m_reorderWallData);
  
  m_TID = 0;
  setParameter("TID",&m_TID);

  m_UID = 0;
  setParameter("UID",&m_UID);

  m_VID = 0;
  setParameter("VID",&m_VID);

  m_WID = 0;
  setParameter("WID",&m_WID);

  m_refLength2D = 1.;
  setParameter("RefLength2D",&m_refLength2D);

  m_refArea = 0.;
  setParameter("RefArea",&m_refArea);

  m_gravityCenter = vector<CFreal>();
  setParameter("GravityCenter",&m_gravityCenter);
  
  m_outputFileConv = "convergence-surf.plt";
  setParameter("OutputFileConv",&m_outputFileConv);
}

//////////////////////////////////////////////////////////////////////////////

AeroForcesFR::~AeroForcesFR()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
AeroForcesFR::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  result.push_back(&socket_states);
  result.push_back(&socket_gradients);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
AeroForcesFR::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::setup()
{
  DataProcessingCom::setup();
  
  if (Common::PE::GetPE().IsParallel()) {
    // there is only one file for the aerodynamic coefficients
    m_nameOutputFileAero = m_nameOutputFileAero + "-P0";
    
    // there is only one file for the wall data
    m_nameOutputFileWall = m_nameOutputFileWall + "-P0";
  }

  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluxReconstructionSolver> fr = spaceMethod.d_castTo<FluxReconstructionSolver>();
  cf_assert(fr.isNotNull());

  m_frData = fr->getData();
  m_updateVarSet = m_frData->getUpdateVar().d_castTo<EulerVarSet>();
  cf_assert(m_updateVarSet.isNotNull());
  m_updateVarSet->getModel()->resizePhysicalData(m_dataState);
  
  // dimensionality and number of equations
  m_dim = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  m_coord.resize(m_dim);
  m_frictionForces.resize(m_dim);
  m_xyzMoment.resize((m_dim == DIM_2D) ? 1 : m_dim);

  m_v01.resize(m_dim);
  m_v02.resize(m_dim);
  m_vCross0102.resize(m_dim);
  m_v12.resize(m_dim);
  m_xg.resize(m_dim);
  m_midFaceNode.resize(m_dim);
  m_rotMat.resize(m_dim,m_dim);
  m_xyzForce.resize(m_dim);
  m_force.resize(m_dim);
  m_aeroForce.resize(m_dim);
  
  // get face builder
  m_faceBuilder = m_frData->getFaceBuilder();

  // get cell builder
  m_cellBuilder = m_frData->getCellBuilder();

  // get some additional data for cell building
  m_isFaceOnBoundary = m_cellBuilder->getGeoBuilder()->getIsFaceOnBoundary();
  m_nghbrCellSide    = m_cellBuilder->getGeoBuilder()->getNeighbrCellSide ();
  m_currCellSide     = m_cellBuilder->getGeoBuilder()->getCurrentCellSide ();
  m_faceOrients      = m_cellBuilder->getGeoBuilder()->getFaceOrient      ();
  m_faceBCIdx        = m_cellBuilder->getGeoBuilder()->getFaceBCIdx       ();

  cf_always_assert(m_uInf > 0.);
  cf_always_assert(m_rhoInf > 0.);
  cf_always_assert(m_pInf > 0.);
  
  const std::string nsp = this->getMethodData().getNamespace();
  if (m_gravityCenter.size() == 0) {
    m_gravityCenter.resize(m_dim, 0.0);
    if (PE::GetPE().GetRank(nsp) == 0) {
      cout << "#### WATCH OUT: gravity center not specified !!! assuming X=0 #####" << endl;
    }
  }

  for (CFuint i =0; i < m_dim; ++i) {
    m_xg[i]= m_gravityCenter[i];
  }

  // pre-computation of the total wet surface
  if (m_dim == DIM_3D) {
     computeWetSurface();
  }
  else {
    cf_always_assert(m_refLength2D > 0.);
    m_wetSurface = m_refLength2D;
  }

  // reference area is set to be the wet area if it is not specified by the user
  m_refArea = (m_refArea > 0.) ? m_refArea : m_wetSurface;
  
  m_outputFileAeroPrepared = false;
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // compute flux point coordinates
  m_flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  m_nbrFaceFlxPnts = m_flxLocalCoords->size();
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
   
  // get the face - flx pnt connectivity per orient
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
	  
  // get the face connectivity per orientation
  m_faceConnPerOrient = frLocalData[0]->getFaceConnPerOrient();
  
  // get the face integration coefficient
  m_faceIntegrationCoefs = frLocalData[0]->getFaceIntegrationCoefs();
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  // get flux point mapped coordinate directions
  m_faceMappedCoordDir = frLocalData[0]->getFaceMappedCoordDir();
  
  // resize m_unitNormalFlxPnts
  m_unitNormalFlxPnts.resize(m_nbrFaceFlxPnts);
  m_cellGradFlxPnt.resize(m_nbrFaceFlxPnts);
  
  // resize m_faceJacobVecSizeFlxPnts
  m_faceJacobVecSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  m_faceJacobVecAbsSizeFlxPnts.resize(m_nbrFaceFlxPnts);
  
  // create internal and ghost states
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt.push_back(new State());
    m_unitNormalFlxPnts[iFlx].resize(m_dim);
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      m_cellGradFlxPnt[iFlx].push_back(new RealVector(m_dim));
    }
  }

  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    m_cellStatesFlxPnt[iFlx]->setLocalID(iFlx);
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::computeWetSurface()
{
//   vector<SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
// 
//   Common::SafePtr<GeometricEntityPool<FaceToCellGEBuilder> >
//     geoBuilder = m_faceBuilder;
// 
//   //SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
//   //geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
//   
//   DataHandle< CFreal> faceAreas = socket_faceAreas.getDataHandle();
//   
//   FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   geoData.isBFace = true;
// 
//   CFreal wSurface = 0.0;
//   for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
//     geoData.trs = trsList[iTRS];
//     
//     const CFuint nbTrsFaces = trsList[iTRS]->getLocalNbGeoEnts();
//     for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
//       // build the GeometricEntity
//       geoData.idx = iFace;
//       GeometricEntity* const currFace = geoBuilder->buildGE();
//       
//       if (currFace->getState(0)->isParUpdatable()) {
//  	wSurface += faceAreas[currFace->getID()]; 
//       }
//       
//       // release the geometric entity
//       geoBuilder->releaseGE();
//     }
//   }
//   
//   // sum all the contributions on this TRS from all processors
//   const std::string nsp = this->getMethodData().getNamespace();
//   MPI_Allreduce(&wSurface, &m_wetSurface, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));     
//   
//   cf_always_assert(m_wetSurface > 0.0);
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // Create the socket source for the boundaryNormals for each TRS to which this command applies
  // This is because we have different normals in FiniteVolume and FluctSplit / FiniteElement
  const CFuint nbTrs = getTrsNames().size();

  for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs) {
    std::string socketName = getTrsName(iTrs) + "-boundaryNormals";
    m_sockets.createSocketSource<const CFreal*>(socketName);
  }

  try {
    setFunction();
  }
  catch (Common::Exception& e) {
    CFout << e.what() << "\n";
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::setFunction()
{
  // some sanity checks on alpha
  std::vector<std::string> m_functionAlphaDef = Common::StringOps::getWords(m_functionAlpha);
  if(m_functionAlphaDef.size() != 1)
    throw BadValueException (FromHere(),"AeroForcesFR::setFuntion(): alpha function wrongly defined.");

  m_functionAlphaParser.Parse(m_functionAlpha, m_vars);

  if (m_functionAlphaParser.ErrorMsg() != 0) {
    std::string msg("ParseError in AeroForcesFR::setFuntion(): ");
    msg += std::string(m_functionAlphaParser.ErrorMsg());
    msg += " Function: " + m_functionAlpha;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }

  // some sanity checks on beta
  std::vector<std::string> m_functionBetaDef = Common::StringOps::getWords(m_functionBeta);
  if(m_functionBetaDef.size() != 1)
    throw BadValueException (FromHere(),"AeroForcesFR::setFuntion(): alpha function wrongly defined.");

  m_functionBetaParser.Parse(m_functionBeta, m_vars);

  if (m_functionBetaParser.ErrorMsg() != 0) {
    std::string msg("ParseError in AeroForcesFR::setFuntion(): ");
    msg += std::string(m_functionBetaParser.ErrorMsg());
    msg += " Function: " + m_functionBeta;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::executeOnTrs()
{
  CFAUTOTRACE;
  
  const std::string nsp = this->getMethodData().getNamespace();
  if (PE::GetPE().GetRank(nsp) == 0) {
    CFLog(VERBOSE, "AeroForcesFR::executeOnTrs() START\n");
  }
  
  //  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute the value of the angle
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionAlphaParser.Eval(&m_eval[0]);
  m_betadeg  = m_functionBetaParser.Eval(&m_eval[0]);
  
  // transform into Radiants
  m_alpha = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;
  m_beta  = m_betadeg*MathTools::MathConsts::CFrealPi()/180;

  // check
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // rotation matrix for 2D with alpha>0 in (x,y) plane
  if (dim == DIM_2D) {
    m_rotMat(XX,XX) = cos(m_alpha);   m_rotMat(XX,YY) = sin(m_alpha);
    m_rotMat(YY,XX) = -sin(m_alpha);  m_rotMat(YY,YY) = cos(m_alpha);
  }

  // rotation matrix for 3D with alpha>0 in (x,z) plane, beta>0 in (x,y) plane
  if (dim == DIM_3D) {
    m_rotMat(XX,XX) = cos(m_alpha)*cos(m_beta);  m_rotMat(XX,YY) = sin(m_beta);   m_rotMat(XX,ZZ) = sin(m_alpha)*cos(m_beta);
    m_rotMat(YY,XX) = -sin(m_beta)*cos(m_alpha); m_rotMat(YY,YY) = cos(m_beta);   m_rotMat(YY,ZZ) = -sin(m_beta)*sin(m_alpha);
    m_rotMat(ZZ,XX) = -sin(m_alpha);             m_rotMat(ZZ,YY) = 0.;            m_rotMat(ZZ,ZZ) = cos(m_alpha);
  }

//   SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
// 
//   // check
//   Common::SafePtr<GeometricEntityPool<FaceCellTrsGeoBuilder> > 
//     geoBuilder = m_frData->getFaceCellTrsGeoBuilder();
// 
//   SafePtr<FaceCellTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
//   geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
//   
//   FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   
//   geoData.faces = currTrs;
//   geoData.isBFace = true;
//   geoData.cells = MeshDataStack::getActive()->getTrs("InnerCells");
//   geoData.allCells = true;
//   // check/
//   
//   const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();

  if (m_valuesMatRes.size() == 0) {initSurfaceResiduals();}
  prepareOutputFileAero();
  prepareOutputFileWall();

//   // this needs to be updated with the latest values
//   m_frData->getPolyReconstructor()->computeGradients();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

//   // get FRMethodData
//   SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
//   SafePtr<FluxReconstructionSolver> fr = spaceMethod.d_castTo<FluxReconstructionSolver>();
//   cf_assert(fr.isNotNull());
//   SafePtr<FluxReconstructionSolverData> frData = fr->getData();

  // get bndFacesStartIdxs from SpectralFDMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = m_frData->getBndFacesStartIdxs();
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

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

//       // SET THE ORIENTATION OF THE FACES
//       m_bndFaceTermComputer->setFaceOrientation(m_orient);
//       m_bndFaceTermComputer->setPointSet(*outputPntsMappedCoord);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
	// build the face GeometricEntity
        geoData.idx = faceID;
        m_currFace = m_faceBuilder->buildGE();

        // GET THE NEIGHBOURING CELL
        m_intCell = m_currFace->getNeighborGeo(0);

        // GET THE STATES IN THE NEIGHBOURING CELL
        m_cellStates = m_intCell->getStates();
	
	// if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // BUILD THE CELL WITH CONNECTIVITY TO ITS FACES
          geoDataCB.idx = m_intCell->getID();
          m_intCell = m_cellBuilder->buildGE();

	  computeWallStatesGrads();

	  computeWall();

          computeAero();
	}

	// release the face and the cell
        m_faceBuilder->releaseGE();
        m_cellBuilder->releaseGE();
	
      } 
    }
  }
  
//   for (m_iFace = 0; m_iFace < nbTrsFaces; ++m_iFace) {
//     CFLog(VERBOSE, "iFace = " << m_iFace << "\n");
// 
//     // build the GeometricEntity
//     geoData.idx = m_iFace;
//     m_currFace = geoBuilder->buildGE();
//     m_fvmccData->getCurrentFace() = m_currFace;
// 
//     // compute the face normal
//     const CFuint faceID = m_currFace->getID();
//     const CFuint startID = faceID*dim;
//     // face normal is pointing outward from the computational domain
//     // which means that is pointing into the body
//     // we consider the opposite normal, pointing inward the domain
//     for (CFuint i = 0; i < dim; ++i) {
//       m_normal[i] = -normals[startID + i];// this IS NOT a unit normal : it is \vec{1}_n \delta S
//     }
//     m_unitNormal = m_normal*(1./m_normal.norm2());// this IS the unit normal
//     
//     // this is fundamental because some algorithms for computing numerical derivatives need this
//     m_fvmccData->getUnitNormal() = -m_unitNormal;  
//     
//     
//     
//     // release the geometric entity
//     geoBuilder->releaseGE();
//   }

  // all data are written on file at once to ease the parallel writing
  updateOutputFileAero(); 

  updateOutputFileWall();

  reorderOutputFileWall();

  computeSurfaceResiduals(); 
  
  if (PE::GetPE().GetRank(nsp) == 0) {
    CFLog(VERBOSE, "AeroForcesFR::executeOnTrs() END\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::computeWallStatesGrads()
{
  CFAUTOTRACE;
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();
  m_cellGrads.resize(nbrStates);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    m_cellGrads[iState] = &gradients[stateID];
  }

  // compute face Jacobian vectors
  vector< RealVector > faceJacobVecs = m_currFace->computeFaceJacobDetVectorAtMappedCoords(*m_flxLocalCoords);
  
  // get face Jacobian vector sizes in the flux points
  DataHandle< vector< CFreal > >
  faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  
  // Loop over flux points to extrapolate the states and gradients to the flux points
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // reset the states in flx pnts
    *(m_cellStatesFlxPnt[iFlxPnt]) = 0.0;

    // reset the grads in flx pnts
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *(m_cellGradFlxPnt[iFlxPnt][iVar]) = 0.0;
    }

    // index of current flx pnt
    const CFuint currFlxIdx = (*m_faceFlxPntConn)[m_orient][iFlxPnt];
    
    // Loop over sol points to add the contributions to each sol pnt
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      *(m_cellStatesFlxPnt[iFlxPnt]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*(*((*m_cellStates)[iSol]));
      
      for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
      {
        *(m_cellGradFlxPnt[iFlxPnt][iVar]) += (*m_solPolyValsAtFlxPnts)[currFlxIdx][iSol]*((*(m_cellGrads[iSol]))[iVar]);
      }
    }

    const CFuint faceID = m_currFace->getID();
    
    // get face Jacobian vector size
    m_faceJacobVecAbsSizeFlxPnts[iFlxPnt] = faceJacobVecSizeFaceFlxPnts[faceID][iFlxPnt];

    // set face Jacobian vector size with sign depending on mapped coordinate direction
    m_faceJacobVecSizeFlxPnts[iFlxPnt] = m_faceJacobVecAbsSizeFlxPnts[iFlxPnt]*(*m_faceMappedCoordDir)[m_orient];

    // set unit normal vector
    m_unitNormalFlxPnts[iFlxPnt] = faceJacobVecs[iFlxPnt]/m_faceJacobVecAbsSizeFlxPnts[iFlxPnt];
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::computeWall()
{
  CFAUTOTRACE;

//   DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
//   const vector<Node*>& nodesInFace = *m_currFace->getNodes();
//   const CFuint nbNodesInFace = nodesInFace.size();
// 
//   *m_avState = 0.0;
//   m_coord = 0.0;
//   for (CFuint i = 0; i < nbNodesInFace; ++i) {
//     *m_avState += nstates[nodesInFace[i]->getLocalID()];
//     m_coord += *nodesInFace[i];
//   }
//   *m_avState /= nbNodesInFace;
//   m_coord /= nbNodesInFace;
// 
//   m_updateVarSet->computePhysicalData(*m_avState, m_dataState);
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::updateOutputFileWall()
{
/*  using namespace boost::filesystem;
  path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
 file = Framework::PathAppender::getInstance().appendAllInfo 
      (file,m_appendIter,m_appendTime,false);  
       
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(file, ios::app);

  const CFreal time = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const CFreal mach = m_dataState[EulerTerm::V]/m_dataState[EulerTerm::A];

  // Output to File
  fout << m_coord    << " "
       << iter << " "
       << time << " "
       << m_p  << " "
       << mach << " "
       << m_Cp << " "
       << m_dataState[EulerTerm::T] << " "
       << m_dataState[EulerTerm::RHO] << "\n";
*/
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::computeAero()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

//   // this is redundant but allows to adapt all the subclasses
//   const vector<Node*>& nodesInFace = *m_currFace->getNodes();
//   const CFuint nbNodesInFace = nodesInFace.size();

//   // compute all physical values corresponding to the given states
//   m_midFaceNode = 0.0;
//   for (CFuint i = 0; i < nbNodesInFace; ++i) {
//     m_midFaceNode += *nodesInFace[i];
//   }
//   m_midFaceNode /= nbNodesInFace;
  
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    // compute coordinates of output point
    m_coord = m_currFace->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
    
    cf_assert(m_updateVarSet.isNotNull());
    m_updateVarSet->computePhysicalData(*(m_cellStatesFlxPnt[iFlx]), m_dataState);

    // compute pressure and Cp (redundant for safety)
    m_p = m_dataState[EulerTerm::P] * m_updateVarSet->getModel()->getPressRef();
    m_Cp = (m_p - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

    // compute shortest distance from GIVEN centre of gravity and segment line passing
    // for inner and ghost states corresponding to the current boundary face
    // from http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    //  m_v01 = m_xg - m_currFace->getState(1)->getCoordinates();
    //   m_v02 = m_xg - m_currFace->getState(0)->getCoordinates();
    //   MathTools::MathFunctions::crossProd(m_v01, m_v02, m_vCross0102);
    //   m_v12 = m_currFace->getState(1)->getCoordinates() - m_currFace->getState(0)->getCoordinates();
    //   const CFreal distanceToForce = m_vCross0102.norm2()/m_v12.norm2();

    // pressure is pointing inward with respect to the body, so it is opposite
    // with respect to the normal (pointing inside the domain)

    // \Delta\vec{F} = (-p \delta_{ij} \cdot \vec{n} + \sigma_{ij} \cdot \vec{n}) \Delta S
    m_xyzForce  = (-m_Cp/m_refArea)*m_unitNormalFlxPnts[iFlx] + m_frictionForces;
    m_aeroForce = m_rotMat*m_xyzForce; // [D C L]^T = {rotation matrix} * [X Y Z]^T

    // vector between center of gravity and centre of the face where the force is applied
    m_v12 = m_coord - m_xg;

    // assemble the contribution of this face only if the internal state is parallel updatable 
    if ((*m_cellStates)[0]->isParUpdatable()) 
    {
      if (m_dim == DIM_2D) 
      {
        // this is the moment with respect to Z axis (Mz > 0 => pitching down)
        m_xyzMoment[0] += (m_v12[XX]*m_xyzForce[YY] - m_v12[YY]*m_xyzForce[XX])/m_refLength2D;
      }
      else 
      {
        MathFunctions::crossProd(m_v12, m_xyzForce, m_vCross0102);
        m_xyzMoment += m_vCross0102/m_refLength2D;
      }

      m_force    += m_xyzForce;
      m_drag     += m_aeroForce[XX];
      m_lateral  += (m_dim == DIM_2D) ? 0. : m_aeroForce[YY];
      m_lift     += m_aeroForce[(m_dim == DIM_2D) ? YY : ZZ];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::unsetup()
{
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    deletePtr(m_cellStatesFlxPnt[iFlx]);
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      deletePtr(m_cellGradFlxPnt[iFlx][iGrad]); 
    }
    m_cellGradFlxPnt[iFlx].clear();
  }
  m_cellStatesFlxPnt.clear();
  
  DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::prepareOutputFileWall()
{
  CFAUTOTRACE;
  
  if (getCurrentTRS()->getLocalNbGeoEnts() > 0) {
    using namespace boost::filesystem;
    path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
    file = Framework::PathAppender::getInstance().appendAllInfo 
      (file,m_appendIter,m_appendTime,false);  
       
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(file);
    
    fout << "TITLE = Unstructured Surface Quantities" << "\n";
    if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
      fout << "VARIABLES = x y z Iter Time Pressure Mach Cp Temperature Density" << "\n";
    }
    else {
      fout << "VARIABLES = x y Iter Time Pressure Mach Cp Temperature Density" << "\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::prepareOutputFileAero()
{
  CFAUTOTRACE;
  
  const std::string nsp = this->getMethodData().getNamespace();
  if (PE::GetPE().GetRank(nsp) == 0 && !m_outputFileAeroPrepared) {
    // preparation of the output
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(m_nameOutputFileAero);
    file = Framework::PathAppender::getInstance().appendAllInfo 
      (file,false,false,false);  
       
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(file);
    
    //Writing the title
    fout << "TITLE = Aerodynamics Coefficients" << "\n";
    
    if (PhysicalModelStack::getActive()->getDim() == DIM_2D) {
      fout << "VARIABLES = Iter PhysTime Alpha Beta CD CC Cx Cy CL CM" << "\n";
    }
    else {
      fout << "VARIABLES = Iter PhysTime Alpha Beta CD CC CL Cx Cy Cz CMx CMy CMz" << "\n";
    }
    fhandle->close(); 
    
    m_outputFileAeroPrepared = true;
  }

  m_drag  = 0.0;
  m_lateral = 0.0;
  m_lift  = 0.0;
  m_xyzMoment  = 0.0;
  m_force  = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::updateOutputFileAero()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  
  // static RealVector ff(0.0, 3);
  // ff += m_frictionForces;
  const CFuint sizeData = 3 + m_force.size() +  m_xyzMoment.size();
  
  CFVec<double> outputData(sizeData); 
  CFVec<double> dragFaeroFxyzMxyz(sizeData);
  
  dragFaeroFxyzMxyz[0] = m_drag;
  dragFaeroFxyzMxyz[1] = m_lateral; 
  dragFaeroFxyzMxyz[2] = m_lift;
  
  CFuint counter = 3;
  for (CFuint i = 0; i < m_force.size(); ++i, ++counter) {
    dragFaeroFxyzMxyz[counter] = m_force[i]; 
  }
  
  for (CFuint i = 0; i < m_xyzMoment.size(); ++i, ++counter) { 
    dragFaeroFxyzMxyz[counter] = m_xyzMoment[i];  
  } 
  
  // here processes must all get here
  const std::string nsp = this->getMethodData().getNamespace();
  MPI_Allreduce(&dragFaeroFxyzMxyz[0], &outputData[0], sizeData, 
		MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));    
  
  if (PE::GetPE().GetRank(nsp) == 0) { 
    // preparation of the output 
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() / 
      boost::filesystem::path(m_nameOutputFileAero); 
    file = Framework::PathAppender::getInstance().appendAllInfo
      (file,false,false,false); 
    
      SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create(); 
      ofstream& fout = fhandle->open(file, ios::app); 
      
      Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive(); 
      
      //Writing the integrated output values
      fout << SubSystemStatusStack::getActive()->getNbIter() << " "
	   << SubSystemStatusStack::getActive()->getCurrentTimeDim() << " "
	   << m_alphadeg << " "
	   << m_betadeg  << " "
	   << outputData << "\n";
      
      // cout << "FRICTION FORCE COEFF = " << ff << endl;
      
      fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::computeSurfaceResiduals()
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  // compute the L2 norm of some quantities of interest
  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
  const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
  if (currTrs->getName() == getTrsList().front()->getName()) {
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      m_valuesMatL2[varID] = 0.;
    } 
  }
  
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) 
  {
    for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
    {
      const CFuint index = m_mapTrsFaceToID.find(currTrs->getLocalGeoID(iFace)*m_nbrFaceFlxPnts+iFlx);
      for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) 
      {
        m_valuesMatL2[varID] += m_valuesMatRes(varID, index)*m_valuesMatRes(varID, index);
      }
    }
  }
  
  cf_assert(m_valuesMatL2.size() > 0);
  
  const std::string nsp = this->getMethodData().getNamespace();
  
  bool writeFile = false;
  if (currTrs->getName() == getTrsList().back()->getName()) {
#ifdef CF_HAVE_MPI
    m_l2Norm = 0.0;
    const CFuint count = m_l2Norm.size();
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    MPI_Allreduce(&m_valuesMatL2[0], &m_l2Norm[0], count, MPI_DOUBLE, MPI_SUM, comm);
    m_valuesMatL2 = m_l2Norm;
#endif
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      //m_valuesMatL2[varID] = (std::abs(m_valuesMatL2[varID]) > 1e-16) ? std::log(std::sqrt(m_valuesMatL2[varID])) : 1e-16;
      m_valuesMatL2[varID] = std::log(std::sqrt(m_valuesMatL2[varID]));
    }
    writeFile = true;
  }
  
  if (PE::GetPE().GetRank(nsp) == 0 && writeFile) {  
    // // preparation of the output
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(m_outputFileConv);
    file = PathAppender::getInstance().appendAllInfo( file, false, false, false );
    
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    
    if (SubSystemStatusStack::getActive()->getNbIter() == 1) {
      ofstream& fout = fhandle->open(file);
      fout << "TITLE = Convergence file for surface quantities\n";
      fout << "VARIABLES = iter ";
      for (CFuint i = dim; i < this->m_varNames.size(); ++i) {
	fout << "Res["<< this->m_varNames[i] << "] ";
      }
      fout << "\n";
      //closing the file
      fhandle->close();
    } 
    
    for (CFuint i = dim; i < this->m_varNames.size(); ++i) {
      cout << "Res["<< this->m_varNames[i] << "] ";
    }
    cout << "\n";
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      cout.precision(8); cout << m_valuesMatL2[varID] << " ";
    }
    cout << endl << endl;
    
    ofstream& fout = fhandle->open(file, ios::app);
    fout << SubSystemStatusStack::getActive()->getNbIter() << " ";
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      fout.precision(8); fout.width(3 + 8); fout << m_valuesMatL2[varID] << " ";
    }
    fout << endl;
    
    //closing the file
    fhandle->close();
  } 
  
  PE::GetPE().setBarrier(nsp);
}
      
//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::initSurfaceResiduals()
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbVariables = m_varNames.size() - dim;
  if (nbVariables < 1) 
  {
    CFLog(WARN, "AeroForcesFR::initSurfaceResiduals() needs nbVariables > 0 => you must use derived class\n");
    cf_always_assert(nbVariables > 0);
  }
  m_valuesMatL2.resize(nbVariables);
  m_l2Norm.resize(nbVariables);
  
  if (m_mapTrsFaceToID.size() == 0) 
  {
    vector< SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
    
    // empty loop just to count faces
    CFuint totalNbFaces = 0;
    for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) 
    {
      totalNbFaces += trsList[iTRS]->getLocalNbGeoEnts();
    }
    
    if (totalNbFaces > 0) 
    {
      const CFuint totalNbFlxPnts = totalNbFaces*m_nbrFaceFlxPnts;

      m_mapTrsFaceToID.reserve(totalNbFlxPnts);
      m_valuesMat.resize(nbVariables, totalNbFlxPnts);
      m_valuesMatRes.resize(nbVariables, totalNbFlxPnts);
      
      CFuint index = 0;  
      for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) 
      {
	SafePtr<TopologicalRegionSet> trs = trsList[iTRS];
	const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) 
	{
	  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
	  {
	    // CFLog(INFO, "faceID = " << trs->getLocalGeoID(iFace) << ", index = " << index << "\n");
	    const CFuint flxPntIdx = (trs->getLocalGeoID(iFace))*m_nbrFaceFlxPnts+iFlx;

	    m_mapTrsFaceToID.insert(flxPntIdx, index++);
	  }
	}
      }
      
      m_mapTrsFaceToID.sortKeys();
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

void AeroForcesFR::reorderOutputFileWall()
{      
  if (m_reorderWallData && PhysicalModelStack::getActive()->getDim() == DIM_2D) {
    SafePtr<TopologicalRegionSet> currTrs = this->getCurrentTRS();
    
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(this->m_nameOutputFileWall + currTrs->getName());
    file = Framework::PathAppender::getInstance().appendAllInfo  
      (file,this->m_appendIter,this->m_appendTime,false);   
        
    const vector<vector<CFuint> >&  trsInfo = MeshDataStack::getActive()->getTotalTRSInfo();
    const vector<std::string>& trsNames = MeshDataStack::getActive()->getTotalTRSNames();
    
    CFuint nbTrsFaces = 0;
    const CFuint nbTRSs = trsInfo.size();
    for(CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      if (trsNames[iTRS] == currTrs->getName()) {
	const CFuint nbTRsInTRS = trsInfo[iTRS].size();
	for(CFuint tr = 0; tr < nbTRsInTRS; ++tr) {
	  nbTrsFaces += trsInfo[iTRS][tr];
	}
      }
    }
    // cout << "nbTrsFaces = " << nbTrsFaces << endl;
    
    const CFuint nbVars = m_varNames.size(); 
    const CFuint nbFlxPnts = nbTrsFaces*m_nbrFaceFlxPnts;
    CFMap<CFreal,CFuint> mmap(nbFlxPnts);
    vector<vector<CFreal> > m(nbVars);
    
    // read the file
    const std::string nsp = this->getMethodData().getNamespace();
    if (PE::GetPE().GetRank(nsp) == 0) {
      Common::SelfRegistPtr<Environment::FileHandlerInput> fhandleIn =
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
      ifstream& fin = fhandleIn->open(file);
      
      const string lastVariableName = m_varNames.back();
      string tmp = "";
      while (tmp != lastVariableName) {
	fin >> tmp;
	CFLog(VERBOSE, "AeroForcesFR::reorderOutputFileWall() => reading " << tmp << "\n");
      }
            
      for (CFuint i = 0; i < nbVars; ++i) {
	m[i].resize(nbFlxPnts);
      }
      
      for (CFuint iFlx = 0; iFlx < nbFlxPnts; ++iFlx) {
	for (CFuint i = 0; i < nbVars; ++i) {
	  fin >> m[i][iFlx];
	}
	mmap.insert(m[0][iFlx],iFlx); //x - iFace
      }
      mmap.sortKeys();
      fin.close();
    }
    
    // write the file header 
    prepareOutputFileWall();
    
    // write te reordered solution data   
    if (PE::GetPE().GetRank(nsp) == 0) {
      SelfRegistPtr<Environment::FileHandlerOutput> fhandleOut =
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
      
      // append to the existing file 
      ofstream& fout = fhandleOut->open(file, ios::app);
      
      for (CFuint iFlx = 0; iFlx < nbFlxPnts; ++iFlx) {
	const CFuint idx = mmap[iFlx];
	for (CFuint i = 0; i < nbVars; ++i) {
	  fout << m[i][idx] << " ";
	}
	fout << "\n";
      }
      fout.close();
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////
     
    } // namespace AeroCoeff

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
