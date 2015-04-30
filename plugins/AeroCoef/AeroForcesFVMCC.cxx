#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"

#include "AeroCoef/AeroCoef.hh"
#include "AeroCoef/AeroForcesFVMCC.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "FiniteVolume/CellCenterFVM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<AeroForcesFVMCC, DataProcessingData, AeroCoefModule>
aeroForcesFVMCCProvider("AeroForcesFVMCC");

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::defineConfigOptions(Config::OptionList& options)
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

AeroForcesFVMCC::AeroForcesFVMCC(const std::string& name) :
  DataProcessingCom(name),
  m_faceTrsGeoBuilder(),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_gstates("gstates"),
  socket_normals("normals"),
  m_updateVarSet(),
  m_fvmccData(CFNULL),
  m_mapTrsFaceToID(),
  m_dataState(),
  m_coord(),
  m_frictionForces(),
  m_unitNormal(),
  m_normal(),
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
  m_avState(CFNULL),
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
  m_beta(0.)
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
  
  _outputFileConv = "convergence-surf.plt";
  setParameter("OutputFileConv",&_outputFileConv);
}

//////////////////////////////////////////////////////////////////////////////

AeroForcesFVMCC::~AeroForcesFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
AeroForcesFVMCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
AeroForcesFVMCC::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::setup()
{
  if (Common::PE::GetPE().IsParallel()) {
    // there is only one file for the aerodynamic coefficients
    m_nameOutputFileAero = m_nameOutputFileAero + "-P0";
    
    // there is only one file for the wall data
    m_nameOutputFileWall = m_nameOutputFileWall + "-P0";
  }

  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());

  m_fvmccData = fvmcc->getData();
  m_updateVarSet = m_fvmccData->getUpdateVar().d_castTo<EulerVarSet>();
  cf_assert(m_updateVarSet.isNotNull());
  m_updateVarSet->getModel()->resizePhysicalData(m_dataState);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_coord.resize(dim);
  m_frictionForces.resize(dim);
  m_unitNormal.resize(dim);
  m_normal.resize(dim);
  m_xyzMoment.resize((dim == DIM_2D) ? 1 : dim);

  m_v01.resize(dim);
  m_v02.resize(dim);
  m_vCross0102.resize(dim);
  m_v12.resize(dim);
  m_xg.resize(dim);
  m_midFaceNode.resize(dim);
  m_rotMat.resize(dim,dim);
  m_xyzForce.resize(dim);
  m_force.resize(dim);
  m_aeroForce.resize(dim);
  
  m_avState = new State();

  cf_always_assert(m_uInf > 0.);
  cf_always_assert(m_rhoInf > 0.);
  cf_always_assert(m_pInf > 0.);

  if (m_gravityCenter.size() == 0) {
    m_gravityCenter.resize(dim, 0.0);
    if (PE::GetPE().GetRank() == 0) {
      cout << "#### WATCH OUT: gravity center not specified !!! assuming X=0 #####" << endl;
    }
  }

  for (CFuint i =0; i < dim; ++i) {
    m_xg[i]= m_gravityCenter[i];
  }

  // pre-computation of the total wet surface
  if (dim == DIM_3D) {
    computeWetSurface();
  }
  else {
    cf_always_assert(m_refLength2D > 0.);
    m_wetSurface = m_refLength2D;
  }

  // reference area is set to be the wet area if it is not specified by the user
  m_refArea = (m_refArea > 0.) ? m_refArea : m_wetSurface;
  
  m_outputFileAeroPrepared = false;
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::computeWetSurface()
{
  vector<SafePtr<TopologicalRegionSet> >& trsList = getTrsList();

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = m_fvmccData->getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;

  CFreal wSurface = 0.0;
  for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
    geoData.trs = trsList[iTRS];
    
    const CFuint nbTrsFaces = trsList[iTRS]->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity* const currFace = geoBuilder->buildGE();
      
      if (currFace->getState(0)->isParUpdatable()) {
	// AL: carefull here, you should never ask a face to compute the volume (in this case it should be fine)
	wSurface += currFace->computeVolume();
      }
      
      // release the geometric entity
      geoBuilder->releaseGE();
    }
  }
  
  // sum all the contributions on this TRS from all processors
  MPI_Allreduce(&wSurface, &m_wetSurface, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator());     
  
  cf_always_assert(m_wetSurface > 0.0);
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::configure ( Config::ConfigArgs& args )
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

void AeroForcesFVMCC::setFunction()
{
  // some sanity checks on alpha
  std::vector<std::string> m_functionAlphaDef = Common::StringOps::getWords(m_functionAlpha);
  if(m_functionAlphaDef.size() != 1)
    throw BadValueException (FromHere(),"AeroForcesFVMCC::setFuntion(): alpha function wrongly defined.");

  m_functionAlphaParser.Parse(m_functionAlpha, m_vars);

  if (m_functionAlphaParser.ErrorMsg() != 0) {
    std::string msg("ParseError in AeroForcesFVMCC::setFuntion(): ");
    msg += std::string(m_functionAlphaParser.ErrorMsg());
    msg += " Function: " + m_functionAlpha;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }

  // some sanity checks on beta
  std::vector<std::string> m_functionBetaDef = Common::StringOps::getWords(m_functionBeta);
  if(m_functionBetaDef.size() != 1)
    throw BadValueException (FromHere(),"AeroForcesFVMCC::setFuntion(): alpha function wrongly defined.");

  m_functionBetaParser.Parse(m_functionBeta, m_vars);

  if (m_functionBetaParser.ErrorMsg() != 0) {
    std::string msg("ParseError in AeroForcesFVMCC::setFuntion(): ");
    msg += std::string(m_functionBetaParser.ErrorMsg());
    msg += " Function: " + m_functionBeta;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::executeOnTrs()
{
  CFAUTOTRACE;
  
  if (PE::GetPE().GetRank() == 0) {
    CFLog(VERBOSE, "AeroForcesFVMCC::executeOnTrs() START\n");
  }
	
  //  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute the value of the angle
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionAlphaParser.Eval(m_eval);
  m_betadeg  = m_functionBetaParser.Eval(m_eval);

  // transform into Radiants
  m_alpha = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;
  m_beta  = m_betadeg*MathTools::MathConsts::CFrealPi()/180;

  DataHandle< CFreal> normals = socket_normals.getDataHandle();
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

  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
  CFLog(VERBOSE, "NavierStokesSkinFrictionHeatFluxCC::computeValues() on TRS "
	<< currTrs->getName() << "\n");

  Common::SafePtr<GeometricEntityPool<FaceCellTrsGeoBuilder> > 
    geoBuilder = m_fvmccData->getFaceCellTrsGeoBuilder();

  SafePtr<FaceCellTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  
  geoData.faces = currTrs;
  geoData.isBFace = true;
  geoData.cells = MeshDataStack::getActive()->getTrs("InnerCells");
  geoData.allCells = true;
  
  const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
  
  if (m_valuesMatRes.size() == 0) {initSurfaceResiduals();}
  prepareOutputFileAero();
  prepareOutputFileWall();
  
  // this needs to be updated with the latest values
  m_fvmccData->getPolyReconstructor()->computeGradients();
  
  for (m_iFace = 0; m_iFace < nbTrsFaces; ++m_iFace) {
    CFLogDebugMed("iFace = " << m_iFace << "\n");

    // build the GeometricEntity
    geoData.idx = m_iFace;
    m_currFace = geoBuilder->buildGE();
    m_fvmccData->getCurrentFace() = m_currFace;

    // compute the face normal
    const CFuint faceID = m_currFace->getID();
    const CFuint startID = faceID*dim;
    // face normal is pointing outward from the computational domain
    // which means that is pointing into the body
    // we consider the opposite normal, pointing inward the domain
    for (CFuint i = 0; i < dim; ++i) {
      m_normal[i] = -normals[startID + i];// this IS NOT a unit normal : it is \vec{1}_n \delta S
    }
    m_unitNormal = m_normal*(1./m_normal.norm2());// this IS the unit normal
    
    // this is fundamental because some algorithms for computing numerical derivatives need this
    m_fvmccData->getUnitNormal() = -m_unitNormal;  
    
    computeWall();
    computeAero();
    
    // release the geometric entity
    geoBuilder->releaseGE();
  }
  
  // all data are written on file at once to ease the parallel writing
  updateOutputFileAero(); 
  updateOutputFileWall();
  reorderOutputFileWall();
  computeSurfaceResiduals(); 
  
  if (PE::GetPE().GetRank() == 0) {
    CFLog(VERBOSE, "AeroForcesFVMCC::executeOnTrs() END\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::computeWall()
{
  CFAUTOTRACE;

  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  const vector<Node*>& nodesInFace = *m_currFace->getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();

  *m_avState = 0.0;
  m_coord = 0.0;
  for (CFuint i = 0; i < nbNodesInFace; ++i) {
    *m_avState += nstates[nodesInFace[i]->getLocalID()];
    m_coord += *nodesInFace[i];
  }
  *m_avState /= nbNodesInFace;
  m_coord /= nbNodesInFace;

  m_updateVarSet->computePhysicalData(*m_avState, m_dataState);
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::updateOutputFileWall()
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

void AeroForcesFVMCC::computeAero()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // this is redundant but allows to adapt all the subclasses
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  const vector<Node*>& nodesInFace = *m_currFace->getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();

  // compute all physical values corresponding to the given states
  *m_avState = 0.0;
  m_midFaceNode = 0.0;
  for (CFuint i = 0; i < nbNodesInFace; ++i) {
    *m_avState += nstates[nodesInFace[i]->getLocalID()];
    m_midFaceNode += *nodesInFace[i];
  }
  *m_avState /= nbNodesInFace;
  m_midFaceNode /= nbNodesInFace;
  
  cf_assert(m_updateVarSet.isNotNull());
  m_updateVarSet->computePhysicalData(*m_avState, m_dataState);

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
  m_xyzForce  = (-m_Cp/m_refArea)*m_normal + m_frictionForces;
  m_aeroForce = m_rotMat*m_xyzForce; // [D C L]^T = {rotation matrix} * [X Y Z]^T

  // vector between center of gravity and centre of the face where the force is applied
  m_v12 = m_midFaceNode - m_xg;
  
  // assemble the contribution of this face only if the internal state is parallel updatable 
  if (m_currFace->getState(0)->isParUpdatable()) {
    if (dim == DIM_2D) {
      // this is the moment with respect to Z axis (Mz > 0 => pitching down)
      m_xyzMoment[0] += (m_v12[XX]*m_xyzForce[YY] - m_v12[YY]*m_xyzForce[XX])/m_refLength2D;
    }
    else {
      MathFunctions::crossProd(m_v12, m_xyzForce, m_vCross0102);
      m_xyzMoment += m_vCross0102/m_refLength2D;
    }
    
    m_force    += m_xyzForce;
    m_drag     += m_aeroForce[XX];
    m_lateral  += (dim == DIM_2D) ? 0. : m_aeroForce[YY];
    m_lift     += m_aeroForce[(dim == DIM_2D) ? YY : ZZ];
  }
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::unsetup()
{
  deletePtr(m_avState);
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::prepareOutputFileWall()
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

void AeroForcesFVMCC::prepareOutputFileAero()
{
  CFAUTOTRACE;
  
  if (PE::GetPE().GetRank() == 0 && !m_outputFileAeroPrepared) {
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

  m_drag  = 0.;
  m_lateral = 0.;
  m_lift  = 0.;
  m_xyzMoment  = 0.;
  m_force  = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::updateOutputFileAero()
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
  dragFaeroFxyzMxyz[2] = m_drag;
  
  CFuint counter = 3;
  for (CFuint i = 0; i < m_force.size(); ++i, ++counter) {
    dragFaeroFxyzMxyz[counter] = m_force[i]; 
  }
  
  for (CFuint i = 0; i < m_xyzMoment.size(); ++i, ++counter) { 
    dragFaeroFxyzMxyz[counter] = m_xyzMoment[i];  
  } 
  
  // here processes must all get here 
  MPI_Allreduce(&dragFaeroFxyzMxyz[0], &outputData[0], sizeData, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator());    
  
  if (PE::GetPE().GetRank() == 0) { 
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

void AeroForcesFVMCC::computeSurfaceResiduals()
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
  
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    const CFuint index = m_mapTrsFaceToID.find(currTrs->getLocalGeoID(iFace));
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      m_valuesMatL2[varID] += m_valuesMatRes(varID, index)*m_valuesMatRes(varID, index);
    }
  }
  
  cf_assert(m_valuesMatL2.size() > 0);
  
  bool writeFile = false;
  if (currTrs->getName() == getTrsList().back()->getName()) {
#ifdef CF_HAVE_MPI
    m_l2Norm = 0.0;
    const CFuint count = m_l2Norm.size();
    MPI_Comm comm = PE::GetPE().GetCommunicator();
    MPI_Allreduce(&m_valuesMatL2[0], &m_l2Norm[0], count, MPI_DOUBLE, MPI_SUM, comm);
    m_valuesMatL2 = m_l2Norm;
#endif
    for (CFuint varID = 0; varID < m_valuesMatL2.size(); ++varID) {
      m_valuesMatL2[varID] = (std::abs(m_valuesMatL2[varID]) > 1e-16) ? std::log(std::sqrt(m_valuesMatL2[varID])) : 1e-16;
    }
    writeFile = true;
  }
  
  if (PE::GetPE().GetRank() == 0 && writeFile) {  
    // // preparation of the output
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(_outputFileConv);
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
  
  PE::GetPE().setBarrier();
}
      
//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::initSurfaceResiduals()
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbVariables = m_varNames.size() - dim;
  cf_always_assert(nbVariables > 0);
  m_valuesMatL2.resize(nbVariables);
  m_l2Norm.resize(nbVariables);
  
  if (m_mapTrsFaceToID.size() == 0) {
    vector< SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
    
    // empty loop just to count faces
    CFuint totalNbFaces = 0;
    for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
      totalNbFaces += trsList[iTRS]->getLocalNbGeoEnts();
    }
    
    if (totalNbFaces > 0) {
      m_mapTrsFaceToID.reserve(totalNbFaces);
      m_valuesMat.resize(nbVariables, totalNbFaces);
      m_valuesMatRes.resize(nbVariables, totalNbFaces);
      
      CFuint index = 0;  
      for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
	SafePtr<TopologicalRegionSet> trs = trsList[iTRS];
	const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  // CFLog(INFO, "faceID = " << trs->getLocalGeoID(iFace) << ", index = " << index << "\n");
	  m_mapTrsFaceToID.insert(trs->getLocalGeoID(iFace), index++);
	}
      }
      cf_assert(index == totalNbFaces);
      
      m_mapTrsFaceToID.sortKeys();
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////

void AeroForcesFVMCC::reorderOutputFileWall()
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
    CFMap<CFreal,CFuint> mmap(nbTrsFaces);
    vector<vector<CFreal> > m(nbVars);
    
    // read the file
    if (PE::GetPE().GetRank() == 0) {
      Common::SelfRegistPtr<Environment::FileHandlerInput> fhandleIn =
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
      ifstream& fin = fhandleIn->open(file);
      
      const string lastVariableName = m_varNames.back();
      string tmp = "";
      while (tmp != lastVariableName) {
	fin >> tmp;
	CFLog(VERBOSE, "AeroForcesFVMCC::reorderOutputFileWall() => reading " << tmp << "\n");
      }
            
      for (CFuint i = 0; i < nbVars; ++i) {
	m[i].resize(nbTrsFaces);
      }
      
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	for (CFuint i = 0; i < nbVars; ++i) {
	  fin >> m[i][iFace];
	}
	mmap.insert(m[0][iFace],iFace); //x - iFace
      }
      mmap.sortKeys();
      fin.close();
    }
    
    // write the file header 
    prepareOutputFileWall();
    
    // write te reordered solution data   
    if (PE::GetPE().GetRank() == 0) {
      SelfRegistPtr<Environment::FileHandlerOutput> fhandleOut =
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
      
      // append to the existing file 
      ofstream& fout = fhandleOut->open(file, ios::app);
      
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	const CFuint idx = mmap[iFace];
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
