#include "Common/PE.hh"

#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/PhysicalModel.hh"

#include "AeroCoef/AeroCoefFS.hh"
#include "AeroCoef/NavierStokes2DConsComputeAeroHO.hh"

#include "FluctSplit/FluctuationSplit.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokes2DConsComputeAeroHO,
                      DataProcessingData,
                      AeroCoefFSModule>
aNavierStokes2DConsComputeAeroHOProvider("NavierStokes2DConsComputeAeroHO");

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Definition of the function defining the angle between the flow and the x-axis.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coefficients.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFreal >("UInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("RhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("PInf","Pressure at infinity.");
   options.addConfigOption< std::vector<CFreal> >("Xref","Reference point x to which the pressure momentum is computed");
   options.addConfigOption< CFreal >("TInf","Temperarure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
   options.addConfigOption< bool >("IsViscous","Is the computation viscous?");
   options.addConfigOption< bool >("NACA0012_Exactnormals","Do you want to use the exact normals of the NACA0012 to compute the aerocoefficient?");

}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DConsComputeAeroHO::NavierStokes2DConsComputeAeroHO(const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_faceNeighCell("faceNeighCell"),
  socket_normals("normals"),
  m_updateVarSet(),
  m_vars("t"),
  m_eval(0.0,1),
  actual_lift(),
  actual_drag(),
  actual_momentum(),
  m_alphadeg(0.),
  m_alpharad(0.),
  m_dimState_face(),
  m_dimState_cell(),
  _cells(CFNULL),
  m_diffVar(),
  _gradients(),
  _avValues(),
  m_fsData(CFNULL),
  m_nodal_cf(),
  m_numb_node(),
  m_dataState(),
  m_qd0(),
  m_qd1(),
  m_wqd(),
  m_qdstates(),
  m_qdnodes(),
  m_mappedCoord(),
  m_mappedCoord_line(),
  m_physicalCoord(),
  m_physicalCoord_line(),
  m_gradients(),
  m_tau(),
  m_QdCf(),
  m_QdCmf()
{
  m_fileWall = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFileWall = "Wall";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_saveRateWall = 1;
  setParameter("SaveRateWall",&m_saveRateWall);

  m_function = std::string("0.");
  setParameter("Alpha",&m_function);

  m_nameOutputFileAero = "AeroCoef.plt";
  setParameter("OutputFileAero",&m_nameOutputFileAero);

  m_saveRateAero = 1;
  setParameter("SaveRateAero",&m_saveRateAero);

  m_uInf = 0.;
  setParameter("UInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("RhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("PInf",&m_pInf);
  m_tInf = 0.;
  setParameter("TInf",&m_tInf);

   m_Xref.resize(2);
   m_Xref[XX] = 0.25;
   m_Xref[YY] = 0.0;
   setParameter("Xref",&m_Xref);


  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

  m_isViscous = true;
  setParameter("IsViscous",&m_isViscous);

  m_exact_norm = false;
  setParameter("NACA0012_Exactnormals",&m_exact_norm);

}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DConsComputeAeroHO::~NavierStokes2DConsComputeAeroHO()
{
	if (isSetup()) unsetup();
}
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokes2DConsComputeAeroHO::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_faceNeighCell);
result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NavierStokes2DConsComputeAeroHO::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::setup()
{

  actual_lift = 0.;
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluctuationSplit> fs = spaceMethod.d_castTo<FluctuationSplit>();

  m_fsData = fs->getData();
  try
  {
    m_diffVar = m_fsData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  }
  catch (FailedCastException& e)
  {
    m_diffVar = NULL;
  }

  m_updateVarSet->getModel()->resizePhysicalData(m_dataState);
  const RealVector& refValues = m_updateVarSet->getModel()->getReferencePhysicalData();
  _cells = MeshDataStack::getActive()->getTrs("InnerCells");
  if (refValues[EulerTerm::V] == 0.   ||
      refValues[EulerTerm::RHO] == 0. ||
      refValues[EulerTerm::P] == 0.   ||
     (refValues[EulerTerm::VX] == 0. && refValues[EulerTerm::VY] == 0.))
  {
    std::string msg = std::string("Some needed reference values not set.");
    throw BadValueException (FromHere(),msg);
  }

  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_gradients.resize(PhysicalModelStack::getActive()->getNbEq());
  for (CFuint i = 0; i< nbEqs; ++i)
  {
    m_gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  _avValues.resize(PhysicalModelStack::getActive()->getNbEq());
  prepareOutputFileAero();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_physicalCoord.resize(dim);
  m_mappedCoord.resize(dim);
  m_physicalCoord_line.resize(dim-1);
  m_mappedCoord_line.resize(dim-1);

  CFuint nbqd = 3;
  //m_qdNormal.resize(nbqd);
  node_px.resize(nbqd);
  node_py.resize(nbqd);
  node_m.resize(nbqd);
  m_tau.resize(nbqd);
  m_QdCf.resize(nbqd);
  m_QdCmf.resize(nbqd);

  for (CFuint i = 0; i< nbqd; ++i)
  {
//    m_qdNormal[i].resize(dim);
    m_QdCf[i].resize(dim);
  }

  // Set up the parameters of the quadrature rule'
  m_qd0.resize(nbqd); // quadrature points per face
  m_qd1.resize(nbqd); // quadrature points per face

  const CFreal s  = std::sqrt( 0.6 );
  const CFreal a0 = -s*0.5;
  const CFreal a1 = s*0.5;

  m_qd0[0] = a0;  m_qd1[0] = a1;
  m_qd0[1] = a1;  m_qd1[1] = a0;
  m_qd0[2] = .0;  m_qd1[2] = .0;

  m_wqd.resize(nbqd);
  m_wqd[0] = 5.0/9.0;
  m_wqd[1] = 5.0/9.0;
  m_wqd[2] = 8.0/9.0;

  m_qdstates.resize(nbqd);
  m_qdnodes.resize(nbqd);
  for (CFuint iState = 0; iState < nbqd; ++ iState){
    m_qdstates[iState] = new State();
    m_qdnodes[iState] = new Node();
    m_qdstates[iState]->setSpaceCoordinates(m_qdnodes[iState]);
  }



    m_qdNormal.resize(dim);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<Physics::NavierStokes::EulerVarSet>();

  // Create the socket source for the boundaryNormals for each TRS to which this command applies
  // This is because we have different normals in FiniteVolume and FluctSplit / FiniteElement
  const CFuint nbTrs = getTrsNames().size();

  for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs)
  {
    std::string socketName = getTrsName(iTrs) + "-boundaryNormals";
    m_sockets.createSocketSource<const CFreal*>(socketName);
  }

  // set the function that computes the alpha
  try
  {
    setFunction();
  }
  catch (Common::Exception& e)
  {
    CFout << e.what() << "\n";
    throw;
  }

  CFLog(WARN, "PInf   : " << m_pInf << "\n");
  CFLog(WARN, "RhoInf : " << m_rhoInf << "\n");
  CFLog(WARN, "UInf   : " << m_uInf << "\n");

  if (MathChecks::isZero(m_pInf))   throw BadValueException (FromHere(),"pInf, the pressure at infinity, is zero");
  if (MathChecks::isZero(m_rhoInf)) throw BadValueException (FromHere(),"rhoInf, the density at infinity, is zero");
  if (MathChecks::isZero(m_uInf))   throw BadValueException (FromHere(),"uInf, the velocity at infinity, is zero");


  CFLog(WARN, "PInf   : " << m_pInf << "\n");
  CFLog(WARN, "RhoInf : " << m_rhoInf << "\n");
  CFLog(WARN, "UInf   : " << m_uInf << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokes2DConsComputeAeroHO::setFuntion(): alpha function wrongly defined.");

  m_functionParser.Parse(m_function, m_vars);

  if (m_functionParser.ErrorMsg() != 0)
  {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(m_functionParser.ErrorMsg());
    msg += " Function: " + m_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::executeOnTrs()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();

  const CFuint iter = ssys_status->getNbIter();

  // compute global integrated coefficients
  if(!(iter % m_saveRateAero)) { computeAero(); }

  // compute wall distribution
  if(!(iter % m_saveRateWall)) { computeWall(); }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file is opened here

  const CFreal R     = m_updateVarSet->getModel()->getRdim();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();
  RealVector dimState;
  dimState.resize(nbEqs);

   //First we compute Cf if the computation is viscous
   if (m_isViscous){

   // get the datahandle containing the boundary Normals
  const std::string TRSName = getCurrentTRS()->getName();
  const std::string socketName = TRSName + "-boundaryNormals";
  DataHandle<const CFreal*> boundaryNormals = m_sockets.
  getSocketSource<const CFreal*>(socketName)->getDataHandle();

    //  building the TRS for faces
    // builder for standard TRS GeometricEntity's taht will be used for the faces
    Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
    geoBuilderFace.setup();
    StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();
    geoDataFace.trs = getCurrentTRS();
    const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();

    // builder for standard TRS GeometricEntity's
    Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
    geoBuilderCell.setup();

    StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

    // nodal_cf contains the value of each Cf at each states of the TRS
    // and numb_node its local ID.
    // we need both because here we get the faces of the TRS to compute the Cf
    // but, after to write the wall aerocoefircient we get the states of the TRS

    m_nodal_cf.resize(statesIdx->size(),2);
    m_numb_node.resize(statesIdx->size());
    RealVector normal(0.0, DIM_2D);

    // n is the position of the last Cf stored in m_nodal_cf
    CFuint n = 0;

    for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
    {
      geoDataFace.idx = iFace;
      GeometricEntity & currFace = *geoBuilderFace.buildGE();
      CFuint FaceID = currFace.getID();
      // get the face normal

      for(CFuint iDim = 0; iDim < DIM_2D; iDim++){
        normal[iDim] = -(boundaryNormals[iFace])[iDim];
      }

      normal.normalize();


      std::vector<State*>& states = *(currFace.getStates());
      cf_assert(states.size() >= 2);
      const CFuint  nbNodesInFace = states.size();

      // handle to the neighbor cell
      DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
        faceNeighCell = socket_faceNeighCell.getDataHandle();

      const CFuint cellID = faceNeighCell[FaceID].first;

      geoDataCell.trs = faceNeighCell[FaceID].third;
      geoDataCell.idx = cellID;
      GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

      CFreal tau = 0.0;
      vector<CFreal > u;
      vector<CFreal > v;
      std::vector<RealVector> m_coords(1);
      m_coords[0].resize(2);

      const CFreal faceArea = currFace.computeVolume();

      vector<Framework::State*> * m_states= neighborCell.getStates();
      vector<Framework::State*>& states_cell = (*m_states);
      const CFuint nbStatesInCell = m_states->size();

      m_dimState_cell.resize(nbStatesInCell);
    for (CFuint s = 0; s < nbStatesInCell; ++s){
      m_dimState_cell[s].resize(nbEqs);
    }
    for (CFuint s = 0; s < nbStatesInCell; ++s){
      m_updateVarSet->setDimensionalValues(*states_cell[s],m_dimState_cell[s]);
    }

      for (CFuint iState = 0; iState < nbNodesInFace; ++iState){
        m_coords[0] = states[iState]->getCoordinates();
        m_updateVarSet->setDimensionalValues(*states[iState],dimState);
      // Compute the shape function gradient at the projected point
      const std::vector<RealMatrix> cellGradients= neighborCell.computeSolutionShapeFunctionGradients(m_coords);

      const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

      // Compute the gradients variables
      RealVector gradients_u(0., 2);
      RealVector gradients_v(0., 2);

      for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){

        for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){

          (*m_gradients[iEq])[iDim] = 0.0;
        }
      }

      for (CFuint iNode = 0; iNode < nbStatesInCell ; ++iNode){
        for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){

          for (CFuint iEq = 0; iEq < nbEqs ; ++iEq)
            (*m_gradients[iEq])[iDim] += (cellGradients[0])(iNode,iDim) *m_dimState_cell[iNode][iEq];
          }

      }

      gradients_u[XX] = (1.0/dimState[0])*((*m_gradients[1])[XX] - dimState[1]*(*m_gradients[0])[XX]);
      gradients_u[YY] = (1.0/dimState[0])*((*m_gradients[1])[YY] - dimState[1]*(*m_gradients[0])[YY]);

      gradients_v[XX] = (1.0/dimState[0])*((*m_gradients[2])[XX] - dimState[2]*(*m_gradients[0])[XX]);
      gradients_v[YY] = (1.0/dimState[0])*((*m_gradients[2])[YY] - dimState[2]*(*m_gradients[0])[YY]);



      CFreal mu = 0.;
      if (m_diffVar.isNotNull()) mu = m_diffVar->getDynViscosity(dimState, _gradients);

      // Need to see if we are on the bottom part or on the top.
      // This is important to see if the axes refered to the boundary is well oriented

      if (normal[YY] > 0.0){
        tau = mu*( normal[YY]*( gradients_u[YY]*normal[YY] - gradients_v[YY]*normal[XX] )
                      + normal[XX]*(gradients_u[XX]*normal[YY] - gradients_v[XX]*normal[XX] ) );
      }

      else {
       tau = mu*( normal[YY]*(-gradients_u[YY]*normal[YY] + gradients_v[YY]*normal[XX] )
                           + normal[XX]*(-gradients_u[XX]*normal[YY] + gradients_v[XX]*normal[XX] ) );
      }


      // m_nodal_cf contain the cf of each node of the TRS.
      // since the gradients are computed per elements we need to do an average
      // between the cf of faces
      // and we store the localID of the node in a vector. Then it is necessary to know if m_nodal_cf
      // has already a part of the value of the node

        CFuint marker = 0;
        for (CFuint j = 0; j< n; ++j){
          if(m_numb_node[j] == states[iState]->getLocalID()){
            m_nodal_cf(j,0) += tau*faceArea;
            m_nodal_cf(j,1) += faceArea;
            marker = 1;
          }
        }
        if (marker == 0){
          m_nodal_cf(n,0) += tau*faceArea;
          m_nodal_cf(n,1) += faceArea;
          m_numb_node[n] =  states[iState]->getLocalID();
          n ++;
        }


    }
      geoBuilderFace.releaseGE();
      geoBuilderCell.releaseGE();
    }
    //adimensionalize nodal_cf and devid by the coef for average

    for ( CFuint i = 0; i < statesIdx->size(); ++i){
        m_nodal_cf(i,0) /= m_nodal_cf(i,1)*(0.5 * m_rhoInf * m_uInf * m_uInf);

}
  }

  // We really start with the wtritting of the wall aerocoefficients
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState){
    const State& state = *(states[(*statesIdx)[iState]]);
    m_updateVarSet->setDimensionalValues(state,dimState);

    const CFreal rho  = dimState[0];
    const CFreal rhoU = dimState[1];
    const CFreal rhoV = dimState[2];
    const CFreal rhoE = dimState[3];

    const CFreal invRho = 1./rho;

    const CFreal rhoU2  = rhoU*rhoU;
    const CFreal rhoV2  = rhoV*rhoV;
    const CFreal U2  = rhoU2*invRho*invRho;
    const CFreal V2  = rhoV2*invRho*invRho;

    const CFreal rhoK2  = 0.5*rho*(U2 + V2);

    // calculate pressure in dimensional values
    const CFreal P = gammaMinus1*(rhoE - rhoK2);

    // calculate pressure coefficient
    const CFreal Cp = (P-m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

    // calculate temperature in dimensional values
    const CFreal T = P / (rho*R);

    //Find the cf corresponding to this sate
    CFreal Cf = 0.;

    //Find the cf corresponding to this sate
    if (m_isViscous){
      for (CFuint i = 0; i < statesIdx->size(); ++ i){
        if (m_numb_node[i] == state.getLocalID())
          Cf = m_nodal_cf(i,0);
      }
      fout
      << (state.getCoordinates())[XX]   << " "
      << (state.getCoordinates())[YY]   << " "
      << (state.getCoordinates())[ZZ]   << " "
      << Cp     << " "
      << Cf     << " "
//      << Mach   << " "
      << P      << " "
      << T      << " "
      << rho    << "\n";
    }
    else{
    fout
    << (state.getCoordinates())[XX]   << " "
    << (state.getCoordinates())[YY]   << " "
    << (state.getCoordinates())[ZZ]   << " "
    << Cp     << " "
//    << Mach   << " "
    << P      << " "
    << T      << " "
    << rho    << "\n";
    }

  }


  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::computeAero()
{
  CFAUTOTRACE;

  // get the datahandle containing the boundary Normals
  const std::string TRSName = getCurrentTRS()->getName();
  const std::string socketName = TRSName + "-boundaryNormals";
  DataHandle<const CFreal*> boundaryNormals = m_sockets.getSocketSource<const CFreal*>(socketName)->getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  cf_assert ( dim == DIM_2D );

  const CFreal gammaMinus1 = m_updateVarSet->getModel()->getGamma() - 1.;
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  vector < CFreal > u;
  vector < CFreal > v;
  RealVector qdValues;
  qdValues.resize(nbEqs);

  RealVector normal(0.0, dim);
  RealVector Cp(0.0, dim);
  RealVector faceCp(0.0, dim);
  RealVector faceCf(0.0, dim);
  CFreal faceCmf = 0.0;
  RealVector Cf(0.0, dim);
  CFreal Cmp = 0.0;
  CFreal Cmf = 0.0;
  CFreal m;

  //  building the TRS for faces
  // builder for standard TRS GeometricEntity's taht will be used for the faces
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
  geoBuilderFace.setup();
  StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();
  geoDataFace.trs = getCurrentTRS();
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  cf_assert(boundaryNormals.size() == nbFaces);

  // builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();

  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    geoDataFace.idx = iFace;

     GeometricEntity & currFace = *geoBuilderFace.buildGE();
    // get the face normal
    const CFuint FaceID = currFace.getID();
    std::vector<State*>& states = *(currFace.getStates());
    std::vector<Node*>& nodes = *(currFace.getNodes());
    const CFuint  nbStatesInFace = states.size();
    m_dimState_face.resize(nbStatesInFace);
     for (CFuint s = 0; s < nbStatesInFace; ++s){
       m_dimState_face[s].resize(nbEqs);
     }
    for (CFuint s = 0; s < nbStatesInFace; ++s){
      m_updateVarSet->setDimensionalValues(*states[s],m_dimState_face[s]);
    }
    cf_assert(nodes.size() == 2);
    RealVector shapeFunctionValues(nbStatesInFace);

    // Check if face is owned by one or more processors and compute the weight accordingly
    CFreal weight = 1.0;
    if (!(states[0]->isParUpdatable())) { weight -= 0.5; }
    if (!(states[1]->isParUpdatable())) { weight -= 0.5; }

      // Compute the normals of the face
      for(CFuint id = 0; id < DIM_2D; ++id)
      normal[id] = -(boundaryNormals[iFace])[id];
      normal.normalize();
      const CFreal faceArea = currFace.computeVolume();

      //Compute Cp and cm on the face using quadrature rule
      CFuint nbqd = m_qdstates.size();
      for (CFuint iqd = 0 ; iqd < nbqd; ++iqd ){
         m_mappedCoord_line[0] = m_qd0[iqd]*(-1.0) + m_qd1[iqd]*(1.0);

        shapeFunctionValues = currFace.computeShapeFunctionAtMappedCoord(m_mappedCoord_line);
        (*m_qdstates[iqd]) = 0.0;

        for (CFuint iState = 0; iState < nbStatesInFace; ++ iState){
          (*m_qdstates[iqd]) += m_dimState_face[iState]*shapeFunctionValues[iState];
        }

        m_physicalCoord = currFace.computeCoordFromMappedCoord(m_mappedCoord_line);
        CFreal rho  = (*m_qdstates[iqd])[0];
        CFreal rhoU = (*m_qdstates[iqd])[1];
        CFreal rhoV = (*m_qdstates[iqd])[2];
        CFreal rhoE = (*m_qdstates[iqd])[3];
        CFreal invRho = 1./rho;
        CFreal rhoK2  = 0.5*(rhoU*rhoU + rhoV*rhoV)*invRho;
        CFreal p = gammaMinus1*(rhoE - rhoK2);

        // If we use the exact normals of the NACA-0012, we compute it
        if (m_exact_norm){

          // Coordinate of the vertex
          const CFreal x0 = m_physicalCoord[XX];
          const CFreal y0 = m_physicalCoord[YY];

          if ((y0 >= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
            m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 -
                                     0.1015*4.0*x0*x0*x0);
            m_qdNormal[YY] = 1.0;
            }
          else if ((y0 <= 0.0) && (x0 != 0.0) && (x0 < 1.0)){
            m_qdNormal[XX] = -0.6*(( 0.2969/(2.0*sqrt(x0))) - 0.1260 - 0.3516*2.0*x0 + 3.0*0.2843*x0*x0 -
                                     0.1015*4.0*x0*x0*x0);
            m_qdNormal[YY] = -1.0;
          }
          else if ((x0 == 0.0)){
            m_qdNormal[XX] = -1.0;
            m_qdNormal[YY] = 0.0;
          }
          else {m_qdNormal = normal;}
          m_qdNormal /= m_qdNormal.norm2();
          m_qdNormal *= normal.norm2();
          }

          else {m_qdNormal = normal; }

        node_px[iqd] =  (m_pInf- p)*m_qdNormal[XX];
        node_py[iqd] =  (m_pInf- p)*m_qdNormal[YY];
        node_m[iqd] = (m_pInf-p)*
                      ((m_physicalCoord[XX]-m_Xref[XX])*m_qdNormal[YY]-(m_physicalCoord[YY]-m_Xref[YY])*m_qdNormal[XX]);


      }
      faceCp = 0.0;
      m = 0.0;
      for (CFuint iqd = 0 ; iqd < nbqd; ++iqd ){
        faceCp[XX] += m_wqd[iqd]*node_px[iqd]*0.5;
        faceCp[YY] += m_wqd[iqd]*node_py[iqd]*0.5;
        m += m_wqd[iqd]*node_m[iqd]*0.5;
      }
    // compute Cp, scaled by wieght, depending on ownership of face to processor
    // press-m_pInf are switched to take into account the minus sign in normal
    Cp[XX] += weight *faceArea* faceCp[XX];
    Cp[YY] += weight *faceArea*faceCp[YY];

    Cmp -= m*weight*faceArea ;
    if(m_isViscous){
    //Compute Cf on the face using quadrature rule
     // Build the neighbour cell
    // handle to the neighbor cell
    DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
      faceNeighCell = socket_faceNeighCell.getDataHandle();

    const CFuint cellID = faceNeighCell[FaceID].first;
    const CFuint iFaceLocal = faceNeighCell[FaceID].second;
    geoDataCell.trs = faceNeighCell[FaceID].third;
    geoDataCell.idx = cellID;
    GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

    vector<Framework::State*> * m_states= neighborCell.getStates();
    vector<Framework::State*>& states_cell = (*m_states);
    const CFuint nbStatesInCell = states_cell.size();
    shapeFunctionValues.resize(nbStatesInCell);
    m_dimState_cell.resize(nbStatesInCell);
    for (CFuint s = 0; s < nbStatesInCell; ++s){
      m_dimState_cell[s].resize(nbEqs);
    }
    for (CFuint s = 0; s < nbStatesInCell; ++s){
      m_updateVarSet->setDimensionalValues(*states_cell[s],m_dimState_cell[s]);
    }

    // Depnding on which face we are mapped coordinates are not the same
    Node X0;
    Node X1;
    switch(iFaceLocal) {

	case 0:
	   X0[XX] = 0.0;
           X0[YY] = 0.0;

	   X1[XX] = 1.0;
           X1[YY] = 0.0;
         break;

	case 1:
	   X0[XX] = 1.0;
           X0[YY] = 0.0;

           X1[XX] = 0.0;
           X1[YY] = 1.0;
	break;
	case 2:
           X0[XX] = 0.0;
           X0[YY] = 1.0;

           X1[XX] = 0.0;
           X1[YY] = 0.0;
        break;

    }
    for (CFuint iqd = 0 ; iqd < nbqd; ++iqd ){

      //Here the mapped coordinates are different because we are on one face of
      // the reference element that is between 0 and 1 (and not -1;1 as the refrence line)
      m_mappedCoord[0] = .5*(1.0 - m_qd0[iqd])*X0[XX]+.5*(1.0 - m_qd1[iqd])*X1[XX];
      m_mappedCoord[1] = .5*(1.0 - m_qd0[iqd])*X0[YY]+.5*(1.0 - m_qd1[iqd])*X1[YY];
      // compute rho, u, v ate qd pts
      shapeFunctionValues = neighborCell.computeShapeFunctionAtMappedCoord(m_mappedCoord);
        (*m_qdstates[iqd]) = 0.0;
        for (CFuint iState = 0; iState < nbStatesInCell; ++ iState){
          (*m_qdstates[iqd]) += m_dimState_cell[iState]*shapeFunctionValues[iState];
        }
      // compute gradients at qd pts
      m_physicalCoord = neighborCell.computeCoordFromMappedCoord(m_mappedCoord);
      std::vector<RealVector> phyCoord;
      phyCoord.resize(1);
      phyCoord[0].resize(dim);
      phyCoord[0] = m_physicalCoord;

      const std::vector<RealMatrix> cellGradients= neighborCell.computeSolutionShapeFunctionGradients(phyCoord);
      for (CFuint id = 0; id < DIM_2D ; ++id){
      for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){
        (*m_gradients[iEq])[id] = 0.0;
      }
     }
      for (CFuint iEq = 0; iEq < nbEqs; ++ iEq){
        for (CFuint id = 0; id < dim; ++ id){
          for (CFuint iState = 0; iState < nbStatesInCell; ++iState){
              (*m_gradients[iEq])[id] +=   (cellGradients[0])(iState,id)*m_dimState_cell[iState][iEq];
          }
        }
      }
      qdValues[0] = (*m_qdstates[iqd])[0];
      qdValues[1] = (*m_qdstates[iqd])[1];
      qdValues[2] = (*m_qdstates[iqd])[2];
      qdValues[3] = (*m_qdstates[iqd])[3];
      CFreal mu = 0.0;
      if (m_diffVar.isNotNull()) mu = m_diffVar->getDynViscosity(qdValues, m_gradients);

      RealVector gradients_u;
      gradients_u.resize(dim);

      RealVector gradients_v;
      gradients_v.resize(dim);

      gradients_u[XX] = (1.0/qdValues[0])*((*m_gradients[1])[XX] - qdValues[1]*(*m_gradients[0])[XX]);
      gradients_u[YY] = (1.0/qdValues[0])*((*m_gradients[1])[YY] - qdValues[1]*(*m_gradients[0])[YY]);

      gradients_v[XX] = (1.0/qdValues[0])*((*m_gradients[2])[XX] - qdValues[2]*(*m_gradients[0])[XX]);
      gradients_v[YY] = (1.0/qdValues[0])*((*m_gradients[2])[YY] - qdValues[2]*(*m_gradients[0])[YY]);



      // Need to see if we are on the bottom part or on the top.
      // This is importante to see if the axes refered to the boundary is well oriented
      CFreal tau;
      if (normal[YY] > 0.0){
        tau = mu*(  normal[YY]*(  gradients_u[YY]*normal[YY] - gradients_v[YY]*normal[XX] )
              + normal[XX]*(  gradients_u[XX]*normal[YY] - gradients_v[XX]*normal[XX] ) );
      }

      else{
        tau = mu*(  normal[YY]*( -gradients_u[YY]*normal[YY] + gradients_v[YY]*normal[XX] )
              + normal[XX]*( -gradients_u[XX]*normal[YY] + gradients_v[XX]*normal[XX] ) );
      }

  /// @todo NV: This still need to be checked
    if (normal[YY] > 0.0){
        m_QdCf[iqd][XX] = tau*normal[YY] ;
        m_QdCf[iqd][YY] = -tau*normal[XX] ;

    }
    else{
        m_QdCf[iqd][XX] = -tau*normal[YY] ;
        m_QdCf[iqd][YY] = tau*normal[XX] ;

        }


      m_physicalCoord = neighborCell.computeCoordFromMappedCoord(m_mappedCoord);
      m_QdCmf[iqd] = (m_physicalCoord[YY]-m_Xref[YY])*m_QdCf[iqd][XX] - (m_physicalCoord[XX]-m_Xref[XX])*m_QdCf[iqd][YY];

    }
    geoBuilderCell.releaseGE();
    faceCf = 0.0;
    faceCmf = 0.0;
    for (CFuint iqd = 0 ; iqd < nbqd; ++iqd ){
      faceCf += m_QdCf[iqd]*m_wqd[iqd]*0.5  ;
      faceCmf += m_QdCmf[iqd]*m_wqd[iqd]*0.5 ;

    }

    faceCf *= weight *faceArea;
    faceCmf *= weight *faceArea;

    Cf += faceCf;
    Cmf += faceCmf;

    }
    geoBuilderFace.releaseGE();
  }
  // Compute the value of m_alpha
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(&m_eval[0]);
  // Transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;
  // adimensionalize Cp, Cm and Cf
  Cp /= (0.5 * m_rhoInf * m_uInf * m_uInf);
  Cmp /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength);
  if (m_isViscous){
    Cf /=  (0.5 * m_rhoInf * m_uInf * m_uInf);
    Cmf /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength) ;
  }
  
  // project Cp
  actual_lift = (sin(m_alpharad)*(Cp[XX] + Cf[XX]) + cos(m_alpharad)*(Cp[YY] + Cf[YY])) / refLength;
  actual_drag_p = (cos(m_alpharad)*Cp[XX] - sin(m_alpharad)*Cp[YY]) / refLength;
  actual_drag_f = 0.0;
  actual_momentum_f = 0.0;
  if (m_isViscous){
    actual_drag_f = (cos(m_alpharad)*Cf[XX] - sin(m_alpharad)*Cf[YY]) / refLength;
    actual_momentum_f = Cmf ;
  }
  actual_drag = actual_drag_p + actual_drag_f;
  actual_momentum_p = Cmp / refLength;
  actual_momentum = Cmp + Cmf;


  /// @TODO NV: Has been copy from RMSJouleHeatSource.cxx where this solution was said to be temporary....
  total_lift = 0.0;
  total_drag_p = 0.0;
  total_drag_f = 0.0;
  total_drag = 0.0;
  total_momentum_p = 0.0;
  total_momentum_f = 0.0;
  total_momentum = 0.0;

  const std::string nsp = getMethodData().getNamespace();
  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
  MPI_Allreduce(&actual_lift, &total_lift, 1, MPI_DOUBLE, MPI_SUM,
        PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_drag_p, &total_drag_p, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_drag_f, &total_drag_f, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_drag, &total_drag, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum_p, &total_momentum_p, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum_f, &total_momentum_f, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum, &total_momentum, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator(nsp));
#endif
  }
  else {

    total_lift = actual_lift;
    total_drag_p = actual_drag_p;
    total_drag_f = actual_drag_f;
    total_drag = actual_drag;
    total_momentum_p = actual_momentum_p;
    total_momentum_f = actual_momentum_f;
    total_momentum = actual_momentum;

  }
  // Output to file the coefficients
  updateOutputFileAero();

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::unsetup()
{

  CFuint nbqd = m_qdstates.size();
  for (CFuint i = 0; i < nbqd; ++i) {
   deletePtr(m_qdstates[i]);
   deletePtr(m_qdnodes[i]);
  }


}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::prepareOutputFileWall()
{
  using boost::filesystem::path;

  cf_assert (!m_fileWall->isopen());

  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileWall);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );

  ofstream& fout = m_fileWall->open(file);

    if (m_isViscous){
      fout << "TITLE  =  Values at the Wall" << "\n";
//      fout << "VARIABLES = x y Cp Cf Mach Pressure Temperature Density" << "\n";
      fout << "VARIABLES = x y z Cp Cf Pressure Temperature Density" << "\n";
    }
    else{
    fout << "TITLE  =  Non viscous Values at the Wall" << "\n";
//    fout << "VARIABLES = x y Cp Mach Pressure Temperature Density" << "\n";
    fout << "VARIABLES = x y z Cp Pressure Temperature Density" << "\n";
    }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::prepareOutputFileAero()
{
  boost::filesystem::path fpath = m_nameOutputFileAero;
  fpath = PathAppender::getInstance().appendAllInfo(fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  if (m_isViscous){
    convergenceFile <<"TITLE  =  Aerodynamic Coefficients"  << "\n";
    convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef_p DragCoef_f DragCoef_tot MomentumCoef_p MomentumCoef_f Momentum_tot" << "\n";
  }
  else{
    convergenceFile <<"TITLE  =  Non viscous Aerodynamic Coefficients"  << "\n";
    convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef MomentumCoef" << "\n";
 }

  convergenceFile.close();

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHO::updateOutputFileAero()
{
  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendAllInfo( fpath, false, false);

  SelfRegistPtr<FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  if (m_isViscous){
    convergenceFile
    << subSysStatus->getNbIter()       << " "
    << subSysStatus->getCurrentTime()  << " "
    << m_alphadeg                      << " "
    << total_lift                      << " "
    << total_drag_p                    << " "
    << total_drag_f                    << " "
    << total_drag                      << " "
    << total_momentum_p                << " "
    << total_momentum_f                << " "
    << total_momentum                  << "\n";
  }

  else{
    convergenceFile
      << subSysStatus->getNbIter()       << " "
      << subSysStatus->getCurrentTime()  << " "
      << m_alphadeg                      << " "
      << total_lift                      << " "
      << total_drag                      << " "
      << total_momentum                  << "\n";
  }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

