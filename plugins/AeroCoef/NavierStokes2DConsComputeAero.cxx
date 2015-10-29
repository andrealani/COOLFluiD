#include "Common/PE.hh"

#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Common/VarRegistry.hh"
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
#include "AeroCoef/NavierStokes2DConsComputeAero.hh"

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

MethodCommandProvider<NavierStokes2DConsComputeAero,
                      DataProcessingData,
                      AeroCoefFSModule>
aNavierStokes2DConsComputeAeroProvider("NavierStokes2DConsComputeAero");

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::defineConfigOptions(Config::OptionList& options)
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

}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DConsComputeAero::NavierStokes2DConsComputeAero( const std::string& name) :
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
  m_dimState(),
  _cells(CFNULL),
  m_diffVar(),
  _gradients(),
  _avValues(),
  m_fsData(CFNULL),
  m_nodal_cf(),
  m_numb_node(),
  m_dataState()
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

}

//////////////////////////////////////////////////////////////////////////////

NavierStokes2DConsComputeAero::~NavierStokes2DConsComputeAero()
{ for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokes2DConsComputeAero::needsSockets()
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
NavierStokes2DConsComputeAero::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::setup()
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
  _gradients.resize(PhysicalModelStack::getActive()->getNbEq());
  for (CFuint i = 0; i< nbEqs; ++i)
  {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  _avValues.resize(PhysicalModelStack::getActive()->getNbEq());
  prepareOutputFileAero();

  // put the values in the subsystem var registry
  Common::SafePtr<VarRegistry> ssys_var_regist =
    SubSystemStatusStack::getActive()->getVarRegistry();

  ssys_var_regist->registVar<CFreal>("CL", new CFreal(0.0));
  ssys_var_regist->registVar<CFreal>("CD", new CFreal(0.0));
  ssys_var_regist->registVar<CFreal>("CM", new CFreal(0.0));
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::configure ( Config::ConfigArgs& args )
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

  // resize vector of dimensional values
  m_dimState.resize(4);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokes2DConsComputeAero::setFuntion(): alpha function wrongly defined.");

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

void NavierStokes2DConsComputeAero::executeOnTrs()
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

void NavierStokes2DConsComputeAero::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file is opened here

  const CFreal R     = m_updateVarSet->getModel()->getRdim();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();


  //If computeAero has not been done before,we have to compte the Cf for each node of the TRS
  if(m_saveRateAero != m_saveRateWall){
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

    m_nodal_cf.resize(statesIdx->size(),2);
    RealVector normal(0.0, DIM_2D);

    CFuint n = 0;

 if (m_isViscous){
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
      u.resize(nbStatesInCell);
      v.resize(nbStatesInCell);

      //compute an average of the coords in the face
      for (CFuint s = 0; s < nbNodesInFace; ++s){
        m_updateVarSet->setDimensionalValues(*states[s],m_dimState);

        m_coords[0] += states[s]->getCoordinates();
      }

      m_coords[0]/= nbNodesInFace;

      // Compute the shape function gradient at the projected point
      const std::vector<RealMatrix> cellGradients= neighborCell.computeSolutionShapeFunctionGradients(m_coords);

      const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

      ///@todo NV: this conversion from rhou, rhov to u,v is not very good....should be done in a better way

      for (CFuint iState = 0; iState < nbStatesInCell ; ++iState){
        u[iState] = (*states_cell[iState])[1]/(*states_cell[iState])[0];
        v[iState] = (*states_cell[iState])[2]/(*states_cell[iState])[0];
      }

      // Compute the gradients variables
      RealVector gradients_u(0., 2);
      RealVector gradients_v(0., 2);

      for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){

        for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){

          (*_gradients[iEq])[iDim] = 0.0;
        }
      }

      for (CFuint iNode = 0; iNode < nbStatesInCell ; ++iNode){
        for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){

          (gradients_u)[iDim] += (cellGradients[0])(iNode,iDim) * u[iNode] ;
          (gradients_v)[iDim] += (cellGradients[0])(iNode,iDim) * v[iNode] ;

          for (CFuint iEq = 0; iEq < nbEqs ; ++iEq)
            (*_gradients[iEq])[iDim] += (cellGradients[0])(iNode,iDim) *(*states_cell[iNode])[iEq];
          }

      }

      CFreal mu = 0.;
      if (m_diffVar.isNotNull()) mu = m_diffVar->getDynViscosity(_avValues, _gradients);

      // Need to see if we are on the bottom part or on the top.
      // This is importante to see if the axes refered to the boundary is well oriented

      if (normal[YY] > 0.0){
        tau = mu*( normal[YY]*( gradients_u[YY]*normal[YY] - gradients_v[YY]*normal[XX] )
                      + normal[XX]*(gradients_u[XX]*normal[YY] - gradients_v[XX]*normal[XX] ) );
      }

      else {
       tau = mu*( normal[YY]*(-gradients_u[YY]*normal[YY] + gradients_v[YY]*normal[XX] )
                           + normal[XX]*(-gradients_u[XX]*normal[YY] + gradients_v[XX]*normal[XX] ) );
      }


      // m_nodal_cf contain the cf of each node of the TRS, we do an average between the cf of faces
      // and we store the localID of the node in a vector. Then it is necessary to know if m_nodal_cf
      // has already a part of the value of the node

      // looking at the first node of the face
      if (m_numb_node[n] == states[0]->getLocalID()){
        m_nodal_cf(n,0) += tau*faceArea;
        m_nodal_cf(n,1) += faceArea;
      }
      else {
        CFuint marker = 0;
        for (CFuint j = 0; j< n; ++j){
          if(m_numb_node[j] == states[0]->getLocalID()){
            m_nodal_cf(j,0) += tau*faceArea;
            m_nodal_cf(j,1) += faceArea;
            marker = 1;
          }
        }
        if (marker == 0){
          m_nodal_cf(n,0) += tau*faceArea;
          m_nodal_cf(n,1) += faceArea;
          m_numb_node[n] =  states[0]->getLocalID();
          n ++;
        }
      }

      //looking at the second node of the face
      CFuint marker=0;
      for (CFuint j = 0; j< n; ++j){
        if(m_numb_node[j] == states[1]->getLocalID()){
          m_nodal_cf(j,0) += tau*faceArea;
          m_nodal_cf(j,1) += faceArea;
          marker = 1;
        }
      }
      if (marker == 0){
        m_nodal_cf(n,0) += tau*faceArea;
        m_nodal_cf(n,1) += faceArea;
        m_numb_node[n] =  states[1]->getLocalID();
        n+= 1;
      }

      //the number that are at the interior of the mesh are met for the first time
      for (CFuint s = 2; s < nbNodesInFace; ++s){
        m_nodal_cf(n,0) += tau*faceArea;
        m_nodal_cf(n,1) += faceArea;
        m_numb_node[n] =  states[s]->getLocalID();
        n+= 1;
      }

      geoBuilderFace.releaseGE();
      geoBuilderCell.releaseGE();
    }

    //adimensionalize nodal_cf and devid by the coef for average

    for ( CFuint i = 0; i < statesIdx->size(); ++i)
        m_nodal_cf(i,0) /= m_nodal_cf(i,1)*(0.5 * m_rhoInf * m_uInf * m_uInf);

  }

  }

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState){
    const State& state = *(states[(*statesIdx)[iState]]);
    m_updateVarSet->setDimensionalValues(state,m_dimState);

    const CFreal rho  = m_dimState[0];
    const CFreal rhoU = m_dimState[1];
    const CFreal rhoV = m_dimState[2];
    const CFreal rhoE = m_dimState[3];

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

    // calculate mach and speed of sound
    const CFreal a2 = gamma*R*T;
    const CFreal Mach = sqrt( (U2 + V2) / a2 );

// unused // const CFreal y = state.getCoordinates()[YY];

    CFreal Cf = 0.;

    //Find the cf corresponding to this sate
    if (m_isViscous){
      for (CFuint i = 0; i < statesIdx->size(); ++ i){
        if (m_numb_node[i] == state.getLocalID())
          Cf = m_nodal_cf(i,0);
      }
  // output to File
    fout
    << (state.getCoordinates())[XX]   << " "
    << (state.getCoordinates())[YY]   << " "
    << Cp     << " "
    << Cf     << " "
    << Mach   << " "
    << P      << " "
    << T      << " "
    << rho    << "\n";
    }
  else{ // output to File
    fout
    << (state.getCoordinates())[XX]   << " "
    << (state.getCoordinates())[YY]   << " "
    << Cp     << " "
    << Mach   << " "
    << P      << " "
    << T      << " "
    << rho    << "\n";
}

  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::computeAero()
{
  CFAUTOTRACE;

  // get the datahandle containing the boundary Normals
  const std::string TRSName = getCurrentTRS()->getName();
  const std::string socketName = TRSName + "-boundaryNormals";
  DataHandle<const CFreal*> boundaryNormals = m_sockets.getSocketSource<const CFreal*>(socketName)->getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_assert ( dim == DIM_2D );

  const CFreal gammaMinus1 = m_updateVarSet->getModel()->getGamma() - 1.;
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  CFreal tau = 0.0;
  vector < CFreal > u;
  vector < CFreal > v;


  RealVector normal(0.0, dim);
  RealVector Cp(0.0, dim);
  RealVector Cf(0.0, dim);
  CFreal Cmp = 0.0;
  CFreal Cmf = 0.0;
  CFreal m;
  CFreal p;

  std::vector<RealVector> m_coords(1);
  m_coords[0].resize(2);

  if (m_isViscous){
  // computing number of states on the TRS. It is used to size the vector with friction coefficient at wall
  const CFuint nb_statesTRS = (getCurrentTRS()->getStatesInTrs())->size();
  m_nodal_cf.resize(nb_statesTRS,2);
  m_numb_node.resize(nb_statesTRS);
  for (CFuint i = 0; i < nb_statesTRS; ++i)
    m_numb_node[i] = std::numeric_limits<CFuint>::max();

  }
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

  CFuint n = 0;

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    geoDataFace.idx = iFace;
    GeometricEntity & currFace = *geoBuilderFace.buildGE();

    // get the face normal
    const CFuint FaceID = currFace.getID();
    for(CFuint iDim = 0; iDim < DIM_2D; ++iDim)
      normal[iDim] = -(boundaryNormals[iFace])[iDim];
    normal.normalize();

    std::vector<State*>& states = *(currFace.getStates());
    cf_assert(states.size() >= 2);
    const CFuint  nbNodesInFace = states.size();

    // Check if face is owned by one or more processors and compute the weight accordingly
    CFreal weight = 1.0;
    if (!(states[0]->isParUpdatable())) { weight -= 0.5; }
    if (!(states[1]->isParUpdatable())) { weight -= 0.5; }

    // Compute the pressure on a face
    CFreal press = 0.;
    _avValues = 0.0;
    m = 0.;
    m_coords[0] = 0.0;

    /// @todo NV: this conversion from rhou, rhov to u,v
    ///           could be done usnign variable transformations
    for (CFuint s = 0; s < nbNodesInFace; ++s){
      m_updateVarSet->setDimensionalValues(*states[s],m_dimState);
      const RealVector x = (*states[s]).getCoordinates();

      const CFreal rho  = m_dimState[0];
      const CFreal rhoU = m_dimState[1];
      const CFreal rhoV = m_dimState[2];
      const CFreal rhoE = m_dimState[3];

      const CFreal invRho = 1./rho;
      const CFreal rhoK2  = 0.5*(rhoU*rhoU + rhoV*rhoV)*invRho;
      press += gammaMinus1*(rhoE - rhoK2);
      p =  gammaMinus1*(rhoE - rhoK2);
      m += (m_pInf - p)*((x[XX]-m_Xref[XX])*normal[YY]-(x[YY]-m_Xref[YY])*normal[XX]);
      _avValues[0] +=  rho;
      _avValues[1] +=  rhoU;
      _avValues[2] +=  rhoV;
      _avValues[3] +=  rhoE;
      m_coords[0] += states[s]->getCoordinates();
    }

    press /= nbNodesInFace;
    m /= nbNodesInFace;
    _avValues /= nbNodesInFace;
    m_coords[0]/= nbNodesInFace;

    const CFreal faceArea = currFace.computeVolume();

    // compute Cp, scaled by wieght, depending on ownership of face to processor
    // press-m_pInf are switched to take into account the minus sign in normal
    Cp[XX] += weight * (m_pInf-press) * faceArea * normal[XX];
    Cp[YY] += weight * (m_pInf-press) * faceArea * normal[YY];

    Cmp -= m*weight*faceArea ;


    if(m_isViscous){
    // Build the neighbour cell
    // handle to the neighbor cell
    DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
      faceNeighCell = socket_faceNeighCell.getDataHandle();

    const CFuint cellID = faceNeighCell[FaceID].first;
    geoDataCell.trs = faceNeighCell[FaceID].third;
    geoDataCell.idx = cellID;
    GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

    vector<Framework::State*> * m_states= neighborCell.getStates();
    vector<Framework::State*>& states_cell = (*m_states);
    const CFuint nbStatesInCell = m_states->size();
    u.resize(nbStatesInCell);
    v.resize(nbStatesInCell);


    // Compute the shape function gradient at the projected point
    const std::vector<RealMatrix> cellGradients= neighborCell.computeSolutionShapeFunctionGradients(m_coords);

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    for (CFuint iState = 0; iState < nbStatesInCell ; ++iState){
      u[iState] = (*states_cell[iState])[1]/(*states_cell[iState])[0];
      v[iState] = (*states_cell[iState])[2]/(*states_cell[iState])[0];
    }


    // Compute the gradients variables
    RealVector gradients_u(0., 2);
    RealVector gradients_v(0., 2);

    for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){
      for (CFuint iEq = 0; iEq < nbEqs ; ++iEq){
        (*_gradients[iEq])[iDim] = 0.0;
      }
    }

    for (CFuint iNode = 0; iNode < nbStatesInCell ; ++iNode){
      for (CFuint iDim = 0; iDim < DIM_2D ; ++iDim){
        (gradients_u)[iDim] += (cellGradients[0])(iNode,iDim) * u[iNode] ;
        (gradients_v)[iDim] += (cellGradients[0])(iNode,iDim) * v[iNode] ;
          for (CFuint iEq = 0; iEq < nbEqs ; ++iEq)
            (*_gradients[iEq])[iDim] += (cellGradients[0])(iNode,iDim) *(*states_cell[iNode])[iEq];
    }

  }

  CFreal mu = 0.;
  if (m_diffVar.isNotNull()) mu = m_diffVar->getDynViscosity(_avValues, _gradients);

    // Need to see if we are on the bottom part or on the top.
    // This is importante to see if the axes refered to the boundary is well oriented

  if (normal[YY] > 0.0){
    tau = mu*(  normal[YY]*(  gradients_u[YY]*normal[YY] - gradients_v[YY]*normal[XX] )
              + normal[XX]*(  gradients_u[XX]*normal[YY] - gradients_v[XX]*normal[XX] ) );
  }

  else{
    tau = mu*(  normal[YY]*( -gradients_u[YY]*normal[YY] + gradients_v[YY]*normal[XX] )
              + normal[XX]*( -gradients_u[XX]*normal[YY] + gradients_v[XX]*normal[XX] ) );
  }

  CFreal coef  = tau*weight*faceArea  ;

  /// @todo NV: This still need to be checked
  if (normal[YY] > 0.0){
    Cf[XX] += coef*normal[YY] ;
    Cf[YY] -= coef*normal[XX] ;

    Cmf += (m_coords[0][YY]-m_Xref[YY])*(coef*normal[YY]) - (m_coords[0][XX]-m_Xref[XX])*(-coef*normal[XX]) ;
  }
  else{
    Cf[XX] -= coef*normal[YY] ;
    Cf[YY] += coef*normal[XX] ;

    Cmf += (m_coords[0][YY]-m_Xref[YY])*(-coef*normal[YY]) - (m_coords[0][XX]-m_Xref[XX])*(coef*normal[XX]) ;
  }

  //m_nodal_cf contain the cf of each node of the TRS, we do an average between the cf of faces
  // and we store the localID of the node in a vector. Then it is necessary to know if m_nodal_cf
  //has already a part of the value of the node

  //looking at the first (left) node of the face
  if (m_numb_node[n] == states[0]->getLocalID()){
    m_nodal_cf(n,0) += tau*faceArea;
    m_nodal_cf(n,1) += faceArea;
  }
  else {
    CFuint marker = 0;
    for (CFuint j = 0; j< n; ++j){
      if(m_numb_node[j] == states[0]->getLocalID()){
        m_nodal_cf(j,0) += tau*faceArea;
        m_nodal_cf(j,1) += faceArea;
        marker = 1;
      }
    }
    if (marker == 0){
      m_nodal_cf(n,0) += tau*faceArea;
      m_nodal_cf(n,1) += faceArea;
      m_numb_node[n] =  states[0]->getLocalID();
      n ++;
    }
  }

  // looking at the second (right) node of the face
  CFuint marker=0;
  for (CFuint j = 0; j< n; ++j){
    if(m_numb_node[j] == states[1]->getLocalID()){
      m_nodal_cf(j,0) += tau*faceArea;
      m_nodal_cf(j,1) += faceArea;
      marker = 1;
    }
  }
  if (marker == 0){
    m_nodal_cf(n,0) += tau*faceArea;
    m_nodal_cf(n,1) += faceArea;
    m_numb_node[n] =  states[1]->getLocalID();
    n+= 1;
  }


  //the number that are at the interior of the mesh can meet for the first time
  for (CFuint s = 2; s < nbNodesInFace; ++s){
    m_nodal_cf(n,0) += tau*faceArea;
    m_nodal_cf(n,1) += faceArea;
    m_numb_node[n] =  states[s]->getLocalID();
    n+= 1;
  }

  geoBuilderFace.releaseGE();
  geoBuilderCell.releaseGE();
}
}

  // Compute the value of m_alpha
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(&m_eval[0]);

  // Transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;

  // adimensionalize Cp
  Cp /= (0.5 * m_rhoInf * m_uInf * m_uInf);

  Cmp /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength);

  if (m_isViscous){
    Cf /=  (0.5 * m_rhoInf * m_uInf * m_uInf);
    Cmf /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength) ;
  }
  // project Cp
  actual_lift = (sin(m_alpharad)*(Cp[XX] + Cf[XX]) + cos(m_alpharad)*(Cp[YY] + Cf[YY])) / refLength;
  actual_drag_p = (cos(m_alpharad)*Cp[XX] - sin(m_alpharad)*Cp[YY]) / refLength;
  if (m_isViscous){
    actual_drag_f = (cos(m_alpharad)*Cf[XX] - sin(m_alpharad)*Cf[YY]) / refLength;
    actual_momentum_f = Cmf ;
  }
  else {
    actual_drag_f = 0.0;
    actual_momentum_f = 0.0;
  }
  actual_drag = actual_drag_p + actual_drag_f;
  actual_momentum_p = Cmp / refLength;
  actual_momentum_f = Cmf / refLength;
  actual_momentum = Cmp + Cmf;

  if (m_isViscous){
  //adimensionalize nodal_cf and devid by the coef for average

  const CFuint nb_statesTRS = m_nodal_cf.nbRows();
  for ( CFuint i = 0; i < nb_statesTRS; ++i)
    m_nodal_cf(i,0) /= m_nodal_cf(i,1)*(0.5 * m_rhoInf * m_uInf * m_uInf);

}

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

  // put the values in the subsystem var registry
  Common::SafePtr<SubSystemStatus> ssys_status =
    SubSystemStatusStack::getActive();
  Common::SafePtr<VarRegistry> ssys_var_regist =
    ssys_status->getVarRegistry();

  ssys_var_regist->setVar<CFreal>("CL", total_lift);
  ssys_var_regist->setVar<CFreal>("CD", total_drag);
  ssys_var_regist->setVar<CFreal>("CM", total_momentum);

  // Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::unsetup()
{
  // put the values in the subsystem var registry
  Common::SafePtr<VarRegistry> ssys_var_regist =
    SubSystemStatusStack::getActive()->getVarRegistry();

  CFreal * ptr_cl = ssys_var_regist->unregistVar<CFreal>("CL");
  CFreal * ptr_cd = ssys_var_regist->unregistVar<CFreal>("CD");
  CFreal * ptr_cm = ssys_var_regist->unregistVar<CFreal>("CM");

  deletePtr(ptr_cl);
  deletePtr(ptr_cd);
  deletePtr(ptr_cm);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::prepareOutputFileWall()
{
  using boost::filesystem::path;

  cf_assert (!m_fileWall->isopen());

  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileWall);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );

  ofstream& fout = m_fileWall->open(file);

  if (m_isViscous){
    fout << "TITLE  =  Values at the Wall" << "\n";
    fout << "VARIABLES = x y Cp Cf Mach Pressure Temperature Density" << "\n";
  }
  else{
    fout << "TITLE  =  Non viscous Values at the Wall" << "\n";
    fout << "VARIABLES = x y Cp Mach Pressure Temperature Density" << "\n";
  }

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAero::prepareOutputFileAero()
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

void NavierStokes2DConsComputeAero::updateOutputFileAero()
{
  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendAllInfo( fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
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

