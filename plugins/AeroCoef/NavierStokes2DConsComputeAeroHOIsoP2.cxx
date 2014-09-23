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
#include "AeroCoef/NavierStokes2DConsComputeAeroHOIsoP2.hh"

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

MethodCommandProvider<NavierStokes2DConsComputeAeroHOIsoP2,
                      DataProcessingData,
                      AeroCoefFSModule>
aNavierStokes2DConsComputeAeroHOIsoP2Provider("NavierStokes2DConsComputeAeroHOIsoP2");

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::defineConfigOptions(Config::OptionList& options)
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

NavierStokes2DConsComputeAeroHOIsoP2::NavierStokes2DConsComputeAeroHOIsoP2(const std::string& name) :
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
  m_dataState(),
  m_subface_center(2),
  m_weights(nbQdPts),
  m_qpPos(nbQdPts)
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

NavierStokes2DConsComputeAeroHOIsoP2::~NavierStokes2DConsComputeAeroHOIsoP2()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokes2DConsComputeAeroHOIsoP2::needsSockets()
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
NavierStokes2DConsComputeAeroHOIsoP2::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::setup()
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

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_qdNormal.resize(dim);

  //Centers of subfaces in reference domain <0,1>
  m_subface_center[0] = 0.25;
  m_subface_center[1] = 0.75;

  //Weights at Gauss points
  m_weights[0] = 5./36.;
  m_weights[1] = 8./36.;
  m_weights[2] = 5./36.;

  const CFreal s = std::sqrt(0.6);

  //Position of quadrature points
  m_qpPos[0] = -0.25*s;
  m_qpPos[1] =  0.0;
  m_qpPos[2] =  0.25*s;

  m_qdState = new State();

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::configure ( Config::ConfigArgs& args )
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

void NavierStokes2DConsComputeAeroHOIsoP2::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokes2DConsComputeAeroHOIsoP2IsoP2::setFuntion(): alpha function wrongly defined.");

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

void NavierStokes2DConsComputeAeroHOIsoP2::executeOnTrs()
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

void NavierStokes2DConsComputeAeroHOIsoP2::computeWall()
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

    fout
    << (state.getCoordinates())[XX]   << " "
    << (state.getCoordinates())[YY]   << " "
    << Cp     << " "
    << Mach   << " "
    << P      << " "
    << T      << " "
    << rho    << "\n";


  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::computeAero()
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



  RealVector normal(0.0, dim);
  RealVector Cp(0.0, dim);
  RealVector faceCp(0.0, dim);
  RealVector Cf(0.0, dim);
  CFreal Cmp = 0.0;
  CFreal Cmf = 0.0;
  CFreal m;

  std::vector<RealVector> m_coords(1);
  m_coords[0].resize(2);

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

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    geoDataFace.idx = iFace;

     GeometricEntity & currFace = *geoBuilderFace.buildGE();

    std::vector<State*>& states = *(currFace.getStates());
    std::vector<Node*>& nodes = *(currFace.getNodes());
    cf_assert(states.size() >= 2);

    // Check if face is owned by one or more processors and compute the weight accordingly
    CFreal weight = 1.0;
    if (!(states[0]->isParUpdatable())) { weight -= 0.5; }
    if (!(states[1]->isParUpdatable())) { weight -= 0.5; }
    // Compute the normals and the jacobian in each node of the face
    // node0

  //Initialize pressure coefficient
  faceCp = 0.0;
  //Initialize pitching moment
  m = 0.0;


  for(CFuint ibdrysubface = 0; ibdrysubface < 2; ++ibdrysubface) {


    for(CFuint iQd = 0; iQd < nbQdPts; ++iQd) {
      const CFreal xiQd = m_qpPos[iQd] + m_subface_center[ibdrysubface];

      const CFreal p0 = (1.0-3.0*xiQd+2.0*xiQd*xiQd);
      const CFreal p1 = (2.0*xiQd*xiQd-xiQd);
      const CFreal p2 = 4.0*xiQd*(1.0-xiQd);

      const CFreal xqd = p0*(*(nodes[0]))[XX] + p1*(*(nodes[1]))[XX] + p2*(*(nodes[2]))[XX];
      const CFreal yqd = p0*(*(nodes[0]))[YY] + p1*(*(nodes[1]))[YY] + p2*(*(nodes[2]))[YY];
      (*m_qdState) = p0*(*states[0]) + p1*(*states[1]) + p2*(*states[2]);

      m_updateVarSet->setDimensionalValues(*m_qdState,m_dimState);

      const CFreal rho  = m_dimState[0];
      const CFreal rhoU = m_dimState[1];
      const CFreal rhoV = m_dimState[2];
      const CFreal rhoE = m_dimState[3];

      const CFreal invRho = 1./rho;
      const CFreal rhoK2  = 0.5*(rhoU*rhoU + rhoV*rhoV)*invRho;
      const CFreal p = gammaMinus1*(rhoE - rhoK2);

      m_CP2N.ComputeBNormal(nodes,xiQd,m_qdNormal);

      m_qdNormal = -1.0 * m_qdNormal;
      const CFreal jacob = m_qdNormal.norm2();
      const CFreal invJacob = 1.0/jacob;
      m_qdNormal = invJacob * m_qdNormal;

//       m += jacob * m_weights[iQd] * (m_pInf - p)*((xqd-m_Xref[XX])*m_qdNormal[YY]-(yqd-m_Xref[YY])*m_qdNormal[XX]);
      m += jacob * m_weights[iQd] * (m_pInf - p)*((xqd-m_Xref[XX])*m_qdNormal[YY]-(yqd-m_Xref[YY])*m_qdNormal[XX]);
      faceCp[XX] += jacob * m_weights[iQd] * (m_pInf- p)*m_qdNormal[XX];
      faceCp[YY] += jacob * m_weights[iQd] * (m_pInf- p)*m_qdNormal[YY];

    }

  } ///Loop over subfaces of the big face


    // compute Cp, scaled by weight, depending on ownership of face to processor
    // press-m_pInf are switched to take into account the minus sign in normal
    Cp[XX] += weight *faceCp[XX];
    Cp[YY] += weight *faceCp[YY];

    Cmp -= m*weight ;

geoBuilderFace.releaseGE();
}  //loop over faces


  // Compute the value of m_alpha
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(m_eval);
  // Transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;
  // adimensionalize Cp
  Cp /= (0.5 * m_rhoInf * m_uInf * m_uInf);
  Cmp /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength);
    // project Cp
  actual_lift = (sin(m_alpharad)*(Cp[XX] + Cf[XX]) + cos(m_alpharad)*(Cp[YY] + Cf[YY])) / refLength;
  actual_drag_p = (cos(m_alpharad)*Cp[XX] - sin(m_alpharad)*Cp[YY]) / refLength;

    actual_drag_f = 0.0;
    actual_momentum_f = 0.0;
  actual_drag = actual_drag_p + actual_drag_f;
  actual_momentum_p = Cmp / refLength;
  actual_momentum_f = Cmf / refLength;
  actual_momentum = Cmp + Cmf;


  /// @TODO NV: Has been copy from RMSJouleHeatSource.cxx where this solution was said to be temporary....
  total_lift = 0.0;
  total_drag_p = 0.0;
  total_drag_f = 0.0;
  total_drag = 0.0;
  total_momentum_p = 0.0;
  total_momentum_f = 0.0;
  total_momentum = 0.0;

  if (PE::GetPE().GetProcessorCount() > 1) {

#ifdef CF_HAVE_MPI
  MPI_Allreduce(&actual_lift, &total_lift, 1, MPI_DOUBLE, MPI_SUM,
        PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_drag_p, &total_drag_p, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_drag_f, &total_drag_f, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_drag, &total_drag, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_momentum_p, &total_momentum_p, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_momentum_f, &total_momentum_f, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
  MPI_Allreduce(&actual_momentum, &total_momentum, 1, MPI_DOUBLE, MPI_SUM,
      PE::GetPE().GetCommunicator());
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

void NavierStokes2DConsComputeAeroHOIsoP2::unsetup()
{
  delete m_qdState;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::prepareOutputFileWall()
{
  using boost::filesystem::path;

  cf_assert (!m_fileWall->isopen());

  path file = Environment::DirPaths::getInstance().getResultsDir() / path (m_nameOutputFileWall);
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );

  ofstream& fout = m_fileWall->open(file);

    fout << "TITLE  =  Non viscous Values at the Wall" << "\n";
    fout << "VARIABLES = x y Cp Mach Pressure Temperature Density" << "\n";

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::prepareOutputFileAero()
{
  boost::filesystem::path fpath = m_nameOutputFileAero;
  fpath = PathAppender::getInstance().appendAllInfo(fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

    convergenceFile <<"TITLE  =  Non viscous Aerodynamic Coefficients"  << "\n";
    convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef MomentumCoef" << "\n";


  convergenceFile.close();

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes2DConsComputeAeroHOIsoP2::updateOutputFileAero()
{
  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendAllInfo( fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

    convergenceFile
      << subSysStatus->getNbIter()       << " "
      << subSysStatus->getCurrentTime()  << " "
      << m_alphadeg                      << " "
      << total_lift                      << " "
      << total_drag                      << " "
      << total_momentum                  << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

