#include "Common/PE.hh"
#include "Common/SafePtr.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "AeroCoef/AeroCoef.hh"
#include "AeroCoef/Euler2DConsComputeAeroFVMCC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Euler2DConsComputeAeroFVMCC, DataProcessingData, AeroCoefModule> Euler2DConsComputeAeroFVMCCProvider("Euler2DConsComputeAeroFVMCC");

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Definition of the function defining the angle between the flow and the airfoil.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coeficients.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFreal >("uInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("rhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("pInf","Pressure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsComputeAeroFVMCC::Euler2DConsComputeAeroFVMCC(const std::string& name) :
  DataProcessingCom(name),
  m_faceTrsGeoBuilder(),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_gstates("gstates"),
  m_updateVarSet(),
  m_vars("t"),
  m_eval(0.0,1),
  m_lift(),
  m_drag(),
  m_alphadeg(0.),
  m_alpharad(0.)
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
  setParameter("uInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("rhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("pInf",&m_pInf);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

Euler2DConsComputeAeroFVMCC::~Euler2DConsComputeAeroFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Euler2DConsComputeAeroFVMCC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Euler2DConsComputeAeroFVMCC::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::setup()
{
  m_lift = 0.;
  m_drag = 0.;

  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<Physics::NavierStokes::EulerVarSet>();
  m_updateVarSet->getModel()->resizePhysicalData(m_dataState);

  cf_assert(m_uInf > 0.);
  cf_assert(m_rhoInf > 0.);
  cf_assert(m_pInf > 0.);

  prepareOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::configure ( Config::ConfigArgs& args )
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

void Euler2DConsComputeAeroFVMCC::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"Euler2DConsComputeAeroFVMCC::setFuntion(): alpha function wrongly defined.");

  m_functionParser.Parse(m_function, m_vars);

  if (m_functionParser.ErrorMsg() != 0) {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(m_functionParser.ErrorMsg());
    msg += " Function: " + m_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute the value of the angle
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(&m_eval[0]);

  // transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;

  // Execute and save file if needed...
  if(!(iter % m_saveRateWall)) {
      computeWall();
  }

  if(!(iter % m_saveRateAero)) {
    computeAeroFVMCC();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file handle is opened here

  // unused // const CFuint nbEqs = PhysicalModel::getNbEq();
  // unused //  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  CFreal invRho;
  CFreal rhoK2;
  CFreal p;
  CFreal Cp;
  CFreal a2;
  CFreal Mach;
  CFreal T;

  Common::SafePtr<std::vector<CFuint> > nodesIdx = getCurrentTRS()->getNodesInTrs();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  const CFreal time = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();

  for (CFuint iNState = 0; iNState < nodesIdx->size(); ++iNState) {

    const RealVector& state = nstates[(*nodesIdx)[iNState]];
    Node *const node = nodes[(*nodesIdx)[iNState]];

    invRho = 1./(state)[0];
    rhoK2 = 0.5*((state)[1]*(state)[1] + (state)[2]*(state)[2])*invRho;
    p = gammaMinus1*((state)[3] - rhoK2);

    Cp = (p - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

    a2 = gamma*p/(state)[0];
    Mach = sqrt((((state)[1]*(state)[1])*invRho + ((state)[2]*(state)[2])*invRho)/a2);
    T = a2/(gamma*R);

    CFreal TDim = T * m_updateVarSet->getModel()->getTempRef();
    CFreal pDim = p * m_updateVarSet->getModel()->getPressRef();

    const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
    const CFreal rhoDim = (state)[0] * rhoRef;

    // Output to File
    fout
    << (*node)[XX]  << " "
    << (*node)[YY]  << " "
    << m_alphadeg   << " "
    << iter         << " "
    << time         << " "
    << pDim         << " "
    << Mach         << " "
    << Cp           << " "
    << TDim         << " "
    << rhoDim       << "\n";
  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::computeAeroFVMCC()
{
  CFAUTOTRACE;

//unused//  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  //Get the datahandle containing the boundary Normals
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  const std::string TRSName = trs->getName();
  const std::string socketName = TRSName + "-boundaryNormals";

  DataHandle<const CFreal*> boundaryNormals = m_sockets.
    getSocketSource<const CFreal*>(socketName)->getDataHandle();
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  const CFuint nbNodesInFace = 2;

  RealVector normal(0.0, dim);
  RealVector Cp(dim);
  Cp = 0.;

  Common::SafePtr<GeometricEntityPool<Framework::FaceTrsGeoBuilder> >
    geoBuilder = &m_faceTrsGeoBuilder;

  //set up the GeometricEntity builders
  m_faceTrsGeoBuilder.setup();

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  Common::SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  geoData.trs = getCurrentTRS();
  geoData.isBFace = true;

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  State tempState;

  cf_assert(boundaryNormals.size() == nbTrsFaces);

  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace)
  {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();

    // get the face normal
    for(CFuint iDim = 0; iDim < DIM_2D; iDim++)
    {
      normal[iDim] = -(boundaryNormals[iFace])[iDim];
    }

    std::vector<Node*>& nodes = *currFace->getNodes();
    cf_assert(nodes.size() == nbNodesInFace);

    vector <CFuint> nodeID(nbNodesInFace);
    for (CFuint s = 0; s < nbNodesInFace; ++s)
    {
      nodeID[s] = (*nodes[s]).getLocalID();
    }

    std::vector<RealVector> nodalStates(nbNodesInFace);
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

    for(CFuint iNode = 0; iNode < nbNodesInFace; iNode++){
      nodalStates[iNode].resize(nbEqs);
      nodalStates[iNode] = nstates[nodeID[iNode]];
    }

    //Compute the pressure on a face
    CFreal press = 0.;
    for (CFuint s = 0; s < nbNodesInFace; ++s)
    {
      tempState = nodalStates[s];
      CFreal invRho = 1./(tempState)[0];
      CFreal rhoK2 = 0.5*((tempState)[1]*(tempState)[1] + (tempState)[2]*(tempState)[2])*invRho;
      CFreal p = gammaMinus1*(tempState[3] - rhoK2);

      press += (m_pInf-p);
    }
    press /= nbNodesInFace;

    //Compute Cp (scaled by m, depending on ownership of face to processor)
    //press-pInf are switched to take into account the minus sign in normal
    Cp[XX] += press * normal[XX];
    Cp[YY] += press * normal[YY];

    // release the face
    geoBuilder->releaseGE();
  }

  // adimensionalize Cp
  Cp /= (0.5 * m_rhoInf * m_uInf * m_uInf);

  // project Cp
  ///@todo modify this for adimensional values!!
  m_lift = (sin(m_alpharad)*Cp[XX] + cos(m_alpharad)*Cp[YY]) / refLength;
  m_drag = (cos(m_alpharad)*Cp[XX] - sin(m_alpharad)*Cp[YY]) / refLength;

  m_xForceCoef = Cp[XX] / refLength;
  m_yForceCoef = Cp[YY] / refLength;

  //Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::prepareOutputFileWall()
{
  CFAUTOTRACE;

  cf_assert (!m_fileWall->isopen());
	
  using namespace boost::filesystem;
  path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
  file = PathAppender::getInstance().appendAllInfo( file );

  ofstream& fout = m_fileWall->open(file);

  fout << "TITLE = Wall_Values" << "\n";
  fout << "VARIABLES = x y Alpha Iter Time Pressure Mach Cp Temperature Density" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::prepareOutputFileAero()
{
  CFAUTOTRACE;

  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  convergenceFile << "TITLE = Aerodynamic_Coeficients"  << "\n";
  convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef FxCoef FyCoef" << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DConsComputeAeroFVMCC::updateOutputFileAero()
{
  CFAUTOTRACE;

  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendAllInfo( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  convergenceFile
  << subSysStatus->getNbIter()        << " "
  << subSysStatus->getCurrentTime()   << " "
  << m_alphadeg                       << " "
  << m_lift                           << " "
  << m_drag                           << " "
  << m_xForceCoef                     << " "
  << m_yForceCoef                     << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoeff

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

