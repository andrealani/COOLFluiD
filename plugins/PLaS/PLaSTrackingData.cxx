
#include "Common/BadValueException.hh"
#include "Common/GlobalReduce.hh"
#include "Environment/DirPaths.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "PLaS/PLaSModule.hh"
#include "PLaS/PLaSTrackingData.hh"
#include "PLaS/StgImplementation.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// Null MethodCommandProvider instance
MethodCommandProvider< NullMethodCommand< PLaSTrackingData >,PLaSTrackingData,PLaSModule > nullPLaSTrackingComProvider("Null");

/// Strategy for implementation, naked pointer (extern from StgImplementation)
StgImplementation* g_stg_implementation;

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::defineConfigOptions(Config::OptionList& options)
{
  // strategies
  options.addConfigOption< std::string >("StgImplementation","PLaS strategy for interface implementation (default \"StgImplementationStd\")");

  // PLaS input parameters
  options.addConfigOption< int >("NumMaxEnt","Maximum number of entities per process (default 10000)");
  options.addConfigOption< int >("NumIniEnt","Number of initially distributed entities (default 0)");
  options.addConfigOption< int >("IniDiamType","Type of diameter distribution (default 0)");
  options.addConfigOption< double >("IniDiam","Diameter of dispersed entities (default 4.0e-4)");
  options.addConfigOption< double >("IniDiamStd","Diameter standard deviation (default 7.55e-5)");
  options.addConfigOption< std::vector< double > >("IniVel","Initial velocity of dispersed entities (default < 0. 0. 0. >");
  options.addConfigOption< double >("IniTempDisp","Temperature of dispersed entities (default 293.15)");
  options.addConfigOption< int >("Material","Flag: Entity material (defines flow type) (default 7)");
  options.addConfigOption< int >("MomentumCoupl","Flag: Momentum coupling (default 0)");
  options.addConfigOption< int >("VolfracCoupl","Flag: Volume fraction coupling (default 0)");
  options.addConfigOption< int >("EnergyCoupl","Flag: Energy coupling (default 0)");
  options.addConfigOption< int >("CollisionModel","Flag: Collision model (default 0)");
  options.addConfigOption< int >("LiftForce","Flag: Lift force (only for bubbly flow) (default 0)");
  options.addConfigOption< int >("EvapModel","Flag: Evaporation model (only for droplet flow) (default 0)");
  options.addConfigOption< int >("SaturModel","Flag: Saturation model (only for bubbly flow) (default 0)");
  options.addConfigOption< int >("PerBnd","Flag: Periodic boundaries for dispersed entities (default 0)");
  options.addConfigOption< std::vector< double > >("GravVec","Gravity vector (default <>)");
  options.addConfigOption< std::string >("WriteStatsFilename","Write output statistics filename (default \"plas.txt\")");
  options.addConfigOption< std::string >("WriteTecplotFilename","Write output tecplot filename (default \"plas.plt\")");
  options.addConfigOption< std::string >("ConfFilename","Write configuration filename (default \"plas.conf\")");

  // PLaS flow solver parameters
  options.addConfigOption< int >("Restart","Restart from the solution provided (0: don't, default)");
  options.addConfigOption< double >("rhoCont","Primary (continuous) phase density (default 995.65, water at 20dC)");
  options.addConfigOption< double >("nuCont","Primary (continuous) phase kinematic viscosity [m2.s-1] (default 0.801e-6, water at 20dC)");
  options.addConfigOption< double >("dtEul","Primary (continuous) phase kinematic viscosity [s] (default 0.)");
  options.addConfigOption< double >("cpCont","Specific heat coefficient of the flow medium [J.g^-1.K^-1] (default 4.1855, water at 25dC -- wikipedia)");
  options.addConfigOption< double >("kCont","Thermal conductivity of the flow medium [W.m^-1.K^-1] (default 0.6, water -- wikipedia)");

  // boundary walls and production domains
  options.addConfigOption< std::vector< std::string > >("BoundaryWallsTRS","Boundary walls TRSs (default < \"\" >)");
  options.addConfigOption< std::vector< std::string > >("ProductionDomainsTRS","Production domains TRSs (default < \"\" >)");
  options.addConfigOption< std::vector< std::string > >("ProductionDomainsType","Production domains types, which can be \"line\", \"rectangle\" and \"ellipse\" (default < \"\" >)");
  options.addConfigOption< std::vector< double > >("ProductionDomainsMFlux","Production domains mass fluxes (default < >)");

  // variables names and default values
  options.addConfigOption< std::vector< std::string > >("VarVelocities","Velocity variables names (default <>)");
  options.addConfigOption< std::string >("VarPressure","Pressure variable name (default \"\")");
  options.addConfigOption< std::string >("VarTemperature","Temperature variable name (default \"\")");
  options.addConfigOption< std::vector< double > >("DefVelocities","Velocity default values (default < 0. 0. 0.>)");
  options.addConfigOption< double >("DefPressure","Pressure default value (default 0.)");
  options.addConfigOption< double >("DefTemperature","Temperature default value (default 298.15)");
}

//////////////////////////////////////////////////////////////////////////////

PLaSTrackingData::PLaSTrackingData(SafePtr< Method > owner) :
  DataProcessingData(owner),
  m_block(true),
  m_blockable(false)
{
  addConfigOptionsTo(this);

  // strategies
  m_stg_implementation_str = "StgImplementationStd";
  setParameter("StgImplementation",&m_stg_implementation_str);

  // initialize plas_data
  m_plasdata.numExtEnt = 0;

  // PLaS input parameters
  PLAS_INPUT_PARAM& ip = m_plasdata.ip;
  ip.numMaxEnt      = 10000;
  ip.numIniEnt      = 0;
  ip.iniDiamType    = 0;
  ip.iniDiam        = 4.0e-4;
  ip.iniDiamStd     = 7.55e-5;
  ip_iniVel.clear();
  ip.iniTempDisp    = 293.15;
  ip.material       = 7;
  ip.momentumCoupl  = 0;
  ip.volfracCoupl   = 0;
  ip.energyCoupl    = 0;
  ip.collisionModel = 0;
  ip.liftForce      = 0;
  ip.evapModel      = 0;
  ip.saturModel	    = 0;
  ip.perBnd         = 0;
  ip_gravVec.clear();
  ip_writeStatsFilename   = "plas.txt";
  ip_writeTecplotFilename = "plas.plt";
  ip_confFilename         = "plas.conf";
  setParameter("NumMaxEnt",&ip.numMaxEnt);
  setParameter("NumIniEnt",&ip.numIniEnt);
  setParameter("IniDiamType",&ip.iniDiamType);
  setParameter("IniDiam",&ip.iniDiam);
  setParameter("IniDiamStd",&ip.iniDiamStd);
  setParameter("IniVel",&ip_iniVel);
  setParameter("IniTempDisp",&ip.iniTempDisp);
  setParameter("Material",&ip.material);
  setParameter("MomentumCoupl",&ip.momentumCoupl);
  setParameter("VolfracCoupl",&ip.volfracCoupl);
  setParameter("EnergyCoupl",&ip.energyCoupl);
  setParameter("CollisionModel",&ip.collisionModel);
  setParameter("LiftForce",&ip.liftForce);
  setParameter("EvapModel",&ip.evapModel);
  setParameter("SaturModel",&ip.saturModel);
  setParameter("PerBnd",&ip.perBnd);
  setParameter("GravVec",&ip_gravVec);
  setParameter("WriteStatsFilename",&ip_writeStatsFilename);
  setParameter("WriteTecplotFilename",&ip_writeTecplotFilename);
  setParameter("ConfFilename",&ip_confFilename);

  // PLaS flow solver parameters
  PLAS_FLOWSOLVER_PARAM& fp = m_plasdata.fp;
  fp.restart = 0;
  fp.flowSolver = 0;     // (1): set at setup
  fp.numDim  = 0;        // (1)
  fp.numUnk  = 0;        // (1)
  fp.numNod  = 0;        // (1)
  fp.numElm  = 0;        // (1)
  fp.numBnd  = 0;        // (1)
  fp.rhoCont = 995.65;
  fp.muCont  = 797.52;   // (1)
  fp.nuCont  = 0.801e-6;
  fp.cpCont  = 4.1855;
  fp.kCont   = 0.6;
  fp.dtEul   = 0.;
  fp.domainVolume = 0.;  // (2): set by StgImplementation, on initialization
  fp.minElmVolume = 0.;  // (2)
  fp.maxElmVolume = 0.;  // (2)
  fp.iter         = 0;   // (3): set by StgImplementation, on every iteration
  fp.time         = 0.;  // (3)
  fp.writeOutput  = 0;   // (3)
  setParameter("Restart",&fp.restart);
  setParameter("rhoCont",&fp.rhoCont);
  setParameter("nuCont",&fp.nuCont);
  setParameter("dtEul",&fp.dtEul);
  setParameter("cpCont",&fp.cpCont);
  setParameter("kCont",&fp.kCont);

  // boundary walls and production domains
  m_boundarywalls_trs.clear();
  m_pdomains_trs.clear();
  m_pdomains_type.clear();
  m_pdomains_mflux.clear();
  setParameter("BoundaryWallsTRS",&m_boundarywalls_trs);
  setParameter("ProductionDomainsTRS",&m_pdomains_trs);
  setParameter("ProductionDomainsType",&m_pdomains_type);
  setParameter("ProductionDomainsMFlux",&m_pdomains_mflux);

  // variables names and default values
  m_vstr.clear();
  m_pstr = "";
  m_tstr = "";
  m_vdef.clear();
  m_pdef = 0.;
  m_tdef = 298.15;
  setParameter("VarVelocities",&m_vstr);
  setParameter("VarPressure",&m_pstr);
  setParameter("VarTemperature",&m_tstr);
  setParameter("DefVelocities",&m_vdef);
  setParameter("DefPressure",&m_pdef);
  setParameter("DefTemperature",&m_tdef);
}

//////////////////////////////////////////////////////////////////////////////

PLaSTrackingData::~PLaSTrackingData()
{
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  DataProcessingData::configure(args);

  // PLaS strategies
  configureStrategy< StgImplementation,PLaSTrackingData >( args,
    m_stg_implementation, m_stg_implementation_str, m_stg_implementation_str,
    SharedPtr< PLaSTrackingData >(this) );
  CFLog(INFO,"PLaSTrackingData: implementation strategy: \"" << m_stg_implementation->getName() << "\"\n");

  // production domains
  m_plasdata.ip.numProdDom = (int) m_pdomains_trs.size();
  if ( m_pdomains_trs.size()!=m_pdomains_type.size() ||
       m_pdomains_trs.size()!=m_pdomains_mflux.size() ) {
    std::string m("PLaSTrackingData: \"ProductionDomainsTRS/Type/MFlux\" not consistent!");
    CFLog(ERROR,m << "\n");
    throw BadValueException(FromHere(),m);
  }
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::setup()
{
  CFAUTOTRACE;
  DataProcessingData::setup();


  // get nodes/states, and number of dimensions/variables
  DataHandle< Node*,GLOBAL >  h_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  DataHandle< State*,GLOBAL > h_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  const CFuint nbDim = h_nodes[0]->size();
  const CFuint nbVar = h_states[0]->size();


  // get 'inner'/'boundary' regions
  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  std::vector< SafePtr< TopologicalRegionSet > > itrs = MeshDataStack::getActive()->getFilteredTrsList("inner");
  cf_always_assert_desc("no appropriate 'inner' TRS found",itrs.size());


  // m_plasdata ip/fp parameters to adjust
  PLAS_INPUT_PARAM&      ip = m_plasdata.ip;
  PLAS_FLOWSOLVER_PARAM& fp = m_plasdata.fp;


  CFLog(INFO,"PLaSTrackingData: set vector parameters...\n");
  if (!ip_gravVec.size())  ip_gravVec.assign(nbDim,0.);
  if (!ip_iniVel.size())   ip_iniVel.assign(nbDim,0.);
  if (!m_vstr.size())      m_vstr.assign(nbDim,"");
  if (!m_vdef.size())      m_vdef.assign(nbDim,0.);
  if (ip_gravVec.size()!=nbDim || ip_iniVel.size()!=nbDim ||
      m_vstr.size()!=nbDim     || m_vdef.size()!=nbDim) {
    std::string m("\"GravVec\", \"IniVel\", \"VarVelocities\" and \"DefVelocities\" sizes should be " + StringOps::to_str(nbDim));
    CFLog(ERROR,m << '\n');
    throw BadValueException(FromHere(),m);
  }
  for (CFuint i=0; i<nbDim; ++i) {
    ip.gravVec[i] = ip_gravVec[i];
    ip.iniVel[i]  = ip_iniVel[i];
  }
  CFLog(VERBOSE,"PLaSTrackingData: set vector parameters.\n");


  CFLog(INFO,"PLaSTrackingData: set variables...\n");
  m_vidx.assign(nbDim,-1);
  for (CFuint i=0; i<nbDim; ++i) {
    m_vidx[i] = setVariableIndex(m_vstr[i],"Velocity");
    CFLog(INFO,"  Velocity: " << (m_vidx[i]<0?
      StringOps::to_str(m_vdef[i]) : std::string("(variable)")) << '\n');
  }
  m_pidx = setVariableIndex(m_pstr,"Pressure");
  CFLog(INFO,"  Pressure: " << (m_pidx<0?
    StringOps::to_str(m_pdef) : std::string("(variable)")) << '\n');
  m_tidx = setVariableIndex(m_tstr,"Temperature");
  CFLog(INFO,"  Temperature: " << (m_tidx<0?
    StringOps::to_str(m_tdef) : std::string("(variable)")) << '\n');
  CFLog(VERBOSE,"PLaSTrackingData: set variables.\n");


  CFLog(INFO,"PLaSTrackingData: set filenames...\n");
  ip_writeStatsFilename   = setFilename(ip_writeStatsFilename,ip.writeStatsFilename);
  ip_writeTecplotFilename = setFilename(ip_writeTecplotFilename,ip.writeTecplotFilename);
  ip_confFilename         = setFilename(ip_confFilename,ip.confFilename);
  CFLog(VERBOSE,"PLaSTrackingData: set filenames.\n");


  CFLog(INFO,"PLaSTrackingData: set boundaries...\n");
  {
    // reset boundary walls and production domains
    m_boundary = btrs;
    m_boundary_is_wall.assign(btrs.size(),false);
    m_boundary_is_pdomain.assign(btrs.size(),false);
    m_pdomains.clear();

    // set TRSs as walls and/or production domains
    for (CFuint b=0; b<btrs.size(); ++b) {

      // check if it is a wall
      for (CFuint t=0; t<m_boundarywalls_trs.size(); ++t)
        if (btrs[b]->getName()==m_boundarywalls_trs[t]) {
          m_boundary_is_wall[b] = true;
          break;
        }

      // check if it is a production domain, get bounding box and set description
      for (CFuint t=0; t<m_pdomains_trs.size(); ++t)
        if (btrs[b]->getName()==m_pdomains_trs[t]) {

          RealVector bb0( 1.e99,nbDim);
          RealVector bb1(-1.e99,nbDim);
          std::vector< CFuint >& nodes = *(btrs[b]->getNodesInTrs());
          for (std::vector< CFuint >::const_iterator in=nodes.begin(); in!=nodes.end(); ++in)
            for (CFuint d=0; d<nbDim; ++d) {
              bb0[d] = std::min(bb0[d],(*h_nodes[*in])[d]);
              bb1[d] = std::max(bb1[d],(*h_nodes[*in])[d]);
            }
          for (CFuint d=0; d<nbDim; ++d) {
            GlobalReduceOperation< GRO_MIN >(&bb0[d],&bb0[d]);
            GlobalReduceOperation< GRO_MAX >(&bb1[d],&bb1[d]);
          }
          PE::GetPE().setBarrier();

          std::ostringstream p;
          p << (m_pdomains_type[t]=="line"?      '1' :
               (m_pdomains_type[t]=="rectangle"? '2' :
               (m_pdomains_type[t]=="ellipse"?   '3' : '0' )))
            << ' ' << bb0 << bb1 << m_pdomains_mflux[t];
          m_pdomains.push_back(p.str());

          m_boundary_is_pdomain[b] = true;
          break;
        }

      // report
      CFLog(INFO,"PLaSTrackingData: boundary \"" << btrs[b]->getName() << '"'
        << (m_boundary_is_wall[b]?    " (wall)":"")
        << (m_boundary_is_pdomain[b]? " (p. domain)":"") << '\n');
    }
  }
  CFLog(VERBOSE,"PLaSTrackingData: set boundaries.\n");


  CFLog(INFO,"PLaSTrackingData: set flow solver parameters...\n");
  fp.flowSolver = FLOWSOLVER_COOLFLUID;
  fp.numDim = (int) nbDim;
  fp.numUnk = (int) nbVar;
  fp.numNod = (int) h_nodes.size();
  fp.numElm = (int) itrs[0]->getLocalNbGeoEnts();
  fp.numBnd = (int) btrs.size();
  fp.muCont = fp.nuCont*fp.rhoCont;
  CFLog(VERBOSE,"PLaSTrackingData: set flow solver parameters.\n");


  if (!PE::GetPE().GetRank()) {
    CFLog(INFO,"PLaSTrackingData: write configuration file...\n");
    writeConfiguration(ip_confFilename,ip);
    CFLog(VERBOSE,"PLaSTrackingData: write configuration file.\n");
  }
  PE::GetPE().setBarrier();
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::setFlowProperties(double _rho, double _nu, double _dt, double _cpCont, double _kCont, int _numunk)
{
  /*
   * note: there is a risk that setting dt=0. will provoke corruption inside
   * PLaS. this is not a problem because dt is only used at execution time
   * (where it is set propperly) but if PLaS used dt at initialization, a fix
   * is necessary.
   */
  PLAS_FLOWSOLVER_PARAM& fp = m_plasdata.fp;
  fp.rhoCont = _rho;
  fp.muCont  = _nu*_rho;
  fp.nuCont  = _nu;
  fp.dtEul   = _dt;  // overlapped at setFlow/IterationProperties
  fp.cpCont  = _cpCont;
  fp.kCont   = _kCont;
  fp.numUnk  = _numunk;
  CFLog(INFO,"PLaSTrackingData: density:              " << _rho    << " [kg/m3]\n");
  CFLog(INFO,"PLaSTrackingData: kinematic viscosity:  " << _nu     << " [m2/s]\n");
  CFLog(INFO,"PLaSTrackingData: eulerian time scale:  " << _dt     << " [s]\n");
  CFLog(INFO,"PLaSTrackingData: specific heat coeff.: " << _cpCont << " [J.g^-1.K^-1]\n");
  CFLog(INFO,"PLaSTrackingData: thermal conductivity: " << _kCont  << " [W.m^-1.K^-1]\n");
  CFLog(INFO,"PLaSTrackingData: number of unknowns:   " << _numunk << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::setIterationProperties(int _iter, double _time, double _dt, int _output, int _numExtEnt, double *_extEntPos, double *_extEntVel, double *_extEntTemp, double *_extEntDiam)
{
  /*
   * TODO: separate imposing external entities, it is difficult due to the
   * different implementation strategies
   */
  PLAS_FLOWSOLVER_PARAM& fp = m_plasdata.fp;
  fp.iter        = _iter;
  fp.time        = _time;
  fp.dtEul       = _dt;  // overlapped at setFlow/IterationProperties
  fp.writeOutput = _output;
  CFLog(INFO,"PLaSTrackingData: iteration:           " << _iter      << '\n');
  CFLog(INFO,"PLaSTrackingData: time:                " << _time      << " [s]\n");
  CFLog(INFO,"PLaSTrackingData: eulerian time scale: " << _dt        << " [s]\n");
  CFLog(INFO,"PLaSTrackingData: output?              " << (_output? "true":"false") << '\n');

  m_plasdata.numExtEnt  = _numExtEnt;
  m_plasdata.extEntPos  = _extEntPos;
  m_plasdata.extEntVel  = _extEntVel;
  m_plasdata.extEntTemp = _extEntTemp;
  m_plasdata.extEntDiam = _extEntDiam;
  CFLog(INFO,"PLaSTrackingData: external entities:   " << _numExtEnt << '\n');
}

//////////////////////////////////////////////////////////////////////////////

int PLaSTrackingData::setVariableIndex(const std::string& var, const std::string& description)
{
  if (!var.length() || var=="?")
    return -1;

  // get variables names (weird stuff!)
  //FIXME hardcoded provider!
  //FIXME isn't this stupid because there will temporarily exist multiple ConvectiveVarSet?
  SafePtr< PhysicalModel > pm = PhysicalModelStack::getActive();
  SafePtr< BaseTerm > pm_ct = pm->getImplementor()->getConvectiveTerm();
  SelfRegistPtr< ConvectiveVarSet > updatevarset = Environment::Factory< ConvectiveVarSet >::getInstance().getProvider("PhysicalModelDummyPrim")->create(pm_ct);
  cf_assert(updatevarset.isNotNull());
  const std::vector< std::string > varnames = updatevarset->getVarNames();

  int idx = -1;
  for (CFuint i=0; i<varnames.size(); ++i)
    if (var==varnames[i]) {
      idx = i;
      break;
    }

  if (idx<0) {
    std::string m("PLaSTrackingData: Variable for " + description + " (\"" + var + "\") not found");
    CFLog(ERROR,m << '\n');
    throw BadValueException(FromHere(),m);
  }

  return idx;
}

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_MPI
void PLaSTrackingData::setParallelDataStructures()
{
  PE::GetPE().setBarrier();
  CFLog(INFO,"PLaSTracking: set global parallel information...\n");

  // get nodes and states DataHandle (to access local/global indices)
  DataHandle< Node*,GLOBAL >  h_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();
  DataHandle< State*,GLOBAL > h_states = MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  CFLog(INFO,"PLaSTracking: 1/3: set local parallel information\n");

  // MPI node rank/nodes involved, communicator and int MPI_Datatype handle
  const int mpi_rank = (int) PE::GetPE().GetRank();
  const int mpi_size = (int) PE::GetPE().GetProcessorCount();
  const MPI_Comm mpi_comm = PE::GetPE().GetCommunicator();
  const MPI_Datatype mpi_int = MPIDataTypeHandler::GetType< int >();

  //  get states parallel information to get ghost send/receive lists
  const std::vector< std::vector< unsigned int > >& gsend = h_states.getGhostSendList();
  const std::vector< std::vector< unsigned int > >& grecv = h_states.getGhostReceiveList();

  // count total number of partition ghost nodes, and sum on all partitions
  int Ngnodes = 0;
  int TNgnodes = 0;
  for (int r=0; r<mpi_size; ++r)
    Ngnodes += grecv[r].size();
  GlobalReduceOperation< GRO_SUM >(&Ngnodes,&TNgnodes);

  // number of elements sent by each process (light communication) and
  // displacements at which to place the incoming data
  std::vector< int > recvcounts(mpi_size,-1);
  std::vector< int > displs(mpi_size,0);
  int recvl = Ngnodes;
  MPI_Allgather(&recvl,1,mpi_int,&recvcounts[0],1,mpi_int,mpi_comm);
  for (int r=1; r<mpi_size; ++r)
    displs[r] = displs[r-1] + recvcounts[r-1];

  // allocate and set local information on global structures
  // (could be local but reusing the global ones is more efficient)
  m_ghosts_srank.resize(TNgnodes);
  m_ghosts_rrank.resize(TNgnodes);
  m_ghosts_sindex.resize(TNgnodes);
  m_ghosts_rindex.resize(TNgnodes);
  int i = displs[mpi_rank];
  for (int r=0; r<mpi_size; ++r) {
    for (unsigned int n=0; n<grecv[r].size(); ++n, ++i) {
      m_ghosts_srank[i] = r;
      m_ghosts_rrank[i] = mpi_rank;
      m_ghosts_sindex[i] = h_nodes[ grecv[r][n] ]->getGlobalID();
      m_ghosts_rindex[i] = grecv[r][n];
    }
  }

  CFLog(INFO,"PLaSTracking: 2/3: broadcast global indices\n");

  // synchronize between all MPI nodes (heavy communication!)
  MPI_Allgatherv( &m_ghosts_srank[displs[mpi_rank]],  Ngnodes, mpi_int, &m_ghosts_srank[0],  &recvcounts[0], &displs[0], mpi_int, mpi_comm );
  MPI_Allgatherv( &m_ghosts_rrank[displs[mpi_rank]],  Ngnodes, mpi_int, &m_ghosts_rrank[0],  &recvcounts[0], &displs[0], mpi_int, mpi_comm );
  MPI_Allgatherv( &m_ghosts_sindex[displs[mpi_rank]], Ngnodes, mpi_int, &m_ghosts_sindex[0], &recvcounts[0], &displs[0], mpi_int, mpi_comm );
  MPI_Allgatherv( &m_ghosts_rindex[displs[mpi_rank]], Ngnodes, mpi_int, &m_ghosts_rindex[0], &recvcounts[0], &displs[0], mpi_int, mpi_comm );

  CFLog(INFO,"PLaSTracking: 3/3: find local indices from global\n");

  // convert receiving node -global- index to sending node -local- index
  // (heavy communication!)
  // *1: this is inefficient but it guarantees everybody has all information
  // *2: only the processor with the node as local does the search
  for (int i=0; i<TNgnodes; ++i) {
    const CFuint gi = (CFuint) m_ghosts_sindex[i];  // to be found
    const int rr = m_ghosts_rrank[i];               // processor receiving this node
    int sindex = -1;
    if (m_ghosts_srank[i]==mpi_rank) {
      for (int n=0; n<(int) gsend[rr].size(); ++n) {
        if (h_nodes[ gsend[rr][n] ]->getGlobalID()==(CFuint) gi) {
          sindex = h_nodes[ gsend[rr][n] ]->getLocalID();
          break;
        }
      }
    }
    GlobalReduceOperation< GRO_MAX >(&sindex,&m_ghosts_sindex[i]);
    cf_assert_desc("global index not found!",m_ghosts_sindex[i]>-1);
  }

  CFLog(VERBOSE,"PLaSTracking: set global parallel information.\n");
  PE::GetPE().setBarrier();
}
#endif

//////////////////////////////////////////////////////////////////////////////

std::string PLaSTrackingData::setFilename(const std::string& f, char*& fchar)
{
  using boost::filesystem::path;
  using boost::filesystem::exists;

  // - build path with ResultsDir and filename, remove if it exists
  const path p = Environment::DirPaths::getInstance().getResultsDir()/f;
  if (exists(p))
    remove(p);

  // - copy to C-style string including null-termination and return as string
  fchar = new char[p.string().length()+1];
  strcpy(fchar,p.string().c_str());
  return p.string();
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::writeConfiguration(const std::string& f, const PLAS_INPUT_PARAM& p)
{
  CFAUTOTRACE;
  using std::endl;

  std::ofstream fp(f.c_str(),std::ios::trunc);
  fp << "Maximum number of entities:" << endl
     << p.numMaxEnt << endl
     << "Number of initially distributed entities:" << endl
     << p.numIniEnt << endl
     << "Number of production domains:" << endl
     << p.numProdDom << endl
     << "Production domains ([1] = line, [2] = rectangle, [3] = ellipse):" << endl
     << "TYPE X0 Y0 Z0 X1 Y1 Z1 MASSFLUX" << endl;
  for (CFuint i=0; i<m_pdomains.size(); ++i)
    fp << m_pdomains[i] << endl;
  fp << "Diameter of dispersed entities ([0] = constant, [1] = normal, [2] = log-normal):" << endl
     << "TYPE MEAN STD" << endl
     << p.iniDiamType << ' ' << p.iniDiam << ' ' <<p.iniDiamStd << endl
     << "Initial relative velocity of dispersed entities:" << endl
     << "U V W" << endl
     << p.iniVel[0] << ' ' << p.iniVel[1] << ' ' << p.iniVel[2] << endl
     << "Initial temperature of dispersed entities:" << endl
     << p.iniTempDisp << endl
     << "Entity material ([1] = Cu, [2] = C8H8, [3] = H2O, [4] = nHeptane, [5] = H2, [6] = O2, [7] = Air):" << endl
     << p.material << endl
     << "Momentum back-coupling ([0] = none, [1] = PIC, [2] = Projection):" << endl
     << p.momentumCoupl << endl
     << "Volume fraction back-coupling ([0] = no, [1] = yes):" << endl
     << p.volfracCoupl << endl
     << "Energy coupling ([0] = no, [1] = one-way):" << endl
     << p.energyCoupl << endl
     << "Collision model ([0] = none, [1] = uncorrelated, [2] = correlated):" << endl
     << p.collisionModel << endl
     << "Slip-shear lift force for bubbly flow ([0] = no, [1] = yes):" << endl
     << p.liftForce << endl
     << "Thin film evaporation model for droplet flow ([0] = no, [1] = yes):" << endl
     << p.evapModel << endl
     << "Saturation model for bubbly flow ([0] = no, [1] = yes):" << endl
     << p.saturModel << endl
     << "Periodic boundaries ([0] = no, [1] = yes):" << endl
     << p.perBnd << endl
     << "Gravity:" << endl
     << "X Y Z" << endl
     << p.gravVec[0] << ' ' << p.gravVec[1] << ' ' << p.gravVec[2] << endl
     << "Write output statistics filename:" << endl
     << p.writeStatsFilename << endl
     << "Write output tecplot filename:" << endl
     << p.writeTecplotFilename << endl
     << "Configuration filename:" << endl
     << p.confFilename << endl
     << endl;
  fp.close();
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::PLaS_Init()
{
  CFAUTOTRACE;

  if (PE::GetPE().IsParallel()) {
    setParallelDataStructures();
  }

  g_stg_implementation = getStgImplementationNaked();
  ::initPLaS(&m_plasdata);
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::PLaS_Run()
{
  CFAUTOTRACE;

  g_stg_implementation = getStgImplementationNaked();
  if (!m_blockable || !m_block)
    ::runPLaS(&m_plasdata);
}

//////////////////////////////////////////////////////////////////////////////

void PLaSTrackingData::PLaS_Terminate()
{
  CFAUTOTRACE;

  g_stg_implementation = getStgImplementationNaked();
  ::terminatePLaS(&m_plasdata);

  delete [] m_plasdata.ip.writeStatsFilename;
  delete [] m_plasdata.ip.writeTecplotFilename;
  delete [] m_plasdata.ip.confFilename;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

