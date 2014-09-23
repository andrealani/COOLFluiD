
#include <cmath>
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/BCBulk.hh"
#include "Muffin/BCElectrode.hh"
#include "Muffin/BCFile.hh"
#include "Muffin/BCFunction.hh"
#include "Muffin/SystemMITReM.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemMITReM,MuffinData,MuffinModule > cSystemMITReMProvider("MITReM");

//////////////////////////////////////////////////////////////////////////////

SystemMITReM::SystemMITReM(const std::string& name) :
    System(name),
    s_mn_volume("NodalVolume"),              // socket sinks
    s_mn_voidfraction("NodalVoidFraction"),  // ...
    s_jx("CurrentDensityX"),              // socket sources
    s_jy("CurrentDensityY"),              // ...
    s_jz("CurrentDensityZ"),              // ...
    s_conductivity("Conductivity"),       // ...
    s_gasonsurf("GasOnSurface"),          // ...
    Nions(0),
    m_velocity_i(-1),
    m_magfield_i(-1)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemMITReM");
  addConfigOptionsTo(this);

  // general
  m_voidfraction = false;
  m_linearization = "Newton";
  m_linearization_switch = 666;
  m_negative_concentrations.clear();
  m_refelectrode_xyz.clear();
  setParameter("UseVoidFraction",&m_voidfraction);
  setParameter("LinearizationSwitchToNewton",&m_linearization_switch);
  setParameter("Linearization",&m_linearization);
  setParameter("NegativeConcentrations",&m_negative_concentrations);
  setParameter("ReferenceElectrodePoint",&m_refelectrode_xyz);

  // mhd
  m_magfield_h = 0.;
  setParameter("HartmannReferenceLength",&m_magfield_h);

  // library options
  m_chemistry_database = "database.ec.xml";
  m_chemistry_label    = "";
  m_chargeconservation = true;
  m_swap               = true;
  m_convection_str     = "N";          // Empty | Galerkin | MDC (2D) | LDA | N
  m_diffusion_str      = "Galerkin";   // Empty | Galerkin | MDC
  m_migration_str      = "Galerkin";   // Empty | Galerkin | MDC
  m_magnetic_str       = "Empty";      // Empty | Galerkin
  m_homreaction_str    = "Galerkin_Diagonal";  // Empty | Galerkin | Galerkin_Diagonal | MDC
  m_electrostatics_str = "Galerkin";   // Empty | Galerkin | MDC
  m_time_str           = "Empty";      // Empty | Galerkin | MDC
  m_elecreaction_str   = "Pointwise";  // Empty | Galerkin | MDC | Pointwise
  m_gasreaction_str    = "Galerkin";   // Empty | Galerkin
  setParameter("ChemistryDatabase",&m_chemistry_database);
  setParameter("ChemistryLabel",&m_chemistry_label);
  setParameter("DoChargeConservation",&m_chargeconservation);
  setParameter("DoSwapEquations",&m_swap);
  setParameter("ConvectionScheme",&m_convection_str);
  setParameter("DiffusionScheme",&m_diffusion_str);
  setParameter("MigrationScheme",&m_migration_str);
  setParameter("MagneticScheme",&m_magnetic_str);
  setParameter("HomReactionScheme",&m_homreaction_str);
  setParameter("ElectrostaticsScheme",&m_electrostatics_str);
  setParameter("TimeScheme",&m_time_str);
  setParameter("ElecReactionScheme",&m_elecreaction_str);
  setParameter("GasReactionScheme",&m_gasreaction_str);
}

//////////////////////////////////////////////////////////////////////////////

SystemMITReM::~SystemMITReM()
{
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

  // general
  options.addConfigOption< bool >("UseVoidFraction","If void fraction is to be used, or not (default)");
  options.addConfigOption< std::string >("Linearization","Jacobian matrix linearization method, which can be \"Picard\" or \"Newton\" (default)");
  options.addConfigOption< CFuint >("LinearizationSwitchToNewton","Newton linearization forced after this number of iterations (default 666)");
  options.addConfigOption< std::string >("NegativeConcentrations","How to treat negative concentrations (one of \"Ignore\" (default), \"SetZero\" or \"SetBulk\")");
  options.addConfigOption< std::vector< CFreal > >("ReferenceElectrodePoint","Reference electrode position to monitor field potential (default < 0. 0. 0. >)");

  // mhd
  options.addConfigOption< double >("HartmannReferenceLength","Reference length (height of channel) to calculate Hartmann number (M) (default 0.)");

  // library options
  options.addConfigOption< std::string >("ChemistryDatabase","MITReM chemistry database file (default \"database.ec.xml\")");
  options.addConfigOption< std::string >("ChemistryLabel","MITReM chemistry system label (default \"\")");
  options.addConfigOption< bool >("DoChargeConservation","If charge flux is to be assembled in place of one mass balance (default \"true\")");
  options.addConfigOption< bool >("DoSwapEquations","If first and last equations are to be swapped (default \"true\")");
  options.addConfigOption< std::string >("ConvectionScheme","Convection scheme (default \"N\")");
  options.addConfigOption< std::string >("DiffusionScheme","Diffusion scheme (default \"Galerkin\")");
  options.addConfigOption< std::string >("MigrationScheme","Migration scheme (default \"Galerkin\")");
  options.addConfigOption< std::string >("MagneticScheme","Magnetic scheme (default \"Empty\")");
  options.addConfigOption< std::string >("HomReactionScheme","Homogeneous reactions scheme (default \"Galerkin\")");
  options.addConfigOption< std::string >("ElectrostaticsScheme","Electrostatics scheme (default \"Galerkin\")");
  options.addConfigOption< std::string >("TimeScheme","Time scheme (default \"Empty\")");
  options.addConfigOption< std::string >("ElecReactionScheme","Electrode reactions scheme (default \"Galerkin\")");
  options.addConfigOption< std::string >("GasReactionScheme","Gas-production reactions scheme (default \"Empty\")");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  System::configure(args);

  // adjust files paths
  {
    using boost::filesystem::path;
    using boost::filesystem::exists;
    path work = Environment::DirPaths::getInstance().getWorkingDir();
    m_chemistry_database = path(work/path(m_chemistry_database)).string();

    // m_chemistry_database must be a valid file name
    if (!exists(m_chemistry_database))
      err("\"ChemistryDatabase \"" + m_chemistry_database + "\" not found");
  }


  log("set MITReM...");
  m_mitrem = new MITReM(m_chemistry_database,m_chemistry_label);
  Nions = (int) m_mitrem->getNIons();
  if (Nions==1 && m_chargeconservation)
    err( "charge conservation equation cannot be applied to 1 ion only");
  ver("set MITReM.");


  log("set bulk concentrations...");
  m_bulk.assign(m_mitrem->getNIons(),0.);
  for (CFuint i=0; i<m_mitrem->getNIons(); ++i) {
    m_bulk[i] = m_mitrem->getIonInletConcentration(i);
    log(m_mitrem->getIonLabel(i) + ": " + StringOps::to_str(m_bulk[i]));
  }
  ver("set bulk concentrations.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::setup()
{
  CFAUTOTRACE;
  System::setup();
  MuffinData& d = getMethodData();


  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size() : 0);

  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");


  // variable types and number of ions
  for (int i=iv; i<iv+Nsys; ++i)
    d.m_vartypes[i] = VCONCENTRATION;
  d.m_vartypes[iv+Nsys-1] = VPOTENTIAL;
  if (Nions!=Nsys-1)
    err( "number of equations in system (" + StringOps::to_str(Nsys) + ") must be one more than number of ions (" + StringOps::to_str(Nions) + ")" );


  log("set ElementMatrixAssembler...");
  m_assembler = new ElementMatrixAssembler(
    (nbDim==2? "2D":"3D"), m_mitrem,
    m_convection_str, m_diffusion_str, m_migration_str, m_magnetic_str,
    m_homreaction_str, m_electrostatics_str, m_time_str, m_elecreaction_str,
    m_gasreaction_str,
    m_chargeconservation, m_swap );
  ver("set ElementMatrixAssembler.");


  log("set electrodes...");
  m_electrodes.clear();
  for (CFuint i=0; i<d.m_vcomm_bc.size(); ++i) {
    if (d.m_vcomm_bc[i]->hasTag("Muffin::BCElectrode")) {
      try {
        SafePtr< BCElectrode > p(d.m_vcomm_bc[i].d_castTo< BCElectrode >());
        m_electrodes.push_back(p);
      }
      catch (FailedCastException& e) {
        log("couldn't cast \"" + d.m_vcomm_bc[i]->getName() + "\" to BCElectrode");
        throw;
      }
    }
  }
  ver("set electrodes.");


  log("set electrode and gas reactions...");
  std::vector< std::string > ereactions;
  std::vector< std::string > greactions;
  for (CFuint r=0; r<m_mitrem->getNElecReactions(); ++r)
    ereactions.push_back(m_mitrem->getElecReactionLabel(r));
  for (CFuint r=0; r<m_mitrem->getNGasReactions(); ++r)
    greactions.push_back(m_mitrem->getGasReactionLabel(r));
  for (CFuint e=0; e<m_electrodes.size(); ++e)
    m_electrodes[e]->setupReactions(*this,ereactions,greactions);
  ver("set electrode and gas reactions.");


  log("set gas on surface tracking...");
  h_gasonsurf.resize(btrs.size());
  for (CFuint i=0; i<btrs.size(); ++i) {
    if (!btrs[i]->hasTag("Muffin::GasProduction"))
      continue;
    const ConnectivityTable< CFuint >& faces = *btrs[i]->getGeo2NodesConn();

    const GasOnSurface gempty = { 0., 0., 0., 0., true };
    h_gasonsurf[i].assign(faces.nbRows(),gempty);
    log("gas on surface elements at \""+btrs[i]->getName()+"\": "+StringOps::to_str(h_gasonsurf[i].size()));

    // set boundary elements size and if they are local
    for (CFuint f=0; f<faces.nbRows(); ++f) {
      const Node& n0 = *h_nodes[faces(f,0)];
      const Node& n1 = *h_nodes[faces(f,1)];
      const Node& n2 = (nbDim>2? *h_nodes[faces(f,2)]:n1);

      h_gasonsurf[i][f].isLocal = n0.isParUpdatable() &&
                                  n1.isParUpdatable() &&
                                  n2.isParUpdatable();
      if (nbDim==3) {
        const RealVector v1 = n1 - n0;
        const RealVector v2 = n2 - n0;
        const double S = 0.5 * (
          v1[1]*v2[2] - v1[2]*v2[1] +
          v1[2]*v2[0] - v1[0]*v2[2] +
          v1[0]*v2[1] - v1[1]*v2[0] );
        h_gasonsurf[i][f].S = (S>0.? S:-S);
      }
      else if (nbDim==2) {
        const RealVector v1 = n1 - n0;
        h_gasonsurf[i][f].S = v1.norm2();
      }
    }

  }
  ver("set gas on surface tracking.");


  // set linearization
  if (m_linearization!="Picard" && m_linearization!="Newton")
    err("\"Linearization\" must be either \"Picard\" or \"Newton\"");


  // get external variables
  m_velocity_i = d.getVariableIndex(VVELOCITY);
  m_magfield_i = d.getVariableIndex(VMAGFIELD);
  log("velocity field: " + (m_velocity_i>=0? d.m_varnames[m_velocity_i] : "not found"));
  log("magnetic field: " + (m_magfield_i>=0? d.m_varnames[m_magfield_i] : "not found"));


  // reset current density and conductivity DataHandles
  DataHandle< std::vector< double > > h_jx = s_jx.getDataHandle();
  DataHandle< std::vector< double > > h_jy = s_jy.getDataHandle();
  DataHandle< std::vector< double > > h_jz = s_jz.getDataHandle();
  DataHandle< double > h_kappa = s_conductivity.getDataHandle();
  const CFuint nbCells = geo2nodes->nbRows();
  h_jx.resize(Nions+1);
  h_jy.resize(Nions+1);
  h_jz.resize(Nions+1);
  for (int i=0; i<Nions+1; ++i) {
    h_jx[i].assign(nbCells,0.);
    h_jy[i].assign(nbCells,0.);
    h_jz[i].assign(nbCells,0.);
  }
  h_kappa.resize(h_nodes.size());
  h_kappa = 0.;


  log("reset current and current density (at the boundaries)...");
  m_vi.assign(btrs.size(),0.);
  m_vj.resize(m_mitrem->getNElecReactions());
  for (unsigned r=0; r<m_mitrem->getNElecReactions(); ++r)
    m_vj[r].assign(h_nodes.size(),0.);
  ver("reset current and current density (at the boundaries).");


  log("set reference electrode node...");
  // set to closest point to coordinates (negative means no monitoring)
  m_refelectrode_node = -1;
  if (m_refelectrode_xyz.size()) {
    if (m_refelectrode_xyz.size()!=nbDim)
      err("\"ReferenceElectrodePoint\" needs to be same size as dimensions");

    // seek closest node
    const RealVector refposition(nbDim,&m_refelectrode_xyz[0]);
    CFreal d2_min = 1.e99;
    for (CFuint n=0; n<h_nodes.size(); ++n) {
      if (h_nodes[n]->isParUpdatable()) {
        const CFreal d2 = RealVector(refposition - *h_nodes[n]).sqrNorm();
        if (d2<d2_min) {
          d2_min = d2;
          m_refelectrode_node = (int) n;
        }
      }
    }
    if (m_refelectrode_node>=0) {
      std::ostringstream msg;
      msg << *h_nodes[m_refelectrode_node];
      log("reference electrode position: " + msg.str());
      log("reference electrode node: " + StringOps::to_str(m_refelectrode_node));
    }
    else {
      log("reference electrode position not found");
    }

  }
  ver("set reference electrode node.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::unsetup()
{
  CFAUTOTRACE;
  System::unsetup();

  delete m_assembler;
  delete m_mitrem;
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::execute()
{
  CFAUTOTRACE;


  log("System iteration cycle...");
  System::execute();
  ver("System iteration cycle.");


  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const CFuint nbStates = h_states.size();
  {
    std::ostringstream msg;
    if (m_refelectrode_node>=0)  msg << (*h_states[m_refelectrode_node])[iv+Nions];
    else                         msg << "not set";
    log("reference electrode potential [V]: " + msg.str());
  }


  log("correcting negative concentrations...");
  const bool setzero = (m_negative_concentrations=="SetZero");
  const bool setbulk = (m_negative_concentrations=="SetBulk");
  std::vector< int > nneg(Nions,0);
  for (CFuint n=0; n<nbStates; ++n) {
    for (int i=0; i<Nions; ++i) {
      if ((*h_states[n])[iv+i]<0.) {
        ++nneg[i];
        if (setzero) {
          for (int j=0; j<Nions; ++j)
            (*h_states[n])[iv+j] = 0.;
          break;
        }
        if (setbulk) {
          for (int j=0; j<Nions; ++j)
            (*h_states[n])[iv+j] = m_bulk[j];
          break;
        }
      }
    }
  }
  for (int j=0; j<Nions; ++j) {
    if (nneg[j]) {
      std::string s;
      for (int i=0; i<Nions; ++i)
        s += " " + StringOps::to_str(nneg[i]);
      log("negative concentration nodes: " + s +
        (setzero? " set to zero":"") +
        (setbulk? " set to setbulk":"") );
      break;
    }
  }
  ver("correcting negative concentrations.");


  log("updating current density and conductivity...");
  updateCurrentDensity();
  writeCurrentDensity();
  ver("updating current density and conductivity.");


  log("updating electrodes current and gas production rate...");
  for (CFuint e=0; e<m_electrodes.size(); ++e)
    m_electrodes[e]->updateCurrentAndGasRate(*this);
  ver("updating electrodes current and gas production rate.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::executeOnTrs()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  const TopologicalRegionSet& trs = *getCurrentTRS();
  const ConnectivityTable< CFuint >& geo2nodes = *trs.getGeo2NodesConn();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();


  if (m_iteration>=m_linearization_switch) {
    m_linearization = "Newton";
  }


  // allocate block accumulator and nodal variables arrays
  double **coordinates    = d.allocate_double_matrix(Nvtcell,Ndim);
  double **velocities     = d.allocate_double_matrix(Nvtcell,Ndim);
  double **concentrations = d.allocate_double_matrix(Nvtcell,Nions);
  double  *potentials     = new double[Nvtcell];
  double  *temperatures   = new double[Nvtcell];
  double  *densities      = new double[Nvtcell];
  double  *voidfractions  = new double[Nvtcell];
  double **bvectors       = d.allocate_double_matrix(Nvtcell,3);
  BlockAccumulator* acc = createBlockAccumulator(Nvtcell,Nvtcell,Nsys);


  // element matrix and residual vector contribution
  RealMatrix emat(Nvtcell*Nsys,Nvtcell*Nsys);
  RealVector eres(Nvtcell*Nsys);


  log("assemble (" + m_linearization + " linearization)...");
  for (CFuint ic=0; ic<trs.getLocalNbGeoEnts(); ++ic) {

    // get/set nodes indices and nodal variables
    for (int inc=0; inc<Nvtcell; ++inc) {
      getNodalVariables( geo2nodes(ic,inc),
        coordinates[inc], velocities[inc], concentrations[inc], potentials[inc],
        temperatures[inc], densities[inc], voidfractions[inc], bvectors[inc] );
      acc->setRowColIndex(inc,geo2nodes(ic,inc));
    }

    // assemble element
    assembleElement(
      emat, eres,
      coordinates, velocities, concentrations, potentials, temperatures,
      densities, voidfractions, bvectors );

    // add to system matrix and residual vector
    acc->setValues(emat);
    matrix->addValues(*acc);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int e=0; e<Nsys; ++e)
        h_rhs(geo2nodes(ic,inc),iv+e,Neqns) += eres[inc*Nsys+e];

  }
  log("assemble.");


  // deallocate block accumulator and nodal variables arrays
  delete acc;
  delete[] bvectors[0];        delete[] bvectors;
  delete[] voidfractions;
  delete[] densities;
  delete[] temperatures;
  delete[] potentials;
  delete[] concentrations[0];  delete[] concentrations;
  delete[] velocities[0];      delete[] velocities;
  delete[] coordinates[0];     delete[] coordinates;
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::assembleElement(
  RealMatrix& emat, RealVector& eres,
  double **coordinates, double **velocities, double **concentrations,
  double *potentials, double *temperatures, double *densities,
  double *voidfractions, double **bvectors )
{
  // Ax = b

  // get element matrix contribution
  DoubleMatrix emat_ = m_assembler->calcElementMat(
    coordinates, velocities, concentrations, potentials,
    temperatures, densities, voidfractions, bvectors );

  // copy into local matrix
  for (int ir=0; ir<Nvtcell*Nsys; ++ir)
    for (int ic=0; ic<Nvtcell*Nsys; ++ic)
      emat(ir,ic) = emat_[ir][ic];

  // element residual vector contribution (-Ax part of r = b-Ax)
  for (int ir=0; ir<Nvtcell*Nsys; ++ir) {
    double sum = 0.;
    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int i=0; i<Nions; ++i)
        sum += emat_[ir][inc*Nsys+i]*concentrations[inc][i];
      sum += emat_[ir][inc*Nsys+Nions]*potentials[inc];
    }
    eres[ir] = sum;
  }

  // Picard finishes here, Newton continues further
  if (m_linearization=="Picard")
    return;

  // K = AMat + AJac - dB/dX (K*dx = r)

  // get element jacobian matrix contribution
  emat_ = m_assembler->calcElementJac(
    coordinates, velocities, concentrations, potentials,
    temperatures, densities, voidfractions, bvectors );

  // add to element contribution
  for (int ir=0; ir<Nvtcell*Nsys; ++ir)
    for (int ic=0; ic<Nvtcell*Nsys; ++ic)
      emat(ir,ic) += emat_[ir][ic];
}

//////////////////////////////////////////////////////////////////////////////

std::vector< double > SystemMITReM::getDiffusivity()
{
  cf_assert_desc(
    "system should have at least one species, not setup yet?",
    Nions>0 );
  std::vector< double > D(Nsys,0.);
  for (int i=0; i<Nions; ++i)
    D[i]= m_mitrem->getIonDiffusionConstant((unsigned) i);
  return D;
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::updateCurrentDensity()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();

  // DataHandles for current density and conductivity
  DataHandle< std::vector< double > > h_jx = s_jx.getDataHandle();
  DataHandle< std::vector< double > > h_jy = s_jy.getDataHandle();
  DataHandle< std::vector< double > > h_jz = s_jz.getDataHandle();
  DataHandle< double > h_kappa = s_conductivity.getDataHandle();

  // variables to set mitrem library state
  std::vector< double > c(Nions);
  double U = 0.;
  const double T       = m_mitrem->getSolutionTemperature();
  const double density = m_mitrem->getSolutionDensity();

  // variables to calculate Hartmann number
  const double r = (double) getMethodData().m_density;
  double kappa = 0.;  // average conductivity (size-weigthed)
  double magB  = 0.;  // average magnitude of B vector (size-weigthed)


  // allocate block accumulator and nodal variables arrays
  double **coordinates    = d.allocate_double_matrix(Nvtcell,Ndim);
  double **velocities     = d.allocate_double_matrix(Nvtcell,Ndim);
  double **concentrations = d.allocate_double_matrix(Nvtcell,Nions);
  double  *potentials     = new double[Nvtcell];
  double  *temperatures   = new double[Nvtcell];
  double  *densities      = new double[Nvtcell];
  double  *voidfractions  = new double[Nvtcell];
  double **bvectors       = d.allocate_double_matrix(Nvtcell,3);


  // calculate current density (iterate over cells)
  for (CFuint ic=0; ic<geo2nodes->nbRows(); ++ic) {

    // get/set nodal variables
    for (int inc=0; inc<Nvtcell; ++inc) {
      getNodalVariables( (*geo2nodes)(ic,inc),
        coordinates[inc], velocities[inc], concentrations[inc], potentials[inc],
        temperatures[inc], densities[inc], voidfractions[inc], bvectors[inc] );
    }

    // get element current densities
    DoubleVectorList j = m_assembler->calcIonCurrentDensities(
      coordinates, velocities, concentrations, potentials, temperatures,
      densities, voidfractions, bvectors );
    for (int i=0; i<=Nions; ++i) {
      h_jx[i][ic] = j[i][0];
      h_jy[i][ic] = j[i][1];
      if (Ndim>2)
        h_jz[i][ic] = j[i][2];
    }

  }


  // calculate conductivity (and Hartmann number) (iterate over nodes)
  for (CFuint n=0; n<h_nodes.size(); ++n) {

    // set state
    for (int i=0; i<Nions; ++i)
      c[i] = (*h_states[n])[iv+i];
    U = (*h_states[n])[iv+Nions];

    // calculate conductivity
    m_mitrem->init(&c[0],U,T,density);
    h_kappa[n] = m_mitrem->calcTransportConductivity();

    if (h_nodes[n]->isParUpdatable()) {
      kappa += h_kappa[n]*h_mn_volume[n];
      if (m_magfield_i>=0) {
        const double Bx = (*h_states[n])[m_magfield_i+0];
        const double By = (*h_states[n])[m_magfield_i+1];
        const double Bz = (*h_states[n])[m_magfield_i+2];
        magB += sqrt(Bx*Bx + By*By + Bz*Bz)*h_mn_volume[n];
      }
    }

  }

  GlobalReduceOperation< GRO_SUM >(&kappa,&kappa);
  GlobalReduceOperation< GRO_SUM >(&magB,&magB);
  kappa /= d.m_volume;
  magB  /= d.m_volume;
  const double M = magB*sqrt(kappa/(nulam*r))*(m_magfield_h/2.);
  log("Weighted conductivity: kappa = " + StringOps::to_str(kappa));
  log("Weighted magnetic field magnitude: |B| = " + StringOps::to_str(magB));
  log("Hartmann number: M = " + StringOps::to_str(M));

  // deallocate nodal variables arrays
  delete[] bvectors[0];        delete[] bvectors;
  delete[] voidfractions;
  delete[] densities;
  delete[] temperatures;
  delete[] potentials;
  delete[] concentrations[0];  delete[] concentrations;
  delete[] velocities[0];      delete[] velocities;
  delete[] coordinates[0];     delete[] coordinates;
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::writeCurrentDensity()
{
  // DataHandles
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< std::vector< double > > h_jx = s_jx.getDataHandle();
  DataHandle< std::vector< double > > h_jy = s_jy.getDataHandle();
  DataHandle< std::vector< double > > h_jz = s_jz.getDataHandle();
  DataHandle< double > h_kappa = s_conductivity.getDataHandle();
  const CFuint nbCells = geo2nodes->nbRows();

  // get filenames and file stream
  const std::string fn1 = getMethodData().getFilename("j",".plt",true);
  const std::string fn2 = getMethodData().getFilename("j","-surf.plt",true);
  std::ofstream f;
  f.precision(16);

  // write volume/species-wise current density
  f.open(fn1.c_str());
  f << "TITLE = \"current densities and conductivity\"" << std::endl;
  if (Ndim>2) {
    f << "VARIABLES = \"x0\" \"x1\" \"x2\"";
    for (unsigned i=0; i<(unsigned) Nions; ++i) {
      const std::string l = m_mitrem->getIonLabel(i);
      f << " \"Jx_" << l << "\" \"Jy_" << l << "\" \"Jz_" << l << "\"";
    }
    f << " \"Jx_total\" \"Jy_total\" \"Jz_total\""
      << " \"kappa\"" << std::endl;
    f << "ZONE T=\"\""
      << ", N=" << h_nodes.size()
      << ", E=" << nbCells
      << ", ZONETYPE=FETETRAHEDRON"
      << ", DATAPACKING=BLOCK";
    f << ", VARLOCATION=(";
    for (unsigned i=0, j=4; i<=(unsigned) Nions; ++i, j+=3)
      f << (!i? "":",")
        << j+0 << "=CELLCENTERED,"
        << j+1 << "=CELLCENTERED,"
        << j+2 << "=CELLCENTERED";
    f << ")" << std::endl;
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << (*h_nodes[n])[0];  if (!(n%500))  f << std::endl; }  f << std::endl;
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << (*h_nodes[n])[1];  if (!(n%500))  f << std::endl; }  f << std::endl;
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << (*h_nodes[n])[2];  if (!(n%500))  f << std::endl; }  f << std::endl;
    for (int i=0; i<=Nions; ++i) {
      for (CFuint c=0; c<nbCells; ++c)  {  f << " " << h_jx[i][c];      if (!(c%500))  f << std::endl; }  f << std::endl;
      for (CFuint c=0; c<nbCells; ++c)  {  f << " " << h_jy[i][c];      if (!(c%500))  f << std::endl; }  f << std::endl;
      for (CFuint c=0; c<nbCells; ++c)  {  f << " " << h_jz[i][c];      if (!(c%500))  f << std::endl; }  f << std::endl;
    }
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << h_kappa[n];        if (!(n%500))  f << std::endl; }  f << std::endl;
  }
  else {
    f << "VARIABLES = \"x0\" \"x1\"";
    for (unsigned i=0; i<(unsigned) Nions; ++i) {
      const std::string l = m_mitrem->getIonLabel(i);
      f << " \"Jx_" << l << "\" \"Jy_" << l << "\"";
    }
    f << " \"Jx_total\" \"Jy_total\""
      << " \"kappa\"" << std::endl;
    f << "ZONE T=\"\""
      << ", N=" << h_nodes.size()
      << ", E=" << nbCells
      << ", ZONETYPE=FETRIANGLE"
      << ", DATAPACKING=BLOCK";
    f << ", VARLOCATION=(";
    for (unsigned i=0, j=3; i<=(unsigned) Nions; ++i, j+=2)
      f << (!i? "":",")
        << j+0 << "=CELLCENTERED,"
        << j+1 << "=CELLCENTERED";
    f << ")" << std::endl;
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << (*h_nodes[n])[0];  if (!(n%500))  f << std::endl; }  f << std::endl;
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << (*h_nodes[n])[1];  if (!(n%500))  f << std::endl; }  f << std::endl;
    for (unsigned i=0; i<=(unsigned) Nions; ++i) {
      for (CFuint c=0; c<nbCells; ++c)  {  f << " " << h_jx[i][c];      if (!(c%500))  f << std::endl; }  f << std::endl;
      for (CFuint c=0; c<nbCells; ++c)  {  f << " " << h_jy[i][c];      if (!(c%500))  f << std::endl; }  f << std::endl;
    }
    for (CFuint n=0; n<h_nodes.size(); ++n)  {  f << " " << h_kappa[n];        if (!(n%500))  f << std::endl; }  f << std::endl;
  }
  for (CFuint c=0; c<nbCells; ++c) {
    for (int d=0; d<Nvtcell; ++d)
      f << " " << (*geo2nodes)(c,d)+1;
    f << "\n";
  }
  f.close();

  // write surface/reaction-wise current density
  f.open(fn2.c_str());
  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  const unsigned Nreac = m_mitrem->getNElecReactions();
  f << "TITLE = \"current densities and conductivity\"" << std::endl
    << "VARIABLES = \"x0\" \"x1\"" << (Ndim>2? " \"x2\"":"");
  for (unsigned r=0; r<Nreac; ++r)
    f << " \"J_" << m_mitrem->getElecReactionLabel(r) << "\"";
  f << std::endl;

  bool nvalueswritten = false;
  for (CFuint b=0; b<btrs.size(); ++b) {
    const ConnectivityTable< CFuint >& faces = *btrs[b]->getGeo2NodesConn();
    f << "ZONE T=\"" << btrs[b]->getName() << "\""
      << ", N=" << h_nodes.size()
      << ", E=" << faces.nbRows()
      << ", ZONETYPE=" << (Ndim>2? "FETRIANGLE":"FELINESEG")
      << ", DATAPACKING=BLOCK";
    if (nvalueswritten)
      f << ", VARSHARELIST=([1-" << Ndim+Nreac << "]=1)";
    f << std::endl;
    if (!nvalueswritten) {
      for (int d=0; d<Ndim; ++d) {
        for (CFuint n=0; n<h_nodes.size(); ++n) {  f << " " << (*h_nodes[n])[d];  if (!(n%500))  f << std::endl; }  f << std::endl;
      }
      for (unsigned r=0; r<Nreac; ++r) {
        for (CFuint n=0; n<h_nodes.size(); ++n) {  f << " " << m_vj[r][n];        if (!(n%500))  f << std::endl; }  f << std::endl;
      }
      nvalueswritten = true;
    }
    for (CFuint c=0; c<faces.nbRows(); ++c) {
      if (Ndim>2)
        f << faces(c,0)+1 << " " << faces(c,1)+1 << " " << faces(c,2)+1 << std::endl;
      else
        f << faces(c,0)+1 << " " << faces(c,1)+1 << std::endl;
    }
  }
  f.close();
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::setInitialSolution()
{
  log("set bulk concentrations and potential field...");
  const DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  for (CFuint n=0; n<h_states.size(); ++n) {
    for (int e=0; e<Nions; ++e)
      (*h_states[n])[iv+e] = m_bulk[e];
    (*h_states[n])[iv+Nions] = 0.;
  }
  ver("set bulk concentrations and potential field.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::getNodalVariables( CFuint n,
  double*& coordinates, double*& velocity, double*& concentrations, double& potential,
  double& temperature, double& density, double& voidfraction, double*& bvector )
{
  const DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const Node&  TheNode  = *h_nodes[n];
  const State& TheState = *h_states[n];
  const CFuint nbDim = TheNode.size();

  // from Node
  for (CFuint i=0; i<nbDim; ++i)
    coordinates[i] = TheNode[i];

  // from State
  for (CFuint i=0; i<nbDim; ++i)
    velocity[i] = (m_velocity_i>=0? TheState[m_velocity_i+i]:0.);
  for (int i=0; i<Nions; ++i)
    concentrations[i] = TheState[iv+i];
  potential = TheState[iv+Nions];
  for (CFuint i=0; i<3; ++i)
    bvector[i] = (m_magfield_i>=0? TheState[m_magfield_i+i]:0.);

  // from MITReM library
  temperature = m_mitrem->getSolutionTemperature();
  density     = m_mitrem->getSolutionDensity();

  // from DataHandle
  voidfraction = (m_voidfraction? s_mn_voidfraction.getDataHandle()[n] : 0.);
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReM::getNodalVariables( CFuint n,
  double*& coordinates, double*& concentrations, double& potential,
  double& temperature, double& density )
{
  const DataHandle< Node*,GLOBAL >  h_nodes  = s_nodes.getDataHandle();
  const DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const Node&  TheNode  = *h_nodes[n];
  const State& TheState = *h_states[n];

  // from Node
  for (CFuint i=0; i<TheNode.size(); ++i)
    coordinates[i] = TheNode[i];

  // from State
  for (int i=0; i<Nions; ++i)
    concentrations[i] = TheState[iv+i];
  potential = TheState[iv+Nions];

  // from MITReM library
  temperature = m_mitrem->getSolutionTemperature();
  density     = m_mitrem->getSolutionDensity();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

