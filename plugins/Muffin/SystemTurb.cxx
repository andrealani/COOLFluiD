
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/BCFile.hh"
#include "Muffin/BCFixVelocity.hh"
#include "Muffin/BCFunction.hh"
#include "Muffin/SystemFlow.hh"
#include "Muffin/SystemTurb.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemTurb,MuffinData,MuffinModule > cSystemTurbProvider("Turb");

//////////////////////////////////////////////////////////////////////////////

SystemTurb::SystemTurb(const std::string& name) :
    System(name),
    s_mn_volume("NodalVolume"),              // socket sinks
    s_mn_bnormal("NodalNormals"),            // ...
    s_mn_walldistance("NodalWallDistance"),  // ...
    s_mn_wallnode("NodalWallNode")           // ...
{
  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemTurb");
  addConfigOptionsTo(this);

  m_turmod_str.clear();
  m_referencevalues.clear();
  m_iteration_init = 5;
  setParameter("Model",&m_turmod_str);
  setParameter("ReferenceValues",&m_referencevalues);
  setParameter("IterationStart",&m_iteration_init);
}

//////////////////////////////////////////////////////////////////////////////

SystemTurb::~SystemTurb()
{
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Model","Tubulence model identifier (default \"\")");
  options.addConfigOption< std::vector< double > >("ReferenceValues","k-e/w reference values vector (default < 0. 0. >)");
  options.addConfigOption< CFuint >("IterationStart","Start real coefficients calculation after this number of iterations (default 5)");
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  System::configure(args);

  // reference values
  if (!m_referencevalues.size()) {
//FIXME this option has to be set, otherwise epsilon at the walls is zero (and some other bad things happen.)
//TODO substitute this by querying MuffinData initial values
//    m_referencevalues.assign(2,0.);
//  else if (m_referencevalues.size()!=2)
    err("\"ReferenceValues\" size must be 2");
  }
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::setup()
{
  CFAUTOTRACE;
  System::setup();
  MuffinData& d = getMethodData();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();

  // variable types
  d.m_vartypes[iv+0] = VTURBK;
  d.m_vartypes[iv+1] = (turmod_ke? VTURBE:VTURBW);


  log("setting turbulence model identifier...");
  m_turmod = ITNULL;
  if (     m_turmod_str=="KE2L") { m_turmod = ITKE2L; }
  else if (m_turmod_str=="KELB") { m_turmod = ITKELB; }
  else if (m_turmod_str=="KENA") { m_turmod = ITKENA; }
  else if (m_turmod_str=="KWHR") { m_turmod = ITKWHR; }
  else if (m_turmod_str=="KWLR") { m_turmod = ITKWLR; }
  else if (m_turmod_str=="KWPD") { m_turmod = ITKWPD; }
  else if (m_turmod_str=="SST")  { m_turmod = ITKWSS; }
  else if (m_turmod_str=="BSL")  { m_turmod = ITKWBS; }
  else
    err("\"Model\" option \"" + m_turmod_str + "\" not recognized");
  log("setting turbulence model identifier.");


  log("setting turbulence model constants...");
  if (m_turmod==ITKE2L) {
    // two-layer model
    Cmu   = 0.09;
    Ceps1 = 1.44;
    Ceps2 = 1.92;
    Sig1  = 1.0;
    Sig2  = 1.3;
  }
  else if (m_turmod==ITKELB) {
    // Lam-Bremhorst model
    Cmu   = 0.09;
    Ceps1 = 1.44;
    Ceps2 = 1.92;
    Sig1  = 1.0;
    Sig2  = 1.3;
    Bmu   = 0.0165;
    Dmu   = 20.5;
    Aeps1 = 0.05;
  }
  else if (m_turmod==ITKENA) {
    // Abe-Kondoh-Nagano model
    Cmu   = 0.09;
    Ceps1 = 1.5;
    Ceps2 = 1.9;
    Sig1  = 1.4;
    Sig2  = 1.4;
  }
  else if (m_turmod==ITKWHR || m_turmod==ITKWLR) {
    // Wilcox k-w model
    Ck   = 0.09;
    Cw1  = 0.5555;
    Cw2  = 0.075;
    Sig1 = 2.0;
    Sig2 = 2.0;
  }
  else if (m_turmod==ITKWBS) {
    // k-w BSL model
    Betas  = 0.09;
    Beta1  = 0.075;
    Sigk1  = 0.5;
    Sigw1  = 0.5;
    Gamma1 = Beta1/Betas-0.168*Sigw1/sqrt(Betas);
    Beta2  = 0.0828;
    Sigk2  = 1.0;
    Sigw2  = 0.856;
    Gamma2 = Beta2/Betas-0.168*Sigw2/sqrt(Betas);
    Cw2    = 0.075;
  }
  else if (m_turmod==ITKWSS) {
    // k-w SST model
    Betas  = 0.09;
    Beta1  = 0.075;
    Sigk1  = 0.85;
    Sigw1  = 0.5;
    Gamma1 = Beta1/Betas-0.168*Sigw1/sqrt(Betas);
    Beta2  = 0.0828;
    Sigk2  = 1.0;
    Sigw2  = 0.856;
    Gamma2 = Beta2/Betas-0.168*Sigw2/sqrt(Betas);
    Cw2    = 0.075;
    A1     = 0.31;
  }
  else if (m_turmod==ITKWPD) {
    // Peng-Davidson-Holmberg k-w model
    Ck   = 0.09;
    Cw1  = 0.42;
    Cw2  = 0.075;
    Cw   = 0.75;
    Sig1 = 0.8;
    Sig2 = 1.35;
  }
  log("setting turbulence model constants.");


  log("reset turbulent viscosity...");
  m_nuturb.assign(h_nodes.size(),0.);
  log("reset turbulent viscosity.");


  log("set turbulence length...");
  m_lenturb.resize(h_nodes.size());

  const double eref = m_referencevalues[1];
  for (CFuint n=0; n<h_nodes.size(); ++n) {
    const double len = pow(1.-h_mn_walldistance[n]/eref,2.);
    m_lenturb[n] = 0.53*eref*(0.14-len*(0.08+0.06*len));
  }

  for (CFuint i=0; i<d.m_walls.size(); ++i) {
    std::vector< CFuint >& wallnodes = (*(d.m_walls)[i]->getNodesInTrs());
    for (CFuint inb=0; inb<wallnodes.size(); ++inb)
      m_lenturb[wallnodes[inb]] = 1.e-20;
  }
  log("set turbulence length.");


  if (!d.m_restart) {
    log("set initial k-e/w field...");

    // get velocity norm
    const double v = d.getReferenceVelocityNorm();

    // k: intensity, e: width/length
    const double kref = m_referencevalues[0];
    const double eref = m_referencevalues[1];
    const double k = 1.5 * kref*kref * v*v;
    const double e = (turmod_ke?
      pow(Cmu,0.75)*pow(k,1.5) / (0.09*eref) :
      sqrt(k) / (0.09*eref) );

    DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
    for (CFuint n=0; n<h_states.size(); ++n) {
      (*h_states[n])[iv+0] = k;
      (*h_states[n])[iv+1] = e;
    }

    log("set initial k-e/w field.");
  }
}


//////////////////////////////////////////////////////////////////////////////

void SystemTurb::executeOnTrs()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();

  int inc_min;
  double vol;
  double nuturb;
  struct local_node_struct No_local[4];

  BlockAccumulator* acc = createBlockAccumulator(Nvtcell,Nvtcell,2);

  /* loop over cells  */
  for (CFuint ic=0; ic<geo2nodes->nbRows(); ++ic) {

    /* reset cell matrix */
    acc->reset();
    for (int inc=0; inc<Nvtcell; ++inc)
      acc->setRowColIndex(inc,(*geo2nodes)(ic,inc));

    /* cell normals and volume */
    d.cellgeom(ic,No_local,&vol,&inc_min);

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = (*h_states[(*geo2nodes)(ic,inc)])[iv];

    /* calculate turbulent viscosity */
    nuturb = turb_viscosity(No_local,vol);

    /* convection and diffusion of turbulent variables */
    if (m_turmod==ITKWBS || m_turmod==ITKWSS) {
      double turb1 = 0.;
      double turb2 = 0.;
      double wdist = 0.;
      double *gradk = new double[3];
      double *gradw = new double[3];
      for (int inc=0; inc<Nvtcell; ++inc) {
        turb1 += No_local[inc].W[iv+0];
        turb2 += No_local[inc].W[iv+1];
        wdist += h_mn_walldistance[No_local[inc].node];
      }
      turb1 /= (double) Nvtcell;
      turb2 /= (double) Nvtcell;
      wdist /= (double) Nvtcell;

      for (int id=0; id<Ndim; ++id) {
        gradk[id] = 0.;
        gradw[id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc) {
          gradk[id] += (*h_states[No_local[inc].node])[iv+0] * No_local[inc].norm[id];
          gradw[id] += (*h_states[No_local[inc].node])[iv+1] * No_local[inc].norm[id];
        }
        gradk[id] /= (double) Nvtfce*vol;
        gradw[id] /= (double) Nvtfce*vol;
      }

      double dkdw=0.;
      for (int id=0; id<Ndim; ++id)
        dkdw += gradk[id]*gradw[id];

      const double F1 = (m_iteration<=m_iteration_init? 1.:F1_function(turb1,turb2,nulam,wdist,dkdw) );

      Sig1 = 1./(F1*Sigk1 + (1.-F1)*Sigk2);
      Sig2 = 1./(F1*Sigw1 + (1.-F1)*Sigw2);

      delete[] gradk;
      delete[] gradw;
    }

    d.scacde(ic,iv+0,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig1,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        acc->setValue(inc,jnc,0,0, No_local[inc].C[jnc]);

    d.scacde(ic,iv+1,No_local,vol,inc_min,scaconv,nulam+nuturb/Sig2,0.,0.,1);
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        acc->setValue(inc,jnc,1,1, No_local[inc].C[jnc]);

    /* add contributions to matrix */
    matrix->addValues(*acc);

    /* add contributions to vector */
    for (int inc=0; inc<Nvtcell; ++inc) {
      h_rhs((*geo2nodes)(ic,inc),iv+0,Neqns) += No_local[inc].Res[iv+0];
      h_rhs((*geo2nodes)(ic,inc),iv+1,Neqns) += No_local[inc].Res[iv+1];
    }

  /* end loop over cells  */
  }
  delete acc;


  log("nodal sources...");
  turb_source_node();
  ver("nodal sources.");


  log("near-wall nodes statistics...");
  ypStatistics();
  ver("near-wall nodes statistics.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::update()
{
  ver("updating solution...");
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // set relaxation vector
  std::vector< double > linrlx(2,1.);
  if (m_rsolution.size())
    linrlx = m_rsolution;

  // negative values relaxation coefficient: relaxation of old solution value:
  // 1.00 use old solution (original), 0.00 set to zero (miotras 0.20)
  const double aneg = 0.20;

  // non-linear update of turbulence variables, correcting negative solution
  // and residual vector (not to confuse norms)
  int nkneg = 0;
  int neneg = 0;
  for (CFuint n=0; n<h_states.size(); ++n) {
    double k = (*h_states[n])[iv+0] - linrlx[0] * h_rhs(n,iv+0,Neqns);
    double e = (*h_states[n])[iv+1] - linrlx[1] * h_rhs(n,iv+1,Neqns);
    if (k<PhysicalConstants::_eps) {
      ++nkneg;
      k = (*h_states[n])[iv+0]*aneg;
      h_rhs(n,iv+0,Neqns) = PhysicalConstants::_eps;
    }
    if (e<PhysicalConstants::_eps) {
      ++neneg;
      e = (*h_states[n])[iv+1]*aneg;
      h_rhs(n,iv+1,Neqns) = PhysicalConstants::_eps;
    }
    (*h_states[n])[iv+0] = k;
    (*h_states[n])[iv+1] = e;
  }
  getMethodData().synchronise();

  log("negative k-e nodes: " + StringOps::to_str(nkneg) +
    " " + StringOps::to_str(neneg));

  ver("updating solution.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::ypStatistics()
{
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< RealVector > h_mn_bnormal = s_mn_bnormal.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();
  DataHandle< CFuint > h_mn_wallnode = s_mn_wallnode.getDataHandle();

  CFuint Nwnode = 0;
  double ypmax = 0.;     double utmax = 0.;  int inu_ypmax = -1;
  double ypmin = 1.e20;  double utmin = 0.;  int inu_ypmin = -1;
  int yp_nnear = 0;
  int yp_nfar  = 0;

  // for all the walls regions, cycle the nodes
  const std::vector< SafePtr< TopologicalRegionSet > >& vwalls = getMethodData().m_walls;
  for (CFuint i=0; i<vwalls.size(); ++i) {
    const std::vector< CFuint >& nodes = *vwalls[i]->getNodesInTrs();
    Nwnode += nodes.size();
    for (CFuint inb=0; inb<nodes.size(); ++inb) {

      const int inu      = nodes[inb];
      const int intw     = h_mn_wallnode[inu];
      const double wdist = h_mn_walldistance[inu];
      const RealVector& normal = h_mn_bnormal[inu];

      double vdotn = 0.;
      for (int id=0; id<Ndim; ++id)
        vdotn += (*h_states[intw])[id+1] * normal[id];

      double vt = 0.;
      for (int id=0; id<Ndim; ++id)
        vt += pow( (*h_states[intw])[id+1] - vdotn*normal[id] - (*h_states[inu])[id+1] ,2.);
      vt = sqrt(vt);

      // tau and y+
      const double utau  = sqrt(nulam*vt/wdist);
      const double yp = wdist*utau/nulam;

      // statistics
      if (yp>ypmax) {
        utmax     = utau;
        ypmax     = yp;
        inu_ypmax = intw;
      }
      if (yp<ypmin && yp>1.e-8) {
        utmin     = utau;
        ypmin     = yp;
        inu_ypmin = intw;
      }
      if (yp>2.)
        ++yp_nfar;
      else if (yp<1.)
        ++yp_nnear;

    }
  }

  log("  y+ max = " + StringOps::to_str(ypmax) + "  node = " + StringOps::to_str(inu_ypmax) + "  utau = " + StringOps::to_str(utmax));
  log("  y+ min = " + StringOps::to_str(ypmin) + "  node = " + StringOps::to_str(inu_ypmin) + "  utau = " + StringOps::to_str(utmin));
  log("  y+ > 2.: " + StringOps::to_str(yp_nfar)  + " nodes (" + StringOps::to_str(100.*(double)yp_nfar /(double) Nwnode) + "%)");
  log("  y+ < 1.: " + StringOps::to_str(yp_nnear) + " nodes (" + StringOps::to_str(100.*(double)yp_nnear/(double) Nwnode) + "%)");
}

//////////////////////////////////////////////////////////////////////////////

double SystemTurb::turb_viscosity(struct local_node_struct *No_loc, double vol)
{
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();

  double fmu = 1.;
  double len = 0.;
  double nuturb = 0.;

  // cell-averaged wall distance and reference turbulence length
  double dwallc = 0.;
  if (turmod_walldistance) {
    for (int inc=0; inc<Nvtcell; ++inc) {
      dwallc += h_mn_walldistance[No_loc[inc].node];
      len    += m_lenturb[No_loc[inc].node];
    }
    dwallc = dwallc/(double) Nvtcell;
    len    = len/(double) Nvtcell;
  }

  // cell-averaged turbulence intensity and length
  double kturb = 0.;
  double eturb = 0.;
  for (int inc=0; inc<Nvtcell; ++inc) {
    kturb += No_loc[inc].W[iv+0];
    eturb += No_loc[inc].W[iv+1];
  }
  kturb = std::max(1.e-10,kturb)/(double) Nvtcell;
  eturb = std::max(1.e-10,eturb)/(double) Nvtcell;

  // cell-averaged turbulent viscosity
  if (m_turmod==ITKE2L) {
    /* two-layer model */
    const double Ry = sqrt(kturb)*dwallc/nulam;
    nuturb = Cmu*kturb*kturb/eturb;
    if (nuturb/nulam<20.0)
      nuturb = pow(Cmu,0.25)*sqrt(kturb)*0.41*dwallc*(1.-exp(-Ry/70.0));
  }
  else if (m_turmod==ITKELB) {
    /* Lam-Bremhorst */
    const double Rt = kturb*kturb/(nulam*eturb);
    const double Ry = sqrt(kturb)*dwallc/nulam;
    fmu = pow(1.-exp(-Bmu*Ry),2.)*(1.+Dmu/(Rt+PhysicalConstants::_eps));
    nuturb = fmu*Cmu*kturb*kturb/eturb;
  }
  else if (m_turmod==ITKENA) {
    /* Abe-Kondoh-Nagano */
    const double ystar = pow(nulam*eturb,0.25)*dwallc/nulam;
    const double Rt = kturb*kturb/(nulam*eturb);
    fmu = pow(1.-exp(-ystar/14.0),2.) * (1.+5.*pow(Rt,-0.75)*exp(-Rt*Rt/4.0e4));
    nuturb = Cmu*fmu*kturb*kturb/eturb;
  }
  else if (m_turmod==ITKWHR || m_turmod==ITKWBS) {
    /* Standard model and BSL model*/
    nuturb = kturb/eturb;
  }
  else if (m_turmod==ITKWSS) {
    /* k-w SST model */
    double gradv[3][3];
    for (int iv=1; iv<=Ndim; ++iv)
      for (int id=0; id<Ndim; ++id) {
        gradv[iv-1][id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc)
          gradv[iv-1][id] += No_loc[inc].W[iv]*No_loc[inc].norm[id];
        gradv[iv-1][id] = gradv[iv-1][id]/((double) Nvtfce*vol);
      }
    double Omega = (gradv[1][0]-gradv[0][1])*(gradv[1][0]-gradv[0][1]);
    if (Ndim==3) {
      Omega += (gradv[2][1]-gradv[1][2])*(gradv[2][1]-gradv[1][2]);
      Omega += (gradv[0][2]-gradv[2][0])*(gradv[0][2]-gradv[2][0]);
    }
    Omega = sqrt(Omega);

    const double F2 = tanh( pow(std::max(2.0*sqrt(kturb)/(0.09*eturb*dwallc),500.0*nulam/(dwallc*dwallc*eturb)),2.) );
    nuturb = A1*kturb/std::max(A1*eturb,Omega*F2);
  }
  else if (m_turmod==ITKWLR) {
    /* Low-Re k-w model */
    const double Rt  = kturb/(nulam*eturb);
    fmu = (0.15+Rt)/(6.0+Rt);
    nuturb = fmu*kturb/eturb;
  }
  else if (m_turmod==ITKWPD) {
    /* Peng-Davidson-Holmberg k-w model */
    const double Rt = kturb/(nulam*eturb);
    fmu = 0.025 + (0.975+(0.001/Rt)*exp(-Rt*Rt/4.e4)) * ( 1.-exp(-pow(0.1*Rt,0.75)) );
    nuturb = fmu*kturb/eturb;
  }

  if (m_iteration<=m_iteration_init)
    nuturb = (turmod_ke?
      1.125*fmu*sqrt(kturb)*len :
      1000.*nulam );

  return(nuturb);
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::turb_source_P(
  double k,
  double turb2,
  double nu_t,
  double nu_l,
  double wd,
  double G,
  double gradkw,
  double *source_k,
  double *source_e )
{
  const double Pk = nu_t*G;

  if (m_turmod==ITKELB) {
    const double Rt = k*k/(nu_l*turb2);
    const double Ry = sqrt(k)*wd/nu_l;
    const double fmu = pow(1.-exp(-Bmu*Ry),2.)*(1.+Dmu/(Rt+PhysicalConstants::_eps));
    const double feps1 = (m_iteration<=m_iteration_init?
      1. :
      1.+pow(Aeps1/(fmu+PhysicalConstants::_eps),3.) );
    *source_k = Pk;
    *source_e = Ceps1*feps1*turb2*Pk/k;
  }
  else if (m_turmod==ITKENA) {
    *source_k = Pk;
    *source_e = Ceps1*turb2*Pk/k;
  }
  else if (m_turmod==ITKWHR) {
    *source_k = Pk;
    *source_e = Cw1*turb2*Pk/k;
  }
  else if (m_turmod==ITKWLR) {
    const double Rt = k/(nu_l*turb2);
    const double fmu = (0.15+Rt)/(6.0+Rt);
    const double fw = (0.27+Rt)/((2.7+Rt)*fmu);
    *source_k = Pk;
    *source_e = fw*Cw1*turb2*Pk/k;
  }
  else if (m_turmod==ITKWBS) {
    const double F1 = (m_iteration<=m_iteration_init?
      1. :
      F1_function(k,turb2,nu_l,wd,gradkw) );
    const double Gamma = F1*Gamma1 + (1.-F1)*Gamma2;
    *source_k = Pk;
    *source_e = (Gamma*Pk/nu_t) + (2.*Sigw2*(1.-F1)*gradkw/turb2);
  }
  else if (m_turmod==ITKWSS) {
    const double F1    = F1_function(k,turb2,nu_l,wd,gradkw);
    const double Gamma = F1*Gamma1 + (1.-F1)*Gamma2;
    *source_k = Pk;
    *source_e = (Gamma*Pk/nu_t) + (2.*Sigw2*(1.-F1)*gradkw/turb2);
  }
  else if (m_turmod==ITKWPD) {
    const double Rt = k/(nu_l*turb2);
    const double fw = 1. + 4.3*exp(-sqrt(0.6667*Rt));
    *source_k = Pk;
    *source_e = fw*Cw1*turb2*Pk/k + Cw*nu_t*gradkw/k;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::turb_source_D(
  double k,
  double turb2,
  double nu_l,
  double wd,
  double len,
  double gradkw,
  double *source_k,
  double *source_e,
  double *deriv_k,
  double *deriv_e,
  double *deriv_ke,
  double *deriv_ek )
{
  if (m_turmod==ITKELB) {
    const double Rt = k*k/(nu_l*turb2);
    const double feps2 = (m_iteration<=m_iteration_init?
      1. :
      1.-exp(-Rt*Rt));

    *source_k = -turb2;
    *source_e = -Ceps2*feps2*turb2*turb2/k;
    *deriv_k  = -2.*turb2/k;
    *deriv_e  = -2.*Ceps2*feps2*turb2/k;
    *deriv_ke = -1.;

    /* Use k-l model for early iterations */
    if (m_iteration<=m_iteration_init) {
      *source_k = -0.08*pow(k,1.5)/len;
      *deriv_k  = -0.12*sqrt(k)/len;
    }
  }
  else if (m_turmod==ITKENA) {
    const double ystar = pow(nu_l*turb2,0.25)*wd/nu_l;
    const double Rt    = k*k/(nu_l*turb2);
    const double feps2 = (m_iteration<=m_iteration_init?
      1. :
      pow(1.-exp(-ystar/3.1),2.)*( 1.-0.3*exp(-Rt*Rt/42.25) ) );

    *source_k = -turb2;
    *source_e = -Ceps2*feps2*turb2*turb2/k;
    *deriv_k  = -2.*turb2/k;
    *deriv_e  = -2.*Ceps2*feps2*turb2/k;
    *deriv_ke = -1.;
    *deriv_ek =  Ceps2*feps2*turb2*turb2/(k*k);

    /* Use k-l model for early iterations */
    if (m_iteration<=m_iteration_init) {
      *source_k = -0.08*pow(k,1.5)/len;
      *deriv_k  = -0.12*sqrt(k)/len;
    }
  }
  else if (m_turmod==ITKWHR) {
    *source_k = -Ck*k*turb2;
    *source_e = -Cw2*turb2*turb2;
    *deriv_k  = -Ck*turb2;
    *deriv_e  = -2.*Cw2*turb2;
    *deriv_ke = -Ck*k;
  }
  else if (m_turmod==ITKWLR) {
    const double Rt = k/(nu_l*turb2);
    const double Rt4 = Rt*Rt*Rt*Rt;
    const double fk = (1137.8+Rt4)/(4096.0+Rt4);

    *source_k = -fk*Ck*k*turb2;
    *source_e = -Cw2*turb2*turb2;
    *deriv_k  = -fk*Ck*turb2;
    *deriv_e  = -2.*Cw2*turb2;
    *deriv_ke = -fk*Ck*k;
  }
  else if (m_turmod==ITKWBS) {
    const double F1 = (m_iteration<=m_iteration_init?
      1. :
      F1_function(k,turb2,nu_l,wd,gradkw) );
    const double Beta = F1*Beta1 + (1.-F1)*Beta2;

    *source_k = -Betas*k*turb2;
    *source_e = -Beta*turb2*turb2;
    *deriv_k  = -Betas*turb2;
    *deriv_e  = -2.*Beta*turb2;
    *deriv_ke = -Betas*k;
  }
  else if (m_turmod==ITKWSS) {
    const double F1 = F1_function(k,turb2,nu_l,wd,gradkw);
    const double Beta = F1*Beta1 + (1.-F1)*Beta2;

    *source_k = -Betas*k*turb2;
    *source_e = -Beta*turb2*turb2;
    *deriv_k  = -Betas*turb2;
    *deriv_e  = -2.*Beta*turb2;
    *deriv_ke = -Betas*k;
  }
  else if (m_turmod==ITKWPD) {
    const double Rt = k/(nu_l*turb2);
    const double fk = 1. - 0.722*exp(-1.e-4*Rt*Rt*Rt*Rt);

    *source_k = -fk*Ck*k*turb2;
    *source_e = -Cw2*turb2*turb2;
    *deriv_k  = -fk*Ck*turb2;
    *deriv_e  = -2.*Cw2*turb2;
    *deriv_ke = -fk*Ck*k;
  }
}

//////////////////////////////////////////////////////////////////////////////

void SystemTurb::turb_source_node()
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();
  DataHandle< CFreal > h_mn_walldistance = s_mn_walldistance.getDataHandle();
  const CFuint nbNodes = s_nodes.getDataHandle().size();

  double gradv[4][3];
  double gradk[3];
  double gradw[3];
  struct local_node_struct No_loc[4];

  double dkdw     = 0.;
  double *No_dkdw = CFNULL;
  double dwall    = 0.;
  double lenturb  = 0.;

  if (m_turmod==ITKWSS || m_turmod==ITKWBS || m_turmod==ITKWPD) {
    No_dkdw = new double[nbNodes];
    for (CFuint n=0; n<nbNodes; ++n)
      No_dkdw[n] = 0.;
  }

  for (CFuint inu=0; inu<nbNodes; ++inu)
    m_nuturb[inu] = 0.;

  /* nodal G function values */
  for (CFuint ic=0; ic<geo2nodes->nbRows(); ++ic) {

    double vol;
    int inc_min;
    getMethodData().cellgeom(ic,No_loc,&vol,&inc_min);

    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_loc[inc].W[iv] = (*h_states[No_loc[inc].node])[iv];

    /* turbulent viscosity */
    const double nuturb = turb_viscosity(No_loc,vol);

    for (int inc=0; inc<Nvtcell; ++inc)
      m_nuturb[No_loc[inc].node] += nuturb*vol/(double) Nvtcell;

    /* average turbulence values */
    double kturb = 0.;
    double turb2 = 0.;
    for (int inc=0; inc<Nvtcell; ++inc) {
      kturb += No_loc[inc].W[iv+0];
      turb2 += No_loc[inc].W[iv+1];
    }
    kturb /= (double) Nvtcell;
    turb2 /= (double) Nvtcell;

    /* velocity gradient and G function */
    for (int v=1; v<=Ndim; ++v)
      for (int id=0; id<Ndim; ++id) {
        gradv[v][id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc)
          gradv[v][id] += No_loc[inc].W[v] * No_loc[inc].norm[id];
        gradv[v][id] = gradv[iv][id]/((double) Nvtfce*vol);
      }
    const double Gfunc = Gfunction(gradv);

    /* average wall distance for cell */
    if (turmod_walldistance) {
      dwall = 0.;
      for (int inc=0; inc<Nvtcell; ++inc)
        dwall += h_mn_walldistance[No_loc[inc].node];

      dwall /= (double) Nvtcell;

      /* for early iterations calculate average turbulence length scale */
      if (m_iteration<=m_iteration_init) {
        lenturb = 0.;
        for (int inc=0; inc<Nvtcell; ++inc)
          lenturb += m_lenturb[No_loc[inc].node];
        lenturb /= (double) Nvtcell;
      }
    }


    /* grad(k).grad(w) for BSL, SST and PDH models */
    if (m_turmod==ITKWSS || m_turmod==ITKWBS || m_turmod==ITKWPD) {
      for (int id=0; id<Ndim; ++id) {
        gradk[id] = 0.;
        gradw[id] = 0.;
        for (int inc=0; inc<Nvtcell; ++inc) {
         gradk[id] += (*h_states[No_loc[inc].node])[iv+0] * No_loc[inc].norm[id];
         gradw[id] += (*h_states[No_loc[inc].node])[iv+1] * No_loc[inc].norm[id];
        }
        gradk[id] /= (double) Nvtfce * vol;
        gradw[id] /= (double) Nvtfce * vol;
      }

      dkdw = 0.;
      for (int id=0; id<Ndim; ++id)
        dkdw += gradk[id] * gradw[id];

      for (int inc=0; inc<Nvtcell; ++inc)
        No_dkdw[No_loc[inc].node] += dkdw * vol / (double) Nvtcell;
    }


    // source terms for production
    double SP1;
    double SP2;
    turb_source_P(kturb,turb2,nuturb,nulam,dwall,Gfunc,dkdw,&SP1,&SP2);

    const double source1 = SP1 * vol / (double) Nvtcell;
    const double source2 = SP2 * vol / (double) Nvtcell;

    /* add cell-based sources and jacobian entries */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int inu = No_loc[inc].node;
      h_rhs(inu,iv+0,Neqns) += source1;
      h_rhs(inu,iv+1,Neqns) += source2;
    }

  }


  for (CFuint inu=0; inu<nbNodes; ++inu)
    m_nuturb[inu] /= h_mn_volume[inu];


  /* add nodal sources and jacobian entries */
  for (CFuint inu=0; inu<nbNodes; ++inu) {

    if (turmod_walldistance) {
      dwall = h_mn_walldistance[inu];
      if (m_iteration<=m_iteration_init)
        lenturb = m_lenturb[inu];
    }

    dkdw = 0.;
    if (m_turmod==ITKWSS || m_turmod==ITKWBS || m_turmod==ITKWPD)
      dkdw = No_dkdw[inu] / h_mn_volume[inu];

    // source terms for dissipation (-SD1)
    double SD1  = 0.;
    double SD2  = 0.;
    double JD1  = 0.;
    double JD2  = 0.;
    double JD12 = 0.;
    double JD21 = 0.;
    turb_source_D(
      (*h_states[inu])[iv+0],
      (*h_states[inu])[iv+1],
      nulam,dwall,lenturb,dkdw,
      &SD1, &SD2,
      &JD1, &JD2, &JD12, &JD21 );

    h_rhs(inu,iv+0,Neqns) += SD1 * h_mn_volume[inu];
    h_rhs(inu,iv+1,Neqns) += SD2 * h_mn_volume[inu];

    matrix->addValue(inu*Nsys+0, inu*Nsys+0, h_mn_volume[inu]*JD1);
    matrix->addValue(inu*Nsys+1, inu*Nsys+1, h_mn_volume[inu]*JD2);
    matrix->addValue(inu*Nsys+0, inu*Nsys+1, h_mn_volume[inu]*JD12);
//    matrix->addValue(inu*Nsys+1, inu*Nsys+0, h_mn_volume[inu]*JD21);  // i added this
//std::cout << "  JD1  = " << JD1  << std::endl;
//std::cout << "  JD2  = " << JD2  << std::endl;
//std::cout << "  JD12 = " << JD12 << std::endl;
//std::cout << "  JD21 = " << JD21 << std::endl;
//exit(0);

  }

  if (m_turmod==ITKWSS || m_turmod==ITKWBS || m_turmod==ITKWPD)
    delete[] No_dkdw;
}

//////////////////////////////////////////////////////////////////////////////

double SystemTurb::Gfunction(double gradv[4][3])
{
  double G = 0.;
  for (int iv=1; iv<=Ndim; ++iv)
    G += 2.*gradv[iv][iv-1]*gradv[iv][iv-1];
  for (int i=1; i<Ndim; ++i)
    for (int j=i+1; j<=Ndim; ++j)
      G += (gradv[j][i-1]+gradv[i][j-1])*(gradv[j][i-1]+gradv[i][j-1]);
  return(G);
}

//////////////////////////////////////////////////////////////////////////////

double SystemTurb::F1_function(double k, double w, double nu, double y, double dkdw)
{
  double arg1;
  const double CDkw = std::max(2.*Sigw2*dkdw/w,1.e-20);
  const double y2 = y*y;
  arg1 = std::max(sqrt(k)/(0.09*w*y),500.0*nu/(y2*w));
  arg1 = std::min(arg1,4.0*Sigw2*k/(CDkw*y2));
  return(tanh(arg1*arg1*arg1*arg1));
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

