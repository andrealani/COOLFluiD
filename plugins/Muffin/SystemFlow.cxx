
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/CC.hh"
#include "Muffin/BC.hh"
#include "Muffin/SystemFlow.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemFlow,MuffinData,MuffinModule > cSystemFlowProvider("Flow");

//////////////////////////////////////////////////////////////////////////////

SystemFlow::SystemFlow(const std::string& name) :
    System(name),
    s_mn_volume("NodalVolume")  // socket sinks
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemFlow");
  addConfigOptionsTo(this);

  m_linearization        = "Newton";
  m_linearization_switch = 666;
  m_vardensity = 207.e-6;
  m_newtoneps  = 1.e-8;
  m_vmax_s     = "-1.";
  setParameter("Linearization",&m_linearization);
  setParameter("LinearizationSwitchToNewton",&m_linearization_switch);
  setParameter("VarDensity",&m_vardensity);
  setParameter("NewtonEps",&m_newtoneps);
  setParameter("Vmax",&m_vmax_s);
}

//////////////////////////////////////////////////////////////////////////////

SystemFlow::~SystemFlow()
{
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Linearization","Jacobian matrix linearization method, which can be \"Picard\" or \"Newton\" (default)");
  options.addConfigOption< CFuint >("LinearizationSwitchToNewton","Newton linearization forced after this number of iterations (default 666)");
  options.addConfigOption< double >("VarDensity","Volumetric thermal expansion coefficient (beta, units 1/K, default = 207.e-6)");
  options.addConfigOption< double >("NewtonEps","Newton linearization perturbation (default 1.e-8)");
  options.addConfigOption< std::string >("Vmax","Maximum velocity norm function (negative means from current solution, default \"-1.\")");
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;
  System::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::setup()
{
  CFAUTOTRACE;
  System::setup();
  MuffinData& d = getMethodData();

  // variable types
  d.m_vartypes[iv] = VPRESSURE;
  for (int i=iv+1; i<iv+Nsys; ++i)
    d.m_vartypes[i] = VVELOCITY;

  // variables indices and if temperature is coupled to flow
  m_iv_temp = d.getVariableIndex(VTEMPERATURE);
  m_iv_turb = d.getVariableIndex(VTURBK);
  m_coupletemp = (m_iv_temp>=iv && m_iv_temp<iv+Nsys);

  // linearization
  if (m_linearization!="Picard" && m_linearization!="Newton")
    err("\"Linearization\" must be either \"Picard\" or \"Newton\"");

  // Vmax evaluation function setup
  std::vector< std::string > v(2);
  v[0] = "i";
  v[1] = "t";
  m_vmax_f.setVariables(v);
  m_vmax_f.setFunctions(std::vector< std::string >(1,m_vmax_s));
  try {
    m_vmax_f.parse();
  }
  catch (ParserException& e) {
    log("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }
}


//////////////////////////////////////////////////////////////////////////////

void SystemFlow::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_volume = s_mn_volume.getDataHandle();

  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size() : 0);


  log("convert pressure state to pressure/density...");
  for (CFuint n=0; n<h_states.size(); ++n)
    (*h_states[n])[iv] /= d.m_density;
  ver("convert pressure state to pressure/density.");


  if (m_iteration>=m_linearization_switch) {
    m_linearization = "Newton";
    m_rresidual.clear();
    m_rsolution.clear();
  }


  if (m_iv_temp>=0) {
    log("calculate bulk temperature...");
    m_T0 = 0.;
    for (CFuint n=0; n<h_nodes.size(); ++n)
      if (h_nodes[n]->isParUpdatable())
        m_T0 += h_mn_volume[n] * (*h_states[n])[m_iv_temp];
    Common::GlobalReduceOperation< Common::GRO_SUM >(&m_T0,&m_T0);
    m_T0 /= d.m_volume;
    log("To = " + StringOps::StringOps::to_str(m_T0));
    ver("calculate bulk temperature.");
  }


  log("find maximum absolute velocity...");
  {
    // calculate from function
    RealVector v(2), r(1);
    v[0] = d.m_status->getNbIter();
    v[1] = d.m_status->getCurrentTimeDim();
    m_vmax_f.evaluate(v,r);
    const CFreal vmax_f = std::max(r[0],PhysicalConstants::_eps);

    // calculate from solution
    CFreal vmax_s = 0.;
    for (CFuint n=0; n<h_states.size(); ++n) {
      RealSliceVector v(h_states[n]->slice(iv+1,nbDim));
      vmax_s = std::max(vmax_s,MathTools::MathFunctions::innerProd(v,v));
    }
    vmax_s = std::max(sqrt(vmax_s),PhysicalConstants::_eps);
    Common::GlobalReduceOperation< Common::GRO_MAX >(&vmax_s,&vmax_s);

    // choose one
    const bool solution = (r[0]<PhysicalConstants::_eps);
    log( "Vmax (f/s) = " + StringOps::StringOps::to_str(vmax_f) + " / " +
                           StringOps::StringOps::to_str(vmax_s) +
                           (solution? " (solution)":" (function)") );
    m_vmax = (solution? vmax_s: vmax_f);
  }
  ver("find maximum absolute velocity.");


  log("System iteration cycle...");
  System::execute();
  ver("System iteration cycle.");


  log("convert pressure/density state to pressure...");
  for (CFuint n=0; n<h_states.size(); ++n)
    (*h_states[n])[iv] *= d.m_density;
  ver("convert pressure/density state to pressure.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::executeOnTrs()
{
  log("assembling (" + m_linearization + " linearization)...");
  if (m_linearization=="Picard")
    assemblePicard(*getCurrentTRS());
  else if (m_linearization=="Newton")
    assembleNewton(*getCurrentTRS());
  ver("assembling.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::assemblePicard(const TopologicalRegionSet& trs)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();
  const ConnectivityTable< CFuint >& geo2nodes = *trs.getGeo2NodesConn();

  DataHandle< CFreal >        h_rhs    = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  struct local_node_struct No_local[4];
  double Cscalar[4][4];
  const double prandtl = d.m_prandtl;

  BlockAccumulator* acc = createBlockAccumulator(Nvtcell,Nvtcell,Nsys);

  /* allocate memory for temporary cell arrays */
  std::vector< double > sbuoy( 1+Ndim+1, 0. );
  std::vector< double > zeta(  1+Ndim,   0. );
  std::vector< double > k(     Nvtcell,  0. );
  double ***A = d.d3tensor(0, Ndim-1,    0, Ndim,   0, Ndim);
  double ***K = d.d3tensor(0, Nvtcell-1, 0, Ndim,   0, Ndim);
  double ***B = d.d3tensor(0, Nvtcell-1, 0, Ndim,   0, Ndim);
  double ***J = d.d3tensor(0, Nvtcell-1, 0, Nsys-1, 0, Nsys-1);


  /* loop over cells */
  for (CFuint ic=0; ic<trs.getLocalNbGeoEnts(); ++ic) {

    /* reset cell matrix */
    acc->reset();
    for (int inc=0; inc<Nvtcell; ++inc)
      acc->setRowColIndex(inc,geo2nodes(ic,inc));

    /* cell normals and volume */
    double vol;
    int inc_min;
    d.cellgeom(ic,No_local,&vol,&inc_min);

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = (*h_states[geo2nodes(ic,inc)])[iv];

    /* convective contribution */
    double betasq = 0.;
    double lwfac = 0.;
    addConvectiveTerms(
      No_local,vol,inc_min,
      &betasq,&lwfac,k,sbuoy,zeta,
      A,K );

    /* distribution matrix B */
    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int e1=0; e1<=Ndim; ++e1) {
        for (int e2=0; e2<=Ndim; ++e2)
          B[inc][e1][e2] = -zeta[e1]*lwfac*K[inc][e1][e2];
        B[inc][e1][e1] += -1./(double) Nvtcell;
      }
    }

    /* viscous contribution */
    const double nueff = nulam + d.getTurbulentViscosity(No_local,vol);
    addViscousTerms(No_local,vol,nueff);

    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int jnc=0; jnc<Nvtcell; ++jnc) {

        double ninj = 0.;
        for (int i=0; i<Ndim; ++i)
          ninj += No_local[inc].norm[i] * No_local[jnc].norm[i];
        ninj /= (double) Nvtfce * vol;

        for (int i=0; i<Ndim; ++i)
          acc->addValue(inc,jnc,1+i,1+i, -nueff*ninj/(double) Nvtfce);

      }
    }
    matrix->addValues( *acc );
    acc->reset();

    /* scalar equation contributions */
    if (m_coupletemp) {
      d.scacde(ic,m_iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,1);
      for (int inc=0; inc<Nvtcell; ++inc)
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          Cscalar[inc][jnc] = No_local[inc].C[jnc];
    }

    /* add residuals to right-hand side vector */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int node_inc = geo2nodes(ic,inc);
      for (int e=0; e<Nsys; ++e)
        h_rhs(node_inc,iv+e,Neqns) += No_local[inc].Res[e];
    }


    /* loop over nodes within cell */
    for (int inc=0; inc<Nvtcell; ++inc) {

      /* zero entries */
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        for (int irow=0; irow<Nsys; ++irow)
          for (int jcol=0; jcol<Nsys; ++jcol)
            J[jnc][irow][jcol] = 0.;

      /* calculate Bi*Kj entries */
      for (int jnc=0; jnc<Nvtcell; ++jnc)
        for (int irow=0; irow<=Ndim; ++irow)
          for (int jcol=0; jcol<=Ndim; ++jcol)
            for (int ij=0; ij<=Ndim; ++ij)
              J[jnc][irow][jcol] += (
                B[inc][irow][ij]*K[jnc][ij][jcol] );

      if (m_coupletemp) {

        /* add dependence of T on T (from scacde) */
        for (int jnc=0; jnc<Nvtcell; ++jnc)
          J[jnc][m_iv_temp][m_iv_temp] = (
            Cscalar[inc][jnc] );

        /* add dependence of T on V */
        if (scaconv==ISSNSC || scaconv==ISSPSI) {
          double Tin = 0.;
          double sumkminus = 0.;
          for (int jnc=0; jnc<Nvtcell; ++jnc) {
            Tin += std::min(0.,k[jnc])*No_local[jnc].W[m_iv_temp];
            sumkminus += std::min(0.,k[jnc]);
          }
          Tin /= sumkminus;
          if (k[inc]>0.)
            for (int jnc=0; jnc<Nvtcell; ++jnc)
              for (int iv=1; iv<=Ndim; ++iv)
                  J[jnc][m_iv_temp][iv] += (
                  -No_local[inc].norm[iv-1]/((double) Nvtfce*(double) Nvtcell) * (No_local[inc].W[m_iv_temp]-Tin) );
        }
        else if (scaconv==ISSLWS) {
          for (int jnc=0; jnc<Nvtcell; ++jnc)
            for (int iv=1; iv<=Ndim; ++iv)
              for (int knc=0; knc<Nvtcell; ++knc) {
                  J[jnc][m_iv_temp][iv] += (
                  -(1./(double) Nvtcell+lwfac*k[inc])*No_local[knc].norm[iv-1]*No_local[knc].W[m_iv_temp]/((double) Nvtcell*(double) Nvtfce) );
                  J[jnc][m_iv_temp][iv] += (
                  -lwfac*No_local[inc].norm[iv-1] * No_local[knc].W[m_iv_temp]*k[knc]/((double) Nvtcell*(double) Nvtfce) );
              }
        }

        /* add dependence of p and V on T (buoyancy) */
        if (buoyancy) {
          for (int jnc=0; jnc<Nvtcell; ++jnc)
            for (int iv=1; iv<=Ndim; ++iv) {
                J[jnc][0][m_iv_temp]  += (
                -(lwfac*betasq*No_local[jnc].norm[iv-1]/(double) Nvtfce)*sbuoy[iv] );
                J[jnc][iv][m_iv_temp] += (
                -(1./(double) Nvtcell+lwfac*k[jnc])*sbuoy[iv] );
            }
        }

      }

      /* add block row to block accumulator */
      for (int j=0; j<Nvtcell; ++j)
        for (int e1=0; e1<Nsys; ++e1)
          for (int e2=0; e2<Nsys; ++e2)
            acc->setValue(inc,j,e1,e2,J[j][e1][e2]);


    /* end loop over nodes within cell */
    }

    /* add contributions to matrix */
    matrix->addValues(*acc);

  /* end loop over cells */
  }

  /* free memory for temporary cell arrays */
  d.free_d3tensor(A, 0, Ndim-1,    0, Ndim,   0, Ndim);
  d.free_d3tensor(K, 0, Nvtcell-1, 0, Ndim,   0, Ndim);
  d.free_d3tensor(B, 0, Nvtcell-1, 0, Ndim,   0, Ndim);
  d.free_d3tensor(J, 0, Nvtcell-1, 0, Nsys-1, 0, Nsys-1);

  delete acc;
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::assembleNewton(const TopologicalRegionSet& trs)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();
  const ConnectivityTable< CFuint >& geo2nodes = *trs.getGeo2NodesConn();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  const CFuint nbDim = (h_nodes.size()? h_nodes[0]->size() : 0);


  int jnc[4];
  struct local_node_struct No_local[4];
  struct local_node_struct No_eps[4];
  const double prandtl = d.m_prandtl;

  BlockAccumulator* acc = createBlockAccumulator(Nvtcell,Nvtcell,Nsys);

  /* allocate memory for temporary cell arrays */
  std::vector< double > sbuoy( 1+nbDim+1, 0. );
  std::vector< double > zeta(  1+nbDim,   0. );
  std::vector< double > k(     Nvtcell,   0. );
  double ***A = d.d3tensor(0, nbDim-1,   0, nbDim,   0, nbDim);
  double ***K = d.d3tensor(0, Nvtcell-1, 0, nbDim,   0, nbDim);


  /* loop over cells  */
  for (CFuint ic=0; ic<trs.getLocalNbGeoEnts(); ++ic) {

    /* reset cell matrix */
    acc->reset();
    for (int inc=0; inc<Nvtcell; ++inc)
      acc->setRowColIndex(inc,geo2nodes(ic,inc));

    /* cell normals and volume */
    double vol;
    int inc_min;
#if 1
    d.cellgeom(ic,No_local,&vol,&inc_min);
#else
std::vector< CFuint > cnodes(nbDim,0);
geo2nodes.setRow(ic,cnodes);
d.cellgeom(cnodes);
#endif

    /* copy nodal values from global to local structure */
    for (int inc=0; inc<Nvtcell; ++inc)
      for (int iv=0; iv<Neqns; ++iv)
        No_local[inc].W[iv] = (*h_states[geo2nodes(ic,inc)])[iv];

    /* convective contribution */
    double betasq;
    double lwfac;
    addConvectiveTerms(
      No_local,vol,inc_min,
      &betasq,&lwfac,k,sbuoy,zeta,
      A,K );

    /* scalar equation contribution */
    if (m_coupletemp)
      d.scacde(ic,m_iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);

    /* viscous contribution */
    const double nueff = nulam + d.getTurbulentViscosity(No_local,vol);
    addViscousTerms(No_local,vol,nueff);

    /* add residuals to right-hand side vector */
    for (int inc=0; inc<Nvtcell; ++inc) {
      const int node_inc = geo2nodes(ic,inc);
      for (int e=0; e<Nsys; ++e)
        h_rhs(node_inc,iv+e,Neqns) += No_local[inc].Res[e];
    }

    // copy unperturbed state
    for (int inc=0; inc<Nvtcell; ++inc) {
      No_eps[inc].node  = No_local[inc].node;
      No_eps[inc].norm2 = No_local[inc].norm2;
      for (CFuint id=0; id<nbDim; ++id)
        No_eps[inc].norm[id] = No_local[inc].norm[id];
      for (int iv=0; iv<Neqns; ++iv)
        No_eps[inc].W[iv] = No_local[inc].W[iv];
    }

    const scalar_convection_type scaconv_save = scaconv;
    if (m_coupletemp && scaconv==ISSPSI) {
      scaconv = ISSNSC;
      d.scacde(ic,m_iv_temp,No_local,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);
    }


    /* loop over nodes within cell and variables */
    for (int inc=0; inc<Nvtcell; ++inc) {
      for (int iv=0; iv<Nsys; ++iv) {

        /* find node numbers of other nodes in cell */
        for (int i=0; i<Nvtcell; ++i)
          jnc[i] = (inc+i)%Nvtcell;

        /* perturb solution value */
        const double Uiv  = No_eps[inc].W[iv];
        const double pUiv = Uiv + m_newtoneps * (Uiv<0.?
           std::min(Uiv,-1.e-3) :
           std::max(Uiv, 1.e-3) );
        //const double pUiv = Uiv + m_newtoneps;
        const double iperturb = 1./(pUiv-Uiv);
        No_eps[inc].W[iv] = pUiv;

        /* convective contribution, perturbed residual */
        addConvectiveTerms(
          No_eps,vol,inc_min,
          &betasq,&lwfac,k,sbuoy,zeta,
          A,K );

        /* viscous contribution, perturbed residual */
        const double nueff = nulam + d.getTurbulentViscosity(No_local,vol);
        addViscousTerms(No_eps,vol,nueff);

        /* scalar equation contribution, perturbed residual */
        if (m_coupletemp)
          d.scacde(ic,m_iv_temp,No_eps,vol,inc_min,scaconv,nulam/prandtl,0.,0.,0);

        /* contributions */
        for (int i=0; i<Nvtcell; ++i)
          for (int jv=0; jv<Nsys; ++jv)
            acc->addValue( jnc[i],inc,jv,iv,
              ( No_eps[jnc[i]].Res[jv] - No_local[jnc[i]].Res[jv] ) * iperturb );

        /* restore unperturbed solution value */
        No_eps[inc].W[iv] = Uiv;

      }
    }
    /* end loop over nodes within cell and variables */

    /* add contributions to matrix */
    matrix->addValues(*acc);

    scaconv = scaconv_save;

  /* end loop over cells */
  }

  /* free memory for temporary cell arrays */
  d.free_d3tensor(A, 0, nbDim-1,    0, nbDim,   0, nbDim   );
  d.free_d3tensor(K, 0, Nvtcell-1,  0, nbDim,   0, nbDim   );

  delete acc;
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::addViscousTerms(struct local_node_struct *No_loc, double vol, double nueff)
{
  // gradients of velocity components
  double gradv[3][3];
  for (int i=0; i<Ndim; ++i)
    for (int j=0; j<Ndim; ++j) {
      gradv[i][j] = 0.;
      for (int inc=0; inc<Nvtcell; ++inc)
        gradv[i][j] += No_loc[inc].W[iv+1+i] * No_loc[inc].norm[j];
      gradv[i][j] /= (double) Nvtfce * vol;
    }

  // viscous tensor
  double tau[3][3];
  for (int i=0; i<Ndim; ++i)
    for (int j=0; j<Ndim; ++j)
      tau[i][j] = nueff * (gradv[i][j] + gradv[j][i]);

  // assemble and distribute viscous residual
  for (int inc=0; inc<Nvtcell; ++inc)
    for (int i=0; i<Ndim; ++i) {
      double viscres = 0.;
      for (int j=0; j<Ndim; ++j)
        viscres += tau[i][j] * No_loc[inc].norm[j];
      No_loc[inc].Res[iv+1+i] -= viscres/(double) Nvtfce;
    }
}

//////////////////////////////////////////////////////////////////////////////

void SystemFlow::addConvectiveTerms(
  struct local_node_struct *No_loc, double vol, int inc_min,
  double *Beta2, double *LWfactor, std::vector< double >& k, std::vector< double >& sbuoy, std::vector< double >& zeta,
  double ***A, double ***K )
{
  /* initialize local nodal residuals to zero */
  for (int inc=0; inc<Nvtcell; ++inc)
    for (int e=0; e<=Ndim; ++e)
      No_loc[inc].Res[e] = 0.;

  /* cell-average velocity */
  double q = 0.;
  double u[3];
  for (int id=0; id<Ndim; ++id) {
    u[id] = 0.;
    for (int inc=0; inc<Nvtcell; ++inc)
      u[id] += No_loc[inc].W[id+1];
    u[id] = u[id]/(double) Nvtcell;
    q += u[id]*u[id];
  }
  q = sqrt(q);

  const double betasq = 0.25*m_vmax*m_vmax;  // C2
  *Beta2 = betasq;

  /* calculate k[i]=a.dot.n[i]/Nvtfce (inflow coefficients) */
  double sumkabs = 0.;
  for (int inc=0; inc<Nvtcell; ++inc) {
    k[inc] = 0.;
    for (int id=0; id<Ndim; ++id)
      k[inc] += u[id]*No_loc[inc].norm[id];
    k[inc] = k[inc]/(double) Nvtfce;
    sumkabs += fabs(k[inc]);
  }


  /*
   * cell length-scale
   * h = Vc/(length or area of smallest edge)
   * h = (length or area of smallest edge)
   * if (q<1.e-10)
   *   h = 2.*vol/h;
   * else
   *   h = 2.*vol*q/sumkabs;
   */
  const double h = sqrt(No_loc[inc_min].norm2);

  /*
   * Lax-Wendroff "factor"
   * lwfac = dt/(2*vol)     : definition
   *       = h/(2*U*vol)    : SUPG (if U is magnitude of cell velocity)
   *       = 1/(2*U*A)      :
   *       = h/(m_vmax*vol) :
   */
  /* use m_vmax/2 for velocity scale
   * lwfac = dt/(2*Vc) = h/(2*U*Vc) = 1/(2*U*A)
   * lwfac = h/(m_vmax*vol);
   */
  const double lwfac = 1./(2.*m_vmax*h);
  *LWfactor = lwfac;

  /*
   * const double nueff = nulam +
   *   getMethodData().getTurbulentViscosity(No_loc,vol);
   *
   * double Reu = m_vmax*h/nueff;
   * zeta[0] = 1.;
   * zeta[0] = Reu/(1.+Reu);
   * Reu = q*h/nueff;
   * for (int iv=1; iv<=Ndim; ++iv)
   *   zeta[iv] = Reu/(1.+Reu);
   */

  /* use nuc=1 for p and nuc=1/2 for velocity (nuc: cell CFL number) */
  zeta[0] = 1.;
  for (int iv=1; iv<=Ndim; ++iv)
    zeta[iv] = 0.5;



  // A (jacobian matrices)
//const double rho = (double) getMethodData().m_density;
  for (int id=0; id<Ndim; ++id) {
    for (int irow=0; irow<=Ndim; ++irow)
      for (int jcol=0; jcol<=Ndim; ++jcol)
        A[id][irow][jcol] = 0.;
    for (int irow=1; irow<=Ndim; ++irow)
      A[id][irow][irow] = u[id];

    // additional terms for conservative form
    /*
     * for (int irow=1; irow<=Ndim; ++irow)
     *   A[id][irow][id+1] += u[irow-1];
     */
    A[id][id+1][0] = 1.;
    A[id][0][id+1] = betasq;
  }

  // Ki (inflow matrices) for each node: Ki = ni.A/Nvtfce
  for (int inc=0; inc<Nvtcell; ++inc) {
    for (int irow=0; irow<=Ndim; ++irow) {
      for (int jcol=0; jcol<=Ndim; ++jcol) {
        K[inc][irow][jcol] = 0.;
        for (int id=0; id<Ndim; ++id)
          K[inc][irow][jcol] += No_loc[inc].norm[id] * A[id][irow][jcol] / (double) Nvtfce;
      }
    }
  }

  /* cell vector residual */
  std::vector< double > Rc(1+Ndim,0.);
  for (int irow=0; irow<=Ndim; ++irow) {
    for (int inc=0; inc<Nvtcell; ++inc) {
      double KiUi = 0.;
      for (int jcol=0; jcol<=Ndim; ++jcol)
        KiUi += K[inc][irow][jcol]*No_loc[inc].W[jcol];
      Rc[irow] -= KiUi;
    }
  }

  // add -2/3*grad(k) term into momentum residuals
  /*
   * if (turmod!=ITNULL && getMethodData().getIteration()>10)
   *   for (int id=0; id<Ndim; ++id) {
   *     double gradk = 0.;
   *     for (int inc=0; inc<Nvtcell; ++inc)
   *       gradk += No_loc[inc].W[m_iv_turb+0]*No_loc[inc].norm[id];
   *     gradk /= (double) Nvtfce;
   *     Rc[id+1] -= 0.66667*gradk;
   *   }
   */

  if (buoyancy) {
    std::vector< double >& grav = getMethodData().m_gravity;

    double Tave = 0.;
    for (int inc=0; inc<Nvtcell; ++inc)
      Tave += No_loc[inc].W[m_iv_temp];
    Tave = Tave/(double) Nvtcell;
    for (int iv=1; iv<=Ndim; ++iv) {
      Rc[iv]    -= vol*grav[iv-1]*m_vardensity*(Tave-m_T0);
      sbuoy[iv]  = vol*grav[iv-1]*m_vardensity/(double) Nvtcell;
    }
    sbuoy[ 0         ] = 0.;
    sbuoy[ m_iv_temp ] = 0.;
  }


  // old-time nodal residuals for convection and pressure contribution
  // Rci = Bci * Kci * Rc
  for (int inc=0; inc<Nvtcell; ++inc) {
    for (int irow=0; irow<=Ndim; ++irow) {
      for (int jcol=0; jcol<=Ndim; ++jcol) {
        No_loc[inc].Res[irow] += zeta[irow]*lwfac*K[inc][irow][jcol]*Rc[jcol];
      }
    }
    for (int irow=0; irow<=Ndim; ++irow)
      No_loc[inc].Res[irow] += Rc[irow]/(double) Nvtcell;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

