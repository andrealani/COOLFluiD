
#include <cmath>
#include "Common/StringOps.hh"
#include "Muffin/BC.hh"
#include "Muffin/CC.hh"
#include "Muffin/System.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

System::System(const std::string& name) :
    MuffinCom(name),
    s_rhs("rhs"),                      // socket sinks
    s_nodes("nodes"),                  // ...
    s_states("states"),                // ...
    s_faceneighcell("faceNeighCell"),  // ...
    m_requires_lss(true),
    m_iteration(0)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  lssname = "";
  m_rresidual.clear();
  m_rsolution.clear();
  scaconv_str = "N-scheme";
  m_restart = false;
  setParameter("LSSName",&lssname);
  setParameter("RelaxResidual",&m_rresidual);
  setParameter("RelaxSolution",&m_rsolution);
  setParameter("ScalarScheme",&scaconv_str);
  setParameter("Restart",&m_restart);
}

//////////////////////////////////////////////////////////////////////////////

System::~System()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void System::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;

  options.addConfigOption< std::string >("LSSName","Name of linear system solver, required to use a system (default \"\")");
  options.addConfigOption< std::vector< double > >("RelaxResidual","Relaxation coefficentsRelaxation fpr right-hand side vector (default < >)");
  options.addConfigOption< std::vector< double > >("RelaxSolution","Relaxation coefficentsRelaxation fpr solution vector (default < >)");
  options.addConfigOption< std::string >("ScalarScheme","Scalar convection scheme, which can be \"none\" (no convection), \"FOU\", \"N-scheme\" (default), \"Galerkin\", \"LDA\", \"LWS\" or \"PSI\"");
  options.addConfigOption< bool >("Restart","Option to restart from the solution provided (default false)");
}

//////////////////////////////////////////////////////////////////////////////

void System::setup()
{
  CFAUTOTRACE;
  MuffinCom::setup();
  MuffinData& d = getMethodData();

  // get system index
  is   = -1;
  Nsys = 0;
  if (m_requires_lss) {
    log("set linear system...");

    MultiMethodHandle< LinearSystemSolver > LSS = d.getLinearSystemSolver();
    for (CFuint i=0; i<LSS.size(); ++i)
      if (LSS[i]->getName() == lssname) {
        is = i;
        break;
      }
    if (is<0)
      err("Didn't find LSS with name \"" + lssname + "\"");

    // get system matrix
    matrix = LSS[is]->getMatrix();

    // get equations and vector index from linear system solver equations mask
    std::valarray< bool >& mask = *(LSS[is]->getMaskArray());
    Neqns = (int) mask.size();
    iv   = -1;
    for (int e=0; e<Neqns; ++e)
      if (mask[e])
        ++Nsys;
    for (int e=0; e<Neqns; ++e)
      if (mask[e]) {
        iv = e;
        break;
      }

    log("set linear system.");
  }


  // setup relaxation vectors
  if (m_rresidual.size()) {
    if (m_rresidual.size()==1 && Nsys>0) {
      m_rresidual.assign(Nsys,m_rresidual[0]);
    }
    else if ((int) m_rresidual.size()!=Nsys) {
      err("incorrect \"RelaxResidual\" vector size");
    }
  }

  if (m_rsolution.size()) {
    if (m_rsolution.size()==1 && Nsys>0) {
      m_rsolution.assign(Nsys,m_rsolution[0]);
    }
    else if ((int) m_rsolution.size()!=Nsys) {
      err("incorrect \"RelaxSolution\" vector size");
    }
  }


  // set scalar convection scheme
  scaconv = (scaconv_str=="FOU"?      ISSFOU :
            (scaconv_str=="N-scheme"? ISSNSC :
            (scaconv_str=="Galerkin"? ISSGAL :
            (scaconv_str=="LDA"?      ISSLDA :
            (scaconv_str=="LWS"?      ISSLWS :
            (scaconv_str=="PSI"?      ISSPSI :
                                      ISSNUL ))))));


  // set initial solution
  if (!d.m_restart || !m_restart) {
    log("set initial solution...");
    setInitialSolution();
    ver("set initial solution.");
  }
}

//////////////////////////////////////////////////////////////////////////////

void System::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();


  // sum of nodes in all partitions
  CFuint LNnode = h_nodes.size();
  CFuint SNnode = 0;
  Framework::GlobalReduceOperation< GRO_SUM >(&LNnode,&SNnode);
  const CFreal Nnode = (CFreal) SNnode;
  const CFreal eps = MathTools::MathConsts::CFrealEps();


  log("reset linear system...");
  matrix->flushAssembly();
  matrix->resetToZeroEntries();
  h_rhs = 0.;
  ver("reset linear system.");


  log("assemble elements residuals...");
  NumericalCommand::execute();
  ver("assemble elements residuals.");


  std::vector< SafePtr< CC > >& cc = d.m_vcomm_cc;
  for (CFuint i=0; i<cc.size(); ++i) {
    if (cc[i]->hasTag("Muffin::System:" + getName())) {

      log("coupling condition \"" + cc[i]->getName() + "\"...");
      matrix->flushAssembly();  // required between setting and adding
      cc[i]->apply(this);
      ver("coupling condition." );

    }
  }


  std::vector< SafePtr< BC > >& bc = d.m_vcomm_bc;
  for (CFuint i=0; i<bc.size(); ++i) {
    if (bc[i]->hasTag("Muffin::System:" + getName())) {

      log("boundary condition \"" + bc[i]->getName() + "\"...");
      matrix->flushAssembly();  // required between setting and adding
      bc[i]->apply(this);
      ver("boundary condition." );

    }
  }


  log("update residual norm (of equations)...");
  for (int e=iv; e<iv+Nsys; ++e) {
    double r2 = 0.;
    for (CFuint n=0; n<h_nodes.size(); ++n)
      if (h_nodes[n]->isParUpdatable()) {
        const double r = fabs(h_rhs(n,e,Neqns));
        r2 += r*r;
      }
    GlobalReduceOperation< GRO_SUM >(&r2,&r2);
    d.m_logl2_rhs[e] = log10(sqrt(r2/Nnode)+eps);
  }
  ver("update residual norm (of equations).");


  if (m_rresidual.size()) {
    log("relax residual...");
    int e = 0;
    int i = 0;
    for (CFuint n=0; n<h_nodes.size(); ++n)
      for (e=iv, i=0; e<iv+Nsys; ++e, ++i)
        h_rhs(n,e,Neqns) *= m_rresidual[i];
    ver("relax residual.");
  }


  log("solve...");
  matrix->finalAssembly();
  d.getLinearSystemSolver()[is]->solveSys();
  ver("solve.");


  log("update solution error norm (of variables)...");
  for (int e=iv; e<iv+Nsys; ++e) {
    double r2 = 0.;
    for (CFuint n=0; n<h_nodes.size(); ++n)
      if (h_nodes[n]->isParUpdatable()) {
        const double r = fabs(h_rhs(n,e,Neqns));
        r2 += r*r;
      }
    GlobalReduceOperation< GRO_SUM >(&r2,&r2);
    d.m_logl2_states[e] = log10(sqrt(r2/Nnode)+eps);
  }
  ver("update solution error norm (of variables).");


  log("update solution...");
  update();
  ver("update solution.");


  // update internal iteration counter
  ++m_iteration;
}

//////////////////////////////////////////////////////////////////////////////

void System::update()
{
  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // apply relaxation if necessary
  int e = 0;
  int i = 0;
  if (m_rsolution.size()) {
    for (CFuint n=0; n<h_states.size(); ++n)
      for (e=iv, i=0; e<iv+Nsys; ++e, ++i)
        (*h_states[n])[e] -= m_rsolution[i] * h_rhs(n,e,Neqns);
  }
  else {
    for (CFuint n=0; n<h_states.size(); ++n)
      for (e=iv; e<iv+Nsys; ++e)
        (*h_states[n])[e] -= h_rhs(n,e,Neqns);
  }

  getMethodData().synchronise();
}

//////////////////////////////////////////////////////////////////////////////

void System::setInitialSolution()
{
  MuffinData& d = getMethodData();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  // set VectorialFunction
  VectorialFunction f;
  f.setFunctions(d.m_initialvalues_def);
  f.setVariables(d.getNodalVariables());
  try {
    f.parse();
  }
  catch (ParserException& e) {
    log("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }

  // evaluate and set State's DataHandle (only variables of this system)
  const CFuint nbVar = (h_states.size()? h_states[0]->size() : 0);
  RealVector vars(0.,d.getNodalVariables().size());
  RealVector eval(0.,nbVar);
  for (CFuint n=0; n<h_states.size(); ++n) {
    d.getNodalValues(n,vars);
    f.evaluate(vars,eval);
    for (int i=iv; i<iv+Nsys; ++i)
      (*h_states[n])[i] = eval[i];
  }

  // synchronise
  d.synchronise();
}


//////////////////////////////////////////////////////////////////////////////

BlockAccumulator* System::createBlockAccumulator(const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize)
{
  return getMethodData().getLinearSystemSolver()[is]->createBlockAccumulator(nbRows,nbCols,subBlockSize);
}

//////////////////////////////////////////////////////////////////////////////

void System::print(std::string basename)
{
  const std::string fn1 = getMethodData().getFilename(basename,"-mat");
  const std::string fn2 = getMethodData().getFilename(basename,"-rhs");
  const std::string fn3 = getMethodData().getFilename(basename,"-sol");
  std::ofstream f;

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();

  log("output matrix to \"" + fn1 + "\"...");
  matrix->finalAssembly();
  matrix->printToFile(fn1.c_str());

  log("output rhs to \"" + fn2 + "\"...");
  f.open(fn2.c_str(),std::ios::trunc );
  for (CFuint n=0; n<h_states.size(); ++n) {
    for (int e=iv; e<iv+Nsys; ++e)
      f << h_rhs(n,e,Neqns) << " ";
    f << "\n";
  }
  f.close();

  log("output solution to \"" + fn3 + "\"...");
  f.open(fn3.c_str(),std::ios::trunc );
  for (CFuint n=0; n<h_states.size(); ++n) {
    for (int e=iv; e<iv+Nsys; ++e)
      f << (*h_states[n])[e] << " ";
    f << "\n";
  }
  f.close();
}

//////////////////////////////////////////////////////////////////////////////

void System::print(RealMatrix& m)
{
  using std::cerr;
  using std::endl;
  cerr<<endl<<"RealMatrix begin:"<<endl;
  for (CFuint i=0; i<m.nbRows(); ++i) {
    for (CFuint j=0; j<m.nbCols(); ++j)
      cerr << " " << m(i,j);
    cerr << endl;
  }
  cerr<<"RealMatrix end"<<endl;
}

//////////////////////////////////////////////////////////////////////////////

void System::print(RealVector& v)
{
  using std::cerr;
  using std::endl;
  cerr<<endl<<"RealVector begin:"<<endl;
  for (CFuint i=0; i<v.size(); ++i)
    cerr << " " << v[i] << endl;
  cerr<<"RealVector end"<<endl;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

