
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemMITReM.hh"
#include "Muffin/BCElectrode.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCElectrode,MuffinData,MuffinModule > cBCElectrodeProvider("Electrode");

//////////////////////////////////////////////////////////////////////////////

BCElectrode::BCElectrode(const std::string& name) :
    BC(name),
    s_gasonsurf("GasOnSurface")  // socket sinks
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCElectrode");
  addConfigOptionsTo(this);

  // set parameters
  m_potential_str = "0.";
  m_ereactions_str.clear();
  m_greactions_str.clear();
  setParameter("Potential",&m_potential_str);
  setParameter("Reactions",&m_ereactions_str);
  setParameter("GasReactions",&m_greactions_str);
}

//////////////////////////////////////////////////////////////////////////////

BCElectrode::~BCElectrode()
{
}

//////////////////////////////////////////////////////////////////////////////

void BCElectrode::defineConfigOptions(Config::OptionList& options)
{
  CFAUTOTRACE;
  options.addConfigOption< std::string >("Potential","Metal potential (default \"0.\")");
  options.addConfigOption< std::vector< std::string > >("Reactions","Vector of electrode reactions labels (from *.elecreactions file) (default <>)");
  options.addConfigOption< std::vector< std::string > >("GasReactions","Vector of gas-producing reactions labels (from *.elecreactions file) (default <>)");
}

//////////////////////////////////////////////////////////////////////////////

void BCElectrode::setup()
{
  CFAUTOTRACE;
  BC::setup();

  // VectorialFunction setup (potential)
  m_potential_f.setVariables(getMethodData().getNodalVariables());
  m_potential_f.setFunctions(std::vector< std::string >(1,m_potential_str));
  try {
    m_potential_f.parse();
  }
  catch (ParserException& e) {
    log("VectorialFunction parsing: " + std::string(e.what()));
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCElectrode::setupReactions(
  SystemMITReM& s,
  const std::vector< std::string >& ereactions,
  const std::vector< std::string >& greactions )
{
  CFAUTOTRACE;

  // setup electrode reactions
  m_ereactions.clear();
  for (CFuint r=0; r<m_ereactions_str.size(); ++r)
    for (CFuint a=0; a<ereactions.size(); ++a)
      if (m_ereactions_str[r]==ereactions[a]) {
        m_ereactions.push_back(a);
        break;
      }
  if (m_ereactions.size()!=m_ereactions_str.size())
    err("not all electrode reactions were found");

  // setup gas reactions
  m_greactions.clear();
  for (CFuint r=0; r<m_greactions_str.size(); ++r)
    for (CFuint a=0; a<greactions.size(); ++a)
      if (m_greactions_str[r]==greactions[a]) {
        m_greactions.push_back(a);
        break;
      }
  if (m_greactions.size()!=m_greactions_str.size())
    err("not all gas reactions were found");

  // mark gas producing regions
  if (m_greactions.size()) {
    std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
    for (CFuint i=0; i<trs.size(); ++i)
      trs[i]->attachTag("Muffin::GasProduction");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCElectrode::applyOnSystemMITReM(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< CFreal > h_rhs = s_rhs.getDataHandle();
  DataHandle< CFuint > h_mn_priority = s_mn_priority.getDataHandle();
  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();

  SystemMITReM& sys = dynamic_cast< SystemMITReM& >(*s);

  // some values to simplify the task
  const int iv   = sys.iv;
  const int Nsys = sys.Nsys;
  const unsigned* p_ereactions = (m_ereactions.size()? &m_ereactions[0]:CFNULL);
  const unsigned* p_greactions = (m_greactions.size()? &m_greactions[0]:CFNULL);

  RealVector surfacegasfractions(Nvtfce);

  // allocate block accumulator and nodal variables arrays
  double **coordinates    = d.allocate_double_matrix(Nvtfce,Ndim);
  double **concentrations = d.allocate_double_matrix(Nvtfce,sys.Nions);
  double  *potentials     = new double[Nvtfce];
  double  *temperatures   = new double[Nvtfce];
  double  *densities      = new double[Nvtfce];
  BlockAccumulator* acc = sys.createBlockAccumulator(Nvtfce,Nvtfce,Nsys);


  // metal potential
  double V;
  double Vavg = 0.;
  RealVector Vvars(0.,d.getNodalVariables().size());
  RealVector Veval(0.,1);


  // cycle boundary elements
  const ConnectivityTable< CFuint >& faces = *t->getGeo2NodesConn();
  const CFuint itrs = getBoundaryTrsID(*t);
  for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {

    // get/set nodes indices and nodal variables
    V = 0.;
    for (int inc=0; inc<Nvtfce; ++inc) {
      sys.getNodalVariables( faces(ifc,inc),
        coordinates[inc], concentrations[inc], potentials[inc],
        temperatures[inc], densities[inc] );
      acc->setRowColIndex(inc,faces(ifc,inc));
      d.getNodalValues(faces(ifc,inc),Vvars);
      m_potential_f.evaluate(Vvars,Veval);
      V += Veval[0]/(double) Nvtfce;
    }
    Vavg += V;
    surfacegasfractions = (m_greactions.size()? h_gasonsurf[itrs][ifc].F : 0.);

    // get face matrix and residual vector contributions
    DoubleMatrix fmat = sys.m_assembler->calcBoundaryElementJac(
      coordinates, concentrations, potentials,
      temperatures, densities, &surfacegasfractions[0],
      m_ereactions.size(), p_ereactions, V,
      m_greactions.size(), p_greactions );
    DoubleVector fres = sys.m_assembler->calcBoundaryElementVec(
      coordinates, concentrations, potentials,
      temperatures, densities, &surfacegasfractions[0],
      m_ereactions.size(), p_ereactions, V,
      m_greactions.size(), p_greactions );

    // add to system matrix, respecting node priority
    for (int iR=0; iR<Nvtfce; ++iR) {
      const bool dontset = h_mn_priority[faces(ifc,iR)]>m_bnpriority;
      for (int ir=0; ir<Nsys; ++ir)
        for (int iC=0; iC<Nvtfce; ++iC)
          for (int ic=0; ic<Nsys; ++ic)
            acc->setValue( iR, iC, ir, ic,
              dontset? 0. : fmat[iR*Nsys+ir][iC*Nsys+ic] );
    }
    sys.matrix->addValues(*acc);

    // add to residual vector (b part of r = b-Ax), respecting node priority
    for (int inc=0; inc<Nvtfce; ++inc)
      if (h_mn_priority[faces(ifc,inc)]<=m_bnpriority)
        for (int e=0; e<Nsys; ++e)
          h_rhs(faces(ifc,inc),iv+e,Neqns) -= fres[inc*Nsys+e];

  }
  Vavg /= (double) std::max((CFuint) 1,faces.nbRows());
  log("TRS \"" + t->getName() + "\" average metal potential:" + StringOps::to_str(Vavg));


  // deallocate block accumulator and nodal variables arrays
  delete acc;
  delete[] densities;
  delete[] temperatures;
  delete[] potentials;
  delete[] concentrations[0];  delete[] concentrations;
  delete[] coordinates[0];     delete[] coordinates;
}

//////////////////////////////////////////////////////////////////////////////

void BCElectrode::updateCurrentAndGasRate(SystemMITReM& sys)
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();

  // some values to simplify the task
  const unsigned* p_ereactions = (m_ereactions.size()? &m_ereactions[0]:CFNULL);

  RealVector surfacegasfractions(Nvtfce);

  // allocate nodal variables arrays
  double **coordinates    = d.allocate_double_matrix(Nvtfce,Ndim);
  double **concentrations = d.allocate_double_matrix(Nvtfce,sys.Nions);
  double  *potentials     = new double[Nvtfce];
  double  *temperatures   = new double[Nvtfce];
  double  *densities      = new double[Nvtfce];


  // metal potential
  double V;
  RealVector Vvars(0.,d.getNodalVariables().size());
  RealVector Veval(0.,1);


  // for all the boundary regions, cycle boundary elements
  std::vector< SafePtr< TopologicalRegionSet > >& trs = getTrsList();
  for (CFuint i=0; i<trs.size(); ++i) {
    const ConnectivityTable< CFuint >& faces = *trs[i]->getGeo2NodesConn();
    const CFuint itrs = getBoundaryTrsID(*trs[i]);
    sys.m_vi[itrs] = 0.;
    for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {

      // get nodes indices and copy variables into the data structures
      V = 0.;
      for (int inc=0; inc<Nvtfce; ++inc) {
        sys.getNodalVariables( faces(ifc,inc),
          coordinates[inc], concentrations[inc], potentials[inc],
          temperatures[inc], densities[inc] );
        d.getNodalValues(faces(ifc,inc),Vvars);
        m_potential_f.evaluate(Vvars,Veval);
        V += Veval[0]/(double) Nvtfce;
      }
      surfacegasfractions = (m_greactions.size()? h_gasonsurf[itrs][ifc].F : 0.);

      // get boundary element current density and current
      DoubleListList j = sys.m_assembler->calcElecReactionCurrentDensities(
        coordinates, concentrations, potentials,
        temperatures, densities, &surfacegasfractions[0],
        m_ereactions.size(), p_ereactions, V );
      for (unsigned r=0; r<sys.m_vj.size(); ++r)
        for (int inc=0; inc<Nvtfce; ++inc)
          sys.m_vj[r][faces(ifc,inc)] = j[inc][r];
      sys.m_vi[itrs] += sys.m_assembler->calcCurrent(
        coordinates, concentrations, potentials,
        temperatures, densities, &surfacegasfractions[0],
        m_ereactions.size(), p_ereactions, V );

      // get boundary element gas production rate
      if (m_greactions.size()) {
        h_gasonsurf[itrs][ifc].dVdt = sys.m_assembler->calcGasGeneration(
          coordinates, concentrations, potentials,
          temperatures, densities, &surfacegasfractions[0],
          m_greactions.size(), &m_greactions[0]);
      }

    }
    GlobalReduceOperation< GRO_SUM >(&sys.m_vi[itrs],&sys.m_vi[itrs]);
    log(trs[i]->getName()+" current: " + StringOps::to_str(sys.m_vi[itrs]) + " A");
  }

  // deallocate nodal variables arrays
  delete[] densities;
  delete[] temperatures;
  delete[] potentials;
  delete[] concentrations[0];  delete[] concentrations;
  delete[] coordinates[0];     delete[] coordinates;
}

//////////////////////////////////////////////////////////////////////////////

CFuint BCElectrode::getBoundaryTrsID(const TopologicalRegionSet& trs)
{
  const std::vector< SafePtr< TopologicalRegionSet > >& btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  for (CFuint i=0; i<btrs.size(); ++i) {
    if (btrs[i]->getName()==trs.getName())
      return i;
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

