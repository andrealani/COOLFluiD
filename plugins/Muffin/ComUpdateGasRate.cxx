
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Element.hh"
#include "Muffin/ComUpdateGasRate.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ComUpdateGasRate,MuffinData,MuffinModule > cComUpdateGasRateProvider("ComUpdateGasRate");

//////////////////////////////////////////////////////////////////////////////

ComUpdateGasRate::ComUpdateGasRate(const std::string& name) :
    MuffinCom(name),
    s_nodes("nodes"),             // socket sinks
    s_gasonsurf("GasOnSurface"),  // ...
    m_mitrem(CFNULL),
    m_issetup(false)
{
  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  // set parameters
  m_mitrem_str.clear();
  m_greactions_str.clear();
  m_issave = false;
  setParameter("MITReM",&m_mitrem_str);
  setParameter("GasReactions",&m_greactions_str);
  setParameter("Save",&m_issave);
}

//////////////////////////////////////////////////////////////////////////////

void ComUpdateGasRate::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("MITReM","MITReM system to calculate gas production rate, name (default \"\")");
  options.addConfigOption< std::vector< std::string > >("GasReactions","Vector of gas-producing reactions labels (default <>)");
  options.addConfigOption< bool >("Save","If gas production rate should be saved to file (default false)");
}

//////////////////////////////////////////////////////////////////////////////

void ComUpdateGasRate::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();

  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbDim = h_nodes[0]->size();


  // this is a delayed setup to wait for SystemMITReM to complete its setup
  if (!m_issetup) {
    m_issetup = true;

    d.log(getName() + ": set MITReM system \"" + m_mitrem_str + "\"...");
    std::vector< SafePtr< System > >& vcomm_sys = d.m_vcomm_sys;
    for (CFuint i=0; i<vcomm_sys.size() && m_mitrem.isNull(); ++i) {
      if (vcomm_sys[i]->getName()==m_mitrem_str) {
        try {
          SafePtr< SystemMITReM > p = vcomm_sys[i].d_castTo< SystemMITReM >();
          // at this point cast is successful. nothing happens when name is set
          // but differs from current one; if name is not set and cast ok, or name
          // was set and is equal to this one's name, set pointer to it and exit.
          m_mitrem = p;
        }
        catch (FailedCastException& e) {
        }
      }
    }
    if (m_mitrem.isNull())
      d.err(getName() + ": no MITReM system \"" + m_mitrem_str + "\" method found");
    d.ver(getName() + ": set MITReM system.");


    d.log(getName() + ": set gas reactions...");
    m_greactions.clear();
    std::vector< std::string > all_greactions_str;
    for (CFuint r=0; r<m_mitrem->m_mitrem->getNGasReactions(); ++r)
      all_greactions_str.push_back(m_mitrem->m_mitrem->getGasReactionLabel(r));

    m_greactions.clear();
    for (CFuint r=0; r<m_greactions_str.size(); ++r)
      for (CFuint a=0; a<all_greactions_str.size(); ++a)
        if (m_greactions_str[r]==all_greactions_str[a]) {
          m_greactions.push_back(a);
          break;
        }
    if (m_greactions.size()!=m_greactions_str.size())
      d.err(getName() + ": not all gas reactions were found");
    if (!m_greactions.size())
      d.err(getName() + ": not gas reactions were found");
    d.ver(getName() + ": set gas reactions.");
  }


  // allocate nodal variables arrays
  double **coordinates    = d.allocate_double_matrix(Nvtfce,nbDim);
  double **concentrations = d.allocate_double_matrix(Nvtfce,m_mitrem->Nions);
  double  *potentials     = new double[Nvtfce];
  double  *temperatures   = new double[Nvtfce];
  double  *densities      = new double[Nvtfce];
  RealVector surfacegasfractions(Nvtfce);


  d.log(getName() + ": update gas production rate...");
  // for all the boundary regions, cycle boundary elements
  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  for (CFuint i=0; i<btrs.size(); ++i) {
    if (!btrs[i]->hasTag("Muffin::GasProduction"))
      continue;
    const ConnectivityTable< CFuint >& faces = *btrs[i]->getGeo2NodesConn();
    double sum_dVdt = 0.;
    for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {

      // get nodes indices and copy variables into the data structures
      for (int inc=0; inc<Nvtfce; ++inc) {
        m_mitrem->getNodalVariables( faces(ifc,inc),
          coordinates[inc], concentrations[inc], potentials[inc],
          temperatures[inc], densities[inc] );
      }

      // get boundary element gas production rate
      surfacegasfractions = h_gasonsurf[i][ifc].F;
      h_gasonsurf[i][ifc].dVdt = m_mitrem->m_assembler->calcGasGeneration(
        coordinates, concentrations, potentials,
        temperatures, densities, &surfacegasfractions[0],
        m_greactions.size(), &m_greactions[0]);
      sum_dVdt += h_gasonsurf[i][ifc].dVdt;

    }
    d.log(getName() + ": gas production rate for \"" + btrs[i]->getName() + "\": " + StringOps::to_str(sum_dVdt) + " [m3 s-1]");
  }
  d.ver(getName() + ": update gas production rate.");

  // deallocate nodal variables arrays
  delete[] densities;
  delete[] temperatures;
  delete[] potentials;
  delete[] coordinates[0];     delete[] coordinates;
  delete[] concentrations[0];  delete[] concentrations;

  if (m_issave) {
    using std::endl;

    // open output file
    const std::string fn = d.getFilename("dVdt","-surf.plt",false,true);
    std::ofstream f(fn.c_str());
    f.precision(16);
    f << "TITLE = \"gas production rate\"" << endl
      << "VARIABLES = \"x0\" \"x1\"" << (nbDim>2? " \"x2\"":"")
      << " \"dVdt [m3 s-1]\" \"V [m3]\" \"F [1]\"" << endl;

    for (CFuint b=0; b<btrs.size(); ++b) {
      const ConnectivityTable< CFuint >& faces = *btrs[b]->getGeo2NodesConn();

      // zone header
      f << "ZONE T=\"" << btrs[b]->getName() << "\""
        << " N=" << h_nodes.size()
        << " E=" << faces.nbRows()
        << " ZONETYPE=" << (nbDim>2? "FETRIANGLE":"FELINESEG")
        << " DATAPACKING=BLOCK"
        << " VARLOCATION=(" << nbDim+1<< "=CELLCENTERED,"
                            << nbDim+2<< "=CELLCENTERED,"
                            << nbDim+3<< "=CELLCENTERED)";
      if (b)
        f <<  " VARSHARELIST=([1-" << nbDim << "]=1)";
      f << endl;

      // variable values: coordinates (only on first zone), current gas
      // production rate, gas volume, surface gas fraction
      if (!b) {
        for (CFuint d=0; d<nbDim; ++d)
          for (CFuint n=0; n<h_nodes.size(); ++n) f << ((n+1)%100? ' ':'\n') << (*h_nodes[n])[d];  f << endl;
      }
      if (btrs[b]->hasTag("Muffin::GasProduction")) {
        for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc)  f << ((ifc+1)%100? ' ':'\n') << h_gasonsurf[b][ifc].dVdt;  f << endl;
        for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc)  f << ((ifc+1)%100? ' ':'\n') << h_gasonsurf[b][ifc].V;     f << endl;
        for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc)  f << ((ifc+1)%100? ' ':'\n') << h_gasonsurf[b][ifc].F;     f << endl;
      }
      else {
        for (CFuint i=0; i<3; ++i) {
          for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc)  f << ((ifc+1)%100? ' ':'\n') << 0.;  f << endl;
        }
      }

      // connectivity
      for (CFuint c=0; c<faces.nbRows(); ++c) {
        for (CFuint i=0; i<faces.nbCols(c); ++i)  f << ' ' << faces(c,i)+1;  f << endl;
      }
    }

    // close output file
    f.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

