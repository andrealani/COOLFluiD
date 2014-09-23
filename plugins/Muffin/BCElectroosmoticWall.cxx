
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Element.hh"
#include "Muffin/System.hh"
#include "Muffin/BCElectroosmoticWall.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< BCElectroosmoticWall,MuffinData,MuffinModule > cBCElectroosmoticWallProvider("ElectroosmoticWall");

//////////////////////////////////////////////////////////////////////////////

BCElectroosmoticWall::BCElectroosmoticWall(const std::string& name) :
    BCWall(name),
    s_mn_bnarea("NodalArea"),
    m_potential_i(-1)
{
  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::BCElectroosmoticWall");
  addConfigOptionsTo(this);

  // set parameters
  m_zeta = 0.;
  m_eps  = 0.;
  setParameter("ZetaPotential",&m_zeta);
  setParameter("MediumPermitivity",&m_eps);
}

//////////////////////////////////////////////////////////////////////////////

void BCElectroosmoticWall::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("ZetaPotential","Zeta potential (default 0.)");
  options.addConfigOption< CFreal >("MediumPermitivity","Medium permitivity (default 0.)");
}

//////////////////////////////////////////////////////////////////////////////

void BCElectroosmoticWall::setup()
{
  CFAUTOTRACE;
  BCWall::setup();

  log("get potential field...");
  const std::vector< var_type >& vartypes = getMethodData().m_vartypes;
  for (CFuint i=0; i<vartypes.size(); ++i)
    if (vartypes[i]==VPOTENTIAL) {
      m_potential_i = i;
      break;
    }
  if (m_potential_i<0)
    err("potential field not found");
  ver("get potential field.");
}

//////////////////////////////////////////////////////////////////////////////

void BCElectroosmoticWall::applyOnSystemFlow(const SafePtr< System > s, const SafePtr< TopologicalRegionSet > t)
{
  CFAUTOTRACE;

  // DataHandles for nodes, states and face neighboring cells
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  DataHandle< CFreal > h_mn_bnarea = s_mn_bnarea.getDataHandle();

  // dynamic viscosity
  const CFreal mu = getMethodData().m_kviscosity * getMethodData().m_density;


  // create AElement
  const CFuint nbDim = h_nodes[0]->size();
  AElement* e(nbDim==2? (AElement*) new ElementLineseg(2) :
                        (AElement*) new ElementTriangle(3) );
  std::vector< Node* > vnodes(e->N);    // 'inner' element nodes
  std::vector< CFreal > vstates(e->S);  // 'inner' element state potential
  RealVector fnormal(nbDim);  // 'face' element normal
  CFreal fsize;               // 'face' element size


  ver("imposing slip-velocities...");
  std::vector< RealVector > vslip(h_nodes.size(),RealVector(0.,nbDim));

  // distribute potential gradient to nodes
  const ConnectivityTable< CFuint >& faces = *t->getGeo2NodesConn();
  for (CFuint ifc=0; ifc<faces.nbRows(); ++ifc) {

    // face properties
    for (CFuint j=0; j<e->N; ++j)  vnodes[j] = h_nodes[faces(ifc,j)];
    for (CFuint j=0; j<e->S; ++j)  vstates[j] = (*h_states[faces(ifc,j)])[m_potential_i];
    fsize = e->element(vnodes);

    // potential gradient in boundary element ("tangent")
    RealVector gradU(nbDim);
    gradU[XX] = e->dd(XX,vstates);
    gradU[YY] = e->dd(YY,vstates);
    if (nbDim>2)
      gradU[ZZ] = e->dd(ZZ,vstates);

    // distribute potential gradient to nodes (nodal area-weighted)
    for (CFuint j=0; j<e->N; ++j) {
      const CFuint n = faces(ifc,j);
      const CFreal w = (fsize/e->N) / h_mn_bnarea[n];
      vslip[n] -= gradU*m_eps*m_zeta/mu * w;  // minus for E=-gradU
    }
  }


  // deallocate AElement
  delete e;


  // apply to nodes
  const std::vector< CFuint >& nodes = *t->getNodesInTrs();
  for (std::vector< CFuint >::const_iterator n=nodes.begin(); n!=nodes.end(); ++n) {
    for (CFuint d=0; d<nbDim; ++d)
      setDirichletCondition( *s,*n,d+1,
        (*h_states[*n])[s->iv+1+d] - vslip[*n][d] );
  }
  ver("imposing slip-velocities.");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

