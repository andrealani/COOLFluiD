
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Element.hh"
#include "Muffin/ComSingleStates.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ComSingleStates,MuffinData,MuffinModule > cComSingleStatesProvider("ComSingleStates");

//////////////////////////////////////////////////////////////////////////////

ComSingleStates::ComSingleStates(const std::string& name) :
    MuffinCom(name),
    s_nodes("nodes"),          // socket sinks
    s_states("states"),        // ...
    s_nvolume("NodalVolume"),  // ...
    s_sol("SolutionAt(t-0)"),     // socket sources
    s_ddx("SolutionAt(t-0)ddx"),  // ...
    s_ddy("SolutionAt(t-0)ddy"),  // ...
    s_ddz("SolutionAt(t-0)ddz")   // ...
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  m_qsolution = 0;
  m_qderivatives = 0;
  setParameter("QueueSolution",&m_qsolution);
  setParameter("QueueDerivatives",&m_qderivatives);
}

//////////////////////////////////////////////////////////////////////////////

void ComSingleStates::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("QueueSolution","If solution is queued (default 0)");
  options.addConfigOption< CFuint >("QueueDerivatives","If solution spacial derivatives are queued (default 0)");
}

//////////////////////////////////////////////////////////////////////////////

void ComSingleStates::setup()
{
  CFAUTOTRACE;
  MuffinCom::setup();
  MuffinData& d = getMethodData();

  // set default Node and State sizes
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
  nbNodes = h_nodes.size();
  nbStates = h_states.size();
  nbDim = h_nodes[0]->size();
  nbVar = h_states[0]->size();
  

  if (m_qsolution) {
    d.log("reset solution...");
    DataHandle< CFreal > h_sol = s_sol.getDataHandle();
    h_sol.resize(nbStates*nbVar);
    h_sol = 0.;
    d.ver("reset solution.");
  }


  if (m_qderivatives) {
    d.log("reset solution derivatives components...");
    DataHandle< CFreal > h_ddx = s_ddx.getDataHandle();
    DataHandle< CFreal > h_ddy = s_ddy.getDataHandle();
    DataHandle< CFreal > h_ddz = s_ddz.getDataHandle();
    switch (nbDim) {
      case 3:  h_ddz.resize(nbStates*nbVar);  h_ddz = 0.;
      case 2:  h_ddy.resize(nbStates*nbVar);  h_ddy = 0.;
      case 1:  h_ddx.resize(nbStates*nbVar);  h_ddx = 0.;
      default: break;
    }
    d.ver("reset solution derivatives components.");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComSingleStates::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();


  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();


  if (m_qsolution) {
    d.log("solution update...");
    DataHandle< CFreal > h_sol = s_sol.getDataHandle();
    for (CFuint i=0; i<nbStates; ++i)
      for (CFuint j=0; j<nbVar; ++j)
        h_sol(i,j,nbVar) = (*h_states[i])[j];
    d.ver("solution update.");
  }


  if (m_qderivatives) {
    DataHandle< CFreal > h_ddx = s_ddx.getDataHandle();
    DataHandle< CFreal > h_ddy = s_ddy.getDataHandle();
    DataHandle< CFreal > h_ddz = s_ddz.getDataHandle();


    // auxiliary variables
    std::vector< CFreal > vdx(nbVar,0.),  // derivatives, x-component
                          vdy(nbVar,0.),  // ... y-component
                          vdz(nbVar,0.);  // ... z-component
    AElement* e(nbDim==2? (AElement*) new ElementTriangle(2) :
                          (AElement*) new ElementTetrahedron(3) );
    std::vector< Node* > vnodes(e->N);
    std::vector< State* > vstates(e->S);


    // reset derivatives
    switch (nbDim) {
      case 3:  h_ddz = 0.;
      case 2:  h_ddy = 0.;
      case 1:  h_ddx = 0.;
      default: break;
    }


    d.log("solution derivatives update...");
    std::vector< SafePtr< TopologicalRegionSet > > vtrs =
      MeshDataStack::getActive()->getFilteredTrsList("inner");
    for (CFuint t=0; t<vtrs.size(); ++t) {  // cycle all 'inner' TRS's
      const ConnectivityTable< CFuint >& elems = *vtrs[t]->getGeo2NodesConn();

      // calculate states derivatives, and distribute weighing with nodal size
      // (per element)
      for (CFuint c=0; c<elems.nbRows(); ++c) {

        for (CFuint j=0; j<e->N; ++j)  vnodes[j]  = h_nodes[elems(c,j)];
        for (CFuint j=0; j<e->N; ++j)  vstates[j] = h_states[elems(c,j)];
        e->element(vnodes);
        switch (nbDim) {
          case 3:  vdz = e->dd(ZZ,vstates);
          case 2:  vdy = e->dd(YY,vstates);
          case 1:  vdx = e->dd(XX,vstates);
          default: break;
        }
        const CFreal w = e->s/CFreal(e->N);
        for (CFuint j=0; j<nbVar; ++j) {
          for (CFuint i=0; i<e->N; ++i) {
            switch (nbDim) {
              case 3:  h_ddz( elems(c,i) ,j,nbVar) += vdz[j]*w;
              case 2:  h_ddy( elems(c,i) ,j,nbVar) += vdy[j]*w;
              case 1:  h_ddx( elems(c,i) ,j,nbVar) += vdx[j]*w;
              default: break;
            }
          }
        }

      }

    }  // cycle all 'inner' TRS's
    d.ver("solution derivatives update.");


    // de-weight the derivatives (per node)
    DataHandle< CFreal > h_nvolume = s_nvolume.getDataHandle();
    for (CFuint n=0; n<h_nvolume.size(); ++n) {
      const CFreal w = h_nvolume[n];
      if (w>PhysicalConstants::_eps) {
        for (CFuint j=0; j<nbVar; ++j) {
          switch (nbDim) {
            case 3:  h_ddz(n,j,nbVar) /= w;
            case 2:  h_ddy(n,j,nbVar) /= w;
            case 1:  h_ddx(n,j,nbVar) /= w;
            default: break;
          }
        }
      }
    }


    // deallocate AElement
    delete e;
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

