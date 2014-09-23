
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Element.hh"
#include "Muffin/ComMultipleStates.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ComMultipleStates,MuffinData,MuffinModule > cComMultipleStatesProvider("ComMultipleStates");

//////////////////////////////////////////////////////////////////////////////

ComMultipleStates::ComMultipleStates(const std::string& name) :
    MuffinCom(name),
    s_nodes("nodes"),          // socket sinks
    s_states("states"),        // ...
    s_nvolume("NodalVolume"),  // ...
    s_tm0_sol("SolutionAt(t-0)"),     // socket sources
    s_tm0_ddx("SolutionAt(t-0)ddx"),  // ...
    s_tm0_ddy("SolutionAt(t-0)ddy"),  // ...
    s_tm0_ddz("SolutionAt(t-0)ddz"),  // ...
    s_tm1_sol("SolutionAt(t-1)"),     // ...
    s_tm1_ddx("SolutionAt(t-1)ddx"),  // ...
    s_tm1_ddy("SolutionAt(t-1)ddy"),  // ...
    s_tm1_ddz("SolutionAt(t-1)ddz"),  // ...
    s_tm2_sol("SolutionAt(t-2)"),     // ...
    s_tm2_ddx("SolutionAt(t-2)ddx"),  // ...
    s_tm2_ddy("SolutionAt(t-2)ddy"),  // ...
    s_tm2_ddz("SolutionAt(t-2)ddz"),  // ...
    m_firstexecute(true)
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  addConfigOptionsTo(this);

  m_qsolution = 3;
  m_qderivatives = 3;
  setParameter("QueueSolution",&m_qsolution);
  setParameter("QueueDerivatives",&m_qderivatives);
}

//////////////////////////////////////////////////////////////////////////////

void ComMultipleStates::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("QueueSolution","Number of solutions to queue (0 to 3, default 3)");
  options.addConfigOption< CFuint >("QueueDerivatives","Number of solution spacial derivatives to queue (0 to 3, default 3)");
}

//////////////////////////////////////////////////////////////////////////////

void ComMultipleStates::setup()
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
    DataHandle< CFreal > h_tm0_sol = s_tm0_sol.getDataHandle();
    DataHandle< CFreal > h_tm1_sol = s_tm1_sol.getDataHandle();
    DataHandle< CFreal > h_tm2_sol = s_tm2_sol.getDataHandle();
    switch (m_qsolution) {
      case 3:  h_tm2_sol.resize(nbStates*nbVar);  h_tm2_sol = 0.;
      case 2:  h_tm1_sol.resize(nbStates*nbVar);  h_tm1_sol = 0.;
      case 1:  h_tm0_sol.resize(nbStates*nbVar);  h_tm0_sol = 0.;
      default: break;
    }
    d.ver("reset solution.");
  }


  if (nbDim>0 && m_qderivatives) {
    d.log("reset solution derivatives x-components...");
    DataHandle< CFreal > h_tm0_ddx = s_tm0_ddx.getDataHandle();
    DataHandle< CFreal > h_tm1_ddx = s_tm1_ddx.getDataHandle();
    DataHandle< CFreal > h_tm2_ddx = s_tm2_ddx.getDataHandle();
    switch (m_qderivatives) {
      case 3:  h_tm2_ddx.resize(nbStates*nbVar);  h_tm2_ddx = 0.;
      case 2:  h_tm1_ddx.resize(nbStates*nbVar);  h_tm1_ddx = 0.;
      case 1:  h_tm0_ddx.resize(nbStates*nbVar);  h_tm0_ddx = 0.;
      default: break;
    }
    d.ver("reset solution derivatives x-components.");
  }


  if (nbDim>1 && m_qderivatives) {
    d.log("reset solution derivatives y-components...");
    DataHandle< CFreal > h_tm0_ddy = s_tm0_ddy.getDataHandle();
    DataHandle< CFreal > h_tm1_ddy = s_tm1_ddy.getDataHandle();
    DataHandle< CFreal > h_tm2_ddy = s_tm2_ddy.getDataHandle();
    switch (m_qderivatives) {
      case 3:  h_tm2_ddy.resize(nbStates*nbVar);  h_tm2_ddy = 0.;
      case 2:  h_tm1_ddy.resize(nbStates*nbVar);  h_tm1_ddy = 0.;
      case 1:  h_tm0_ddy.resize(nbStates*nbVar);  h_tm0_ddy = 0.;
      default: break;
    }
    d.ver("reset solution derivatives y-components.");
  }


  if (nbDim>2 && m_qderivatives) {
    d.log("reset solution derivatives z-components...");
    DataHandle< CFreal > h_tm0_ddz = s_tm0_ddz.getDataHandle();
    DataHandle< CFreal > h_tm1_ddz = s_tm1_ddz.getDataHandle();
    DataHandle< CFreal > h_tm2_ddz = s_tm2_ddz.getDataHandle();
    switch (m_qderivatives) {
      case 3:  h_tm2_ddz.resize(nbStates*nbVar);  h_tm2_ddz = 0.;
      case 2:  h_tm1_ddz.resize(nbStates*nbVar);  h_tm1_ddz = 0.;
      case 1:  h_tm0_ddz.resize(nbStates*nbVar);  h_tm0_ddz = 0.;
      default: break;
    }
    d.ver("reset solution derivatives z-components.");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComMultipleStates::execute()
{
  CFAUTOTRACE;
  MuffinData& d = getMethodData();


  // solution queues
  if (m_qsolution) {
    DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
    DataHandle< CFreal > h_tm0_sol = s_tm0_sol.getDataHandle();
    DataHandle< CFreal > h_tm1_sol = s_tm1_sol.getDataHandle();
    DataHandle< CFreal > h_tm2_sol = s_tm2_sol.getDataHandle();
    switch (m_qsolution) {
      case 3: {
        d.log("solution update queue 3...");
        for (CFuint i=0; i<nbStates*nbVar; ++i)
          h_tm2_sol[i] = h_tm1_sol[i];
        d.ver("solution update queue 3.");
      }
      case 2: {
        d.log("solution update queue 2...");
        for (CFuint i=0; i<nbStates*nbVar; ++i)
          h_tm1_sol[i] = h_tm0_sol[i];
        d.ver("solution update queue 2.");
      }
      case 1: {
        d.log("solution update queue 1...");
        // new input comes from this point
        for (CFuint i=0; i<nbStates; ++i)
          for (CFuint j=0; j<nbVar; ++j)
            h_tm0_sol(i,j,nbVar) = (*h_states[i])[j];
        d.ver("solution update queue 1.");
      }
      default: break;
    }
  }
  

  // derivatives queues
  if (m_qderivatives) {
    DataHandle< CFreal > h_tm0_ddx = s_tm0_ddx.getDataHandle();
    DataHandle< CFreal > h_tm0_ddy = s_tm0_ddy.getDataHandle();
    DataHandle< CFreal > h_tm0_ddz = s_tm0_ddz.getDataHandle();
    DataHandle< CFreal > h_tm1_ddx = s_tm1_ddx.getDataHandle();
    DataHandle< CFreal > h_tm1_ddy = s_tm1_ddy.getDataHandle();
    DataHandle< CFreal > h_tm1_ddz = s_tm1_ddz.getDataHandle();
    DataHandle< CFreal > h_tm2_ddx = s_tm2_ddx.getDataHandle();
    DataHandle< CFreal > h_tm2_ddy = s_tm2_ddy.getDataHandle();
    DataHandle< CFreal > h_tm2_ddz = s_tm2_ddz.getDataHandle();
    switch (m_qsolution) {
      case 3: {
        d.log("solution derivatives update queue 3...");
        switch (nbDim) {
          case 3:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddz[i] = h_tm1_ddz[i];
          case 2:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddy[i] = h_tm1_ddy[i];
          default:  for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddx[i] = h_tm1_ddx[i];
        }
        d.ver("solution derivatives update queue 3.");
      }
      case 2: {
        d.log("solution derivatives update queue 2...");
        switch (nbDim) {
          case 3:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddz[i] = h_tm0_ddz[i];
          case 2:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddy[i] = h_tm0_ddy[i];
          default:  for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddx[i] = h_tm0_ddx[i];
        }
        d.ver("solution derivatives update queue 2.");
      }
      case 1: {
        d.log("solution derivatives update queue 1...");
        // new input comes from this point
        std::vector< SafePtr< TopologicalRegionSet > > vtrs =
          MeshDataStack::getActive()->getFilteredTrsList("inner");
        for (CFuint i=0; i<vtrs.size(); ++i)
          setStatesDerivatives(*vtrs[i],h_tm0_ddx,h_tm0_ddy,h_tm0_ddz);
        d.ver("solution derivatives update queue 1.");
      }
      default: break;
    }
  }


  // on first execution, copy unused queues from first solution queue
  if (m_firstexecute && m_qsolution) {
    DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();
    DataHandle< CFreal > h_tm0_sol = s_tm0_sol.getDataHandle();
    DataHandle< CFreal > h_tm1_sol = s_tm1_sol.getDataHandle();
    DataHandle< CFreal > h_tm2_sol = s_tm2_sol.getDataHandle();
    switch (m_qsolution) {
      case 3: {
        d.log("solution update queue 3 from 1...");
        for (CFuint i=0; i<nbStates*nbVar; ++i)
          h_tm2_sol[i] = h_tm0_sol[i];
        d.ver("solution update queue 3 from 1.");
      }
      case 2: {
        d.log("solution update queue 2 from 1...");
        for (CFuint i=0; i<nbStates*nbVar; ++i)
          h_tm1_sol[i] = h_tm0_sol[i];
        d.ver("solution update queue 2 from 1.");
      }
      default: break;
    }
  }


  // on first execution, copy unused queues from first derivatives queues
  if (m_firstexecute && m_qderivatives) {
    DataHandle< CFreal > h_tm0_ddx = s_tm0_ddx.getDataHandle();
    DataHandle< CFreal > h_tm0_ddy = s_tm0_ddy.getDataHandle();
    DataHandle< CFreal > h_tm0_ddz = s_tm0_ddz.getDataHandle();
    DataHandle< CFreal > h_tm1_ddx = s_tm1_ddx.getDataHandle();
    DataHandle< CFreal > h_tm1_ddy = s_tm1_ddy.getDataHandle();
    DataHandle< CFreal > h_tm1_ddz = s_tm1_ddz.getDataHandle();
    DataHandle< CFreal > h_tm2_ddx = s_tm2_ddx.getDataHandle();
    DataHandle< CFreal > h_tm2_ddy = s_tm2_ddy.getDataHandle();
    DataHandle< CFreal > h_tm2_ddz = s_tm2_ddz.getDataHandle();
    switch (m_qsolution) {
      case 3: {
        d.log("solution derivatives update queue 3 from 1...");
        switch (nbDim) {
          case 3:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddz[i] = h_tm0_ddz[i];
          case 2:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddy[i] = h_tm0_ddy[i];
          default:  for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm2_ddx[i] = h_tm0_ddx[i];
        }
        d.ver("solution derivatives update queue 3 from 1.");
      }
      case 2: {
        d.log("solution derivatives update queue 2 from 1...");
        switch (nbDim) {
          case 3:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddz[i] = h_tm0_ddz[i];
          case 2:   for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddy[i] = h_tm0_ddy[i];
          default:  for (CFuint i=0; i<nbStates*nbVar; ++i)  h_tm1_ddx[i] = h_tm0_ddx[i];
        }
        d.ver("solution derivatives update queue 2 from 1.");
      }
      default: break;
    }
  }


  // next executions won't be the first
  m_firstexecute = false;
#if 0
  d.log("output?");
  {
    DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
    DataHandle< CFreal > h_tm0_sol = s_tm0_sol.getDataHandle();
    DataHandle< CFreal > h_tm0_ddx = s_tm0_ddx.getDataHandle();
    DataHandle< CFreal > h_tm0_ddy = s_tm0_ddy.getDataHandle();
    DataHandle< CFreal > h_tm0_ddz = s_tm0_ddz.getDataHandle();
    const ConnectivityTable< CFuint >& elems = *(MeshDataStack::getActive()->getFilteredTrsList("inner")[0])->getGeo2NodesConn();
    std::ofstream f(d.getFilename("debug-solqueue.plt").c_str(),std::ios::trunc);
    f << "VARIABLES = x y z P U V W" << std::endl;

    // solution
    f << "ZONE T=\"sol\" DATAPACKING=POINT ZONETYPE=FETETRAHEDRON  N=" << nbNodes << " E=" << elems.nbRows() << std::endl;

    // point cloud
    for (CFuint n=0; n<nbNodes; ++n) {
      f << *h_nodes[n];
      for (CFuint j=0; j<nbVar; ++j)
        f << ' ' << h_tm0_sol(n,j,nbVar);
      f << std::endl;
    }

    // connectivity
    for (CFuint ic=0; ic<elems.nbRows(); ++ic) {
      for (CFuint i=0; i<elems.nbCols(ic); ++i)
        f << ' ' << elems(ic,i)+1;
      f << std::endl;
    }

    // ddx
    f << "ZONE T=\"ddx\" DATAPACKING=POINT ZONETYPE=FETETRAHEDRON  N=" << nbNodes << " E=" << elems.nbRows() << std::endl;

    // point cloud
    for (CFuint n=0; n<nbNodes; ++n) {
      f << *h_nodes[n];
      for (CFuint j=0; j<nbVar; ++j)
        f << ' ' << h_tm0_ddx(n,j,nbVar);
      f << std::endl;
    }

    // connectivity
    for (CFuint ic=0; ic<elems.nbRows(); ++ic) {
      for (CFuint i=0; i<elems.nbCols(ic); ++i)
        f << ' ' << elems(ic,i)+1;
      f << std::endl;
    }

    // ddy
    f << "ZONE T=\"ddy\" DATAPACKING=POINT ZONETYPE=FETETRAHEDRON  N=" << nbNodes << " E=" << elems.nbRows() << std::endl;

    // point cloud
    for (CFuint n=0; n<nbNodes; ++n) {
      f << *h_nodes[n];
      for (CFuint j=0; j<nbVar; ++j)
        f << ' ' << h_tm0_ddy(n,j,nbVar);
      f << std::endl;
    }

    // connectivity
    for (CFuint ic=0; ic<elems.nbRows(); ++ic) {
      for (CFuint i=0; i<elems.nbCols(ic); ++i)
        f << ' ' << elems(ic,i)+1;
      f << std::endl;
    }

    // ddz
    f << "ZONE T=\"ddz\" DATAPACKING=POINT ZONETYPE=FETETRAHEDRON  N=" << nbNodes << " E=" << elems.nbRows() << std::endl;

    // point cloud
    for (CFuint n=0; n<nbNodes; ++n) {
      f << *h_nodes[n];
      for (CFuint j=0; j<nbVar; ++j)
        f << ' ' << h_tm0_ddz(n,j,nbVar);
      f << std::endl;
    }

    // connectivity
    for (CFuint ic=0; ic<elems.nbRows(); ++ic) {
      for (CFuint i=0; i<elems.nbCols(ic); ++i)
        f << ' ' << elems(ic,i)+1;
      f << std::endl;
    }

    f.close();
  }
  d.ver("output!");
CF_DEBUG_EXIT;
#endif
}

//////////////////////////////////////////////////////////////////////////////

void ComMultipleStates::setStatesDerivatives(TopologicalRegionSet& trs, DataHandle< CFreal >& ddx, DataHandle< CFreal >& ddy, DataHandle< CFreal >& ddz)
{
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< State*,GLOBAL > h_states = s_states.getDataHandle();


  // auxiliary variables
  std::vector< CFreal > vdx(nbVar,0.),  // derivatives, x-component
                        vdy(nbVar,0.),  // ... y-component
                        vdz(nbVar,0.);  // ... z-component
  std::vector< CFuint >& nodes = *trs.getNodesInTrs();
  const ConnectivityTable< CFuint >& elems = *trs.getGeo2NodesConn();
  AElement* e(nbDim==2? (AElement*) new ElementTriangle(2) :
                        (AElement*) new ElementTetrahedron(3) );
  std::vector< Node* > vnodes(e->N);
  std::vector< State* > vstates(e->S);


  // reset states derivatives (per node)
  for (CFuint i=0; i<nodes.size(); ++i) {
    const CFuint n = nodes[i];
    for (CFuint j=0; j<nbVar; ++j) {
      switch (nbDim) {
        case 3:  ddz(n,j,nbVar) = 0.;
        case 2:  ddy(n,j,nbVar) = 0.;
        case 1:  ddx(n,j,nbVar) = 0.;
        default: break;
      }
    }
  }


  // calculate states derivatives, and distribute weighing with nodal size
  // (per element)
  for (CFuint c=0; c<elems.nbRows(); ++c) {

    for (CFuint j=0; j<e->N; ++j)  vnodes[j]  = h_nodes[elems(c,j)];
    for (CFuint j=0; j<e->S; ++j)  vstates[j] = h_states[elems(c,j)];
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
          case 3:  ddz( elems(c,i) ,j,nbVar) += vdz[j]*w;
          case 2:  ddy( elems(c,i) ,j,nbVar) += vdy[j]*w;
          case 1:  ddx( elems(c,i) ,j,nbVar) += vdx[j]*w;
          default: break;
        }
      }
    }

  }


  // de-weight the derivatives (per node)
  DataHandle< CFreal > h_nvolume = s_nvolume.getDataHandle();
  for (CFuint i=0; i<nodes.size(); ++i) {
    const CFuint n = nodes[i];
    const CFreal w = h_nvolume[n];
    if (w>PhysicalConstants::_eps) {
      for (CFuint j=0; j<nbVar; ++j) {
        switch (nbDim) {
          case 3:  ddz(n,j,nbVar) /= w;
          case 2:  ddy(n,j,nbVar) /= w;
          case 1:  ddx(n,j,nbVar) /= w;
          default: break;
        }
      }
    }
  }


  // deallocate AElement
  delete e;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

