
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/Element.hh"
#include "Muffin/SystemPLaS.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemPLaS,MuffinData,MuffinModule > cSystemPLaSProvider("PLaS");

//////////////////////////////////////////////////////////////////////////////

SystemPLaS::SystemPLaS(const std::string& name) :
    System(name),
    s_mn_voidfraction("NodalVoidFraction"),  // socket sources
    s_evolume("ElementVolumes"),             // ...
    s_ienormals("ElementNormals"),           // ...
    s_benormals("BoundariesNormals")         // ...
{
  CFAUTOTRACE;

  m_requires_lss = false;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemPLaS");
  addConfigOptionsTo(this);

  m_plas_name = "";
  setParameter("Method",&m_plas_name);
}

//////////////////////////////////////////////////////////////////////////////

void SystemPLaS::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("Method","PLaSTracking DataProcessingMethod name (default \"\", first found");
}

//////////////////////////////////////////////////////////////////////////////

void SystemPLaS::setup()
{
  CFAUTOTRACE;
  System::setup();
  MuffinData& d = getMethodData();

  log("get PLaS DataProcessingMethod...");
  std::vector< SafePtr< DataProcessingMethod > > m = MethodRegistry::getInstance().getAllMethods< DataProcessingMethod >(d.getNamespace());
  for (CFuint i=0; i<m.size() && m_plas.isNull(); ++i) {
    if (!m_plas_name.length() || m[i]->getName()==m_plas_name) {
      try {
        SafePtr< PLaS::PLaSTracking > p = m[i].d_castTo< PLaS::PLaSTracking >();
        // at this point cast is successful. nothing happens when name is set
        // but differs from current one; if name is not set and cast ok, or name
        // was set and is equal to this one's name, set pointer to it and exit.
        m_plas = p;
      }
      catch (FailedCastException& e) {
      }
    }
  }
  if (m_plas.isNull())
    err("DataProcessingMethod: no PLaSTracking method found!");
  else
    log("DataProcessingMethod: " + std::string(m_plas->getName()));
  ver("get PLaS DataProcessingMethod.");


  // get DataHandle sinks
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< std::pair< CFuint,CFuint > > h_faceneighcell = s_faceneighcell.getDataHandle();
  const CFuint nbNodes = h_nodes.size();
  const CFuint nbDim = h_nodes[0]->size();


  log("set initial void fraction...");
  DataHandle< CFreal > h_mn_voidfraction = s_mn_voidfraction.getDataHandle();
  h_mn_voidfraction.resize(nbNodes);
  h_mn_voidfraction = 0.;
  ver("set initial void fraction.");


  log("set flow properties...");
  m_plas->setFlowProperties(
    d.m_density, d.m_kviscosity, d.m_status->getDTDim(),
    d.m_cp, d.m_k, nbDim+1 );
  ver("set flow properties.");


  log("set 'inner' elements volumes and normals...");

  DataHandle< CFreal > h_evolume = s_evolume.getDataHandle();
  DataHandle< std::vector< RealVector > > h_ienormals = s_ienormals.getDataHandle();
  DataHandle< std::vector< RealVector > > h_benormals = s_benormals.getDataHandle();
  const std::vector< SafePtr< TopologicalRegionSet > > itrs = MeshDataStack::getActive()->getFilteredTrsList("inner");
  const std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");
  cf_always_assert_desc("no appropriate 'inner' TRS found",itrs.size());

  // create AElement
  AElement* e(nbDim==2? (AElement*) new ElementTriangle(nbDim) :
                        (AElement*) new ElementTetrahedron(nbDim) );
  std::vector< Node* > vnodes(e->N);

  // allocate DataHandles
  h_evolume.resize(itrs[0]->getLocalNbGeoEnts());
  h_ienormals.resize(itrs[0]->getLocalNbGeoEnts());
  for (CFuint i=0; i<h_ienormals.size(); ++i)
    h_ienormals[i].resize(nbDim);
  h_benormals.resize(btrs.size());
  for (CFuint i=0; i<btrs.size(); ++i)
    h_benormals[i].assign(nbDim,RealVector(0.,btrs[i]->getLocalNbGeoEnts()));

  // for all elements, set element and nodal volumes and faces normals
  const ConnectivityTable< CFuint >& elems = *itrs[0]->getGeo2NodesConn();
  for (CFuint c=0; c<elems.nbRows(); ++c) {
    for (CFuint j=0; j<e->N; ++j)
      vnodes[j] = h_nodes[elems(c,j)];
    h_evolume[c] = e->element(vnodes);
    switch (nbDim) {
      case 3:  {  h_ienormals[c][ZZ].resize(e->F,0.);  }
      case 2:  {  h_ienormals[c][YY].resize(e->F,0.);  }
      case 1:  {  h_ienormals[c][XX].resize(e->F,0.);  }
    }
    for (CFuint f=0; f<e->F; ++f) {
      switch (nbDim) {
        case 3:  {  h_ienormals[c][ZZ][f] = e->normal[f][ZZ];  }
        case 2:  {  h_ienormals[c][YY][f] = e->normal[f][YY];  }
        case 1:  {  h_ienormals[c][XX][f] = e->normal[f][XX];  }
      }
    }
  }

  // deallocate AElement
  delete e;

  ver("set 'inner' elements volumes and normals.");


  log("set 'boundary' elements volumes and normals...");
  // set faces normals based on element normals, converting face number from
  // coolfluid/gambit to PLaS-style
  for (CFuint ib=0; ib<h_benormals.size(); ++ib) {
    for (CFuint idx=0; idx<h_benormals[ib][XX].size(); ++idx) {
      const CFuint id = btrs[ib]->getLocalGeoID(idx);
      const CFuint ic = h_faceneighcell[id].first;
      int f = (int) h_faceneighcell[id].second;
      f = (nbDim==2?
        (f==0? 2 : (f==1? 0 : (f==2? 1 : -1 ))) :
        (f==0? 3 : (f==1? 2 : (f==2? 0 : (f==3? 1 : -1 )))) );
      cf_assert_desc("problem converting face numbering!",f>=0);
      for (CFuint i=0; i<nbDim; ++i)
        h_benormals[ib][i][idx] = h_ienormals[ic][i][f];
    }
  }
  ver("set 'boundary' elements volumes and normals.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemPLaS::execute()
{
  CFAUTOTRACE;


  // check there is a timestep set
  SafePtr< SubSystemStatus > status = getMethodData().m_status;
  const CFreal dt = status->getDTDim();
  if (dt<MathTools::MathConsts::CFrealEps()) {
    err("dt is 0., Loop is likely not correct");
  }


  log("process...");
  m_plas->setIterationProperties(
    status->getNbIter(),
    status->getCurrentTimeDim(),
    status->getDTDim(),
    1 );
  m_plas->unblock();  // m_plas is mine, all mine!
  m_plas->process();  // do it, yeah
  m_plas->block();    // block intruders
  ver("process.");


  log("update void fraction...");
  const PLAS_PHASE_DATA* m_plas_phase = m_plas->getPhaseData();
  DataHandle< CFreal > h_mn_voidfraction = s_mn_voidfraction.getDataHandle();
  for (CFuint n=0; n<h_mn_voidfraction.size(); ++n) {
    const double vf = m_plas_phase[n].volFrac;
    h_mn_voidfraction[n] = (CFreal) std::min(1.,std::max(0.,vf));
  }
  ver("update void fraction.");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

