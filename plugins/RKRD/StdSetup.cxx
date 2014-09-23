#include "RKRD/RungeKuttaRD.hh"
#include "RKRD/StdSetup.hh"

#include "MathTools/RealVector.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/GeometricEntityPool.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, RKRDData, RKRDModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(const std::string& name) :
  RKRDCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_kstates("kstates"),
  socket_median_areas("median_areas"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;



  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbstates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  rhs.resize(nbstates*nbEqs);
  rhs = 0.0;



  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  updateCoeff.resize(nbstates);



  const CFuint order = getMethodData().getOrder();


  DataHandle<RealMatrix> kstates = socket_kstates.getDataHandle();
  kstates.resize(nbstates);
  for ( CFuint i = 0; i < nbstates; ++i)
  {
    kstates[i].resize(nbEqs,order);
  }



  // compute median dual cell areas
  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  DataHandle<CFreal> median_areas = socket_median_areas.getDataHandle();
  median_areas.resize(nbstates);
  median_areas = 0.;



  // builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilder;
  geoBuilder.setup();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;



  for(CFuint iGeoEnt = 0; iGeoEnt < nbCells; ++iGeoEnt)
  {
    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity& cell = *geoBuilder.buildGE();

    std::vector<State*>& states = *cell.getStates();
    const CFuint cell_nbstates = states.size();

    const CFreal med_area = cell.computeVolume() / cell_nbstates;

    for ( CFuint s = 0; s < cell_nbstates; ++s)
      median_areas[ states[s]->getLocalID() ] += med_area;

    //release the GeometricEntity
    geoBuilder.releaseGE();
  }



  geoBuilder.unsetup();


}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSetup::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_kstates);
  result.push_back(&socket_median_areas);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD
