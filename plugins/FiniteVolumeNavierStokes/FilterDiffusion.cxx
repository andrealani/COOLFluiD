#include "Common/PE.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "MathTools/MathConsts.hh"
#include "FiniteVolumeNavierStokes/FilterDiffusion.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/EulerTerm.hh"
#include "FiniteVolume/CellCenterFVM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FilterDiffusion, DataProcessingData, FiniteVolumeNavierStokesModule>
filterDiffusionProvider("FilterDiffusion");

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusion::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("FreeStreamTotalEnthalpy","Free stream total enthalpy.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("TotalEnthalpyDeviation","Relative error (in percentage) from free stream total enthalpy.");
}

//////////////////////////////////////////////////////////////////////////////
      
FilterDiffusion::FilterDiffusion(const std::string& name) :
  DataProcessingCom(name),
  socket_activeDiffusion("activeDiffusion"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  m_geoBuilder(),
  m_physicalData()
{
  addConfigOptionsTo(this);
  
  m_freeStreamH = 0.;
  setParameter("FreeStreamTotalEnthalpy",&m_freeStreamH);
  
  m_deviationH = 0.;
  setParameter("TotalEnthalpyDeviation",&m_deviationH);
}
      
//////////////////////////////////////////////////////////////////////////////

FilterDiffusion::~FilterDiffusion()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FilterDiffusion::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_activeDiffusion);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_gstates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusion::setup()
{
  CFAUTOTRACE;
  
  // geometry builder setup
  m_geoBuilder.setup();
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_physicalData);
  cf_assert(m_physicalData.size() > 0);
}

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusion::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusion::execute()
{
  CFAUTOTRACE;
  
  if (SubSystemStatusStack::getActive()->getNbIter() > getMethodData().getStartIter()) {
    CFLog(VERBOSE, "FilterDiffusion::execute()" << "\n");  
    DataHandle<CFreal> activeDiffusion = socket_activeDiffusion.getDataHandle();
    
    //  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");  
    //     CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
    //     // AL: it's better to set the cells here because of the namespace is switched during the computation
    //     // this would point to the wrong TRS !!
    //     geoData.trs = cells;  
    //     const CFuint nbCells = cells->getLocalNbGeoEnts();
    
    // suppose that just one space method is available
    SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
    cf_assert(fvmcc.isNotNull());
    SafePtr<ConvectiveVarSet> updateVarSet = fvmcc->getData()->getUpdateVar();
    cf_assert(updateVarSet.isNotNull());
    
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    const CFuint nbCells = states.size();
    
    const CFreal factorH = 100./m_freeStreamH;
    cf_assert(m_freeStreamH > 0.);
    
    for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
      CFLog( VERBOSE, "iCell "<< iCell << "\n");
      activeDiffusion[iCell] = 1.;
      
      // set the builder data and build the GeometricEntity
      // geoData.idx = iCell;
      // GeometricEntity* currCell = m_geoBuilder.buildGE();
      //updateVarSet->computePhysicalData(*currCell->getState(0), m_physicalData);
      
      updateVarSet->computePhysicalData(*states[iCell], m_physicalData);
      
      const CFreal errorH = std::abs(m_physicalData[EulerTerm::H] - m_freeStreamH)*factorH;
      if (errorH < m_deviationH) {
	activeDiffusion[iCell] = 0.;
      }
      
      // m_geoBuilder.releaseGE();
    }
    
    CFLog(VERBOSE, "FilterDiffusion::execute() END\n");  
  }
}

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusion::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
