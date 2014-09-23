#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/FilterDiffusionByTotEnthalpy.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<FilterDiffusionByTotEnthalpy,
                       CellCenterFVMData,
                       EquationFilter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
filterDiffusionByTotEnthalpyProv("FilterDiffusionByTotEnthalpy");

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusionByTotEnthalpy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("FreeStreamTotalEnthalpy","Free stream total enthalpy.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("TotalEnthalpyDeviation","Relative error (in percentage) from free stream total enthalpy.");
  options.addConfigOption< bool >("FromPostprocessing","Only interpret data coming from postprocessing, do not compute.");
}
      
//////////////////////////////////////////////////////////////////////////////

FilterDiffusionByTotEnthalpy::FilterDiffusionByTotEnthalpy(const std::string& name) :
  FVMCC_EquationFilter(name),
  m_factorH(0.)
{
  addConfigOptionsTo(this);
  
  m_freeStreamH = 0.;
  setParameter("FreeStreamTotalEnthalpy",&m_freeStreamH);
  
  m_deviationH = 0.;
  setParameter("TotalEnthalpyDeviation",&m_deviationH);
  
  m_fromPostProcessing = false;
  setParameter("FromPostprocessing",&m_fromPostProcessing);
}
      
//////////////////////////////////////////////////////////////////////////////

FilterDiffusionByTotEnthalpy::~FilterDiffusionByTotEnthalpy()
{
}

//////////////////////////////////////////////////////////////////////////////

void FilterDiffusionByTotEnthalpy::setup()
{
  if (!m_fromPostProcessing) {
    cf_always_assert(m_freeStreamH > 0.);
    m_factorH = 100./m_freeStreamH;
    
    cf_always_assert(m_deviationH > 0.);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FilterDiffusionByTotEnthalpy::unsetup()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void FilterDiffusionByTotEnthalpy::configure ( Config::ConfigArgs& args )
{
  FVMCC_EquationFilter::configure(args);
}

//////////////////////////////////////////////////////////////////////////////
      
void FilterDiffusionByTotEnthalpy::reset()
{
 //  if (!m_fromPostProcessing) {
//     socket_activeDiffusion.getDataHandle() = 1.;
//   }
}

//////////////////////////////////////////////////////////////////////////////

bool FilterDiffusionByTotEnthalpy::filterOnGeo(Framework::GeometricEntity *const geo)
{
  DataHandle<CFreal> activeDiffusion = socket_activeDiffusion.getDataHandle();
  const State* stateL = geo->getState(0);
  const State* stateR = geo->getState(1);
  
  if (!m_fromPostProcessing) {
    bool result = true;
    if (SubSystemStatusStack::getActive()->getNbIter() > m_startIter) {
      // geo is a face here
      // flux for the right and left state
      CellCenterFVMData& data = getMethodData(); 
      SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
      vector<RealVector>& pdata = polyRec->getExtrapolatedPhysicaData();
      
      activeDiffusion[stateL->getLocalID()] = 1.;
      const CFreal erroHL = std::abs(pdata[0][EulerTerm::H] - m_freeStreamH)*m_factorH;
      if (erroHL < m_deviationH) {
	activeDiffusion[stateL->getLocalID()] = 0.;
	result = false;
      }
      
      if (!stateR->isGhost()) {
	activeDiffusion[stateR->getLocalID()] = 1.;
	const CFreal errorHR = std::abs(pdata[1][EulerTerm::H] - m_freeStreamH)*m_factorH;
	if (errorHR < m_deviationH) {
	  activeDiffusion[stateR->getLocalID()] = 0.;
	  result = false;
	}
      }
    }
    return result;
  }
  
  return (activeDiffusion[stateL->getLocalID()] > 0. || 
	  (!stateR->isGhost() && activeDiffusion[stateR->getLocalID()] > 0.));
}
      
//////////////////////////////////////////////////////////////////////////////
     
} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
