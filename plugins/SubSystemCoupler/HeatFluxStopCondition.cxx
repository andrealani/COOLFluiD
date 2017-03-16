#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/VarRegistry.hh"

#include "SubSystemCoupler/HeatFluxStopCondition.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

      namespace SubSystemCoupler {

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< HeatFluxStopCondition,
                             Framework::StopCondition,
			     SubSystemCouplerModule,
                             1>
aHeatFluxStopConditionProvider("HeatFlux");

//////////////////////////////////////////////////////////////////////////////

void HeatFluxStopCondition::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("RelativeNorm","Convergence threshold for L2 norm of the boundary residual.");
}

//////////////////////////////////////////////////////////////////////////////

HeatFluxStopCondition::HeatFluxStopCondition(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_first_norm = MathTools::MathConsts::CFrealMax();

  m_conv_norm = -3.0;
  setParameter("RelativeNorm",&m_conv_norm);
}

//////////////////////////////////////////////////////////////////////////////

HeatFluxStopCondition::~HeatFluxStopCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void HeatFluxStopCondition::configure (const Config::ConfigArgs& args)
{
  StopCondition::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool HeatFluxStopCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool HeatFluxStopCondition::isAchieved(const Framework::ConvergenceStatus& status)
{
  bool return_value = false; // by default we are not converged

  // skip evaluation before doing any iteration
  if ( status.iter == 0 ) return return_value;

  // get the subsystem status
  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();
  Common::SafePtr<VarRegistry> ssys_var_regist = ssys_status->getVarRegistry();

  // get latest the coefficients
  CFreal& SqrDeltaFlux = ssys_var_regist->getVar<CFreal>("SqrDeltaFlux");

  // compute the norm
  const CFreal l2norm = std::log10 ( std::sqrt(SqrDeltaFlux) );

  // store the first norm on the first iter
  if (status.iter == 1)
      m_first_norm = l2norm;

  // check for convergence
  if ( m_first_norm - l2norm < m_conv_norm )
      return_value = true;

  // reset the norm for BC to accumulate
  SqrDeltaFlux = 0.;

  return return_value;
}

//////////////////////////////////////////////////////////////////////////////

      } // end of namespace SubSystemCoupler

    } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

