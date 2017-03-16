#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/VarRegistry.cxx"

#include "AeroCoef/AeroCoef.hh"
#include "AeroCoef/AeroCoefStopCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AeroCoefStopCondition,
                            Framework::StopCondition,
                            AeroCoefModule,
                            1>
aAeroCoefStopConditionProvider("AeroCoef");

//////////////////////////////////////////////////////////////////////////////

void AeroCoefStopCondition::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("ConvCL","Convergence threshold for CL.");
   options.addConfigOption< CFreal >("ConvCD","Convergence threshold for CD.");
   options.addConfigOption< CFreal >("ConvCM","Convergence threshold for CM.");

   options.addConfigOption< bool >  ("CheckCL","Convergence threshold for CL.");
   options.addConfigOption< bool >  ("CheckCD","Convergence threshold for CD.");
   options.addConfigOption< bool >  ("CheckCM","Convergence threshold for CM.");

   options.addConfigOption< CFuint >("NbIters","Number of iterations to check that the coefficients do not change their value within each threshold.");
   options.addConfigOption< CFuint >("MaxIter","Maximum number of iterations.");
}

//////////////////////////////////////////////////////////////////////////////

AeroCoefStopCondition::AeroCoefStopCondition(const std::string& name) :
  StopCondition(name)
{
  addConfigOptionsTo(this);

  m_check_cl = true;
  setParameter("CheckCL",&m_check_cl);

  m_check_cd = true;
  setParameter("CheckCD",&m_check_cd);

  m_check_cm = true;
  setParameter("CheckCM",&m_check_cm);

  m_conv_cl = 1E-3;
  setParameter("ConvCL",&m_conv_cl);

  m_conv_cd = 0.0;
  setParameter("ConvCD",&m_conv_cd);

  m_conv_cm = 0.0;
  setParameter("ConvCM",&m_conv_cm);

  m_nb_iter = 3;
  setParameter("NbIters",&m_nb_iter);

  m_max_iter = 10000;
  setParameter("MaxIter",&m_max_iter);
}

//////////////////////////////////////////////////////////////////////////////

AeroCoefStopCondition::~AeroCoefStopCondition()
{
}

//////////////////////////////////////////////////////////////////////////////

void AeroCoefStopCondition::configure ( Config::ConfigArgs& args )
{
  StopCondition::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

bool AeroCoefStopCondition::IsGlobal () const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

bool AeroCoefStopCondition::isAchieved(const Framework::ConvergenceStatus& status)
{
  Common::SafePtr<SubSystemStatus> ssys_status =
    SubSystemStatusStack::getActive();

  Common::SafePtr<VarRegistry> ssys_var_regist =
    ssys_status->getVarRegistry();

  // get latest the coefficients
  StoreCoefs last_coeffs;
  last_coeffs.CL = ssys_var_regist->getVar<CFreal>("CL");
  last_coeffs.CD = ssys_var_regist->getVar<CFreal>("CD");
  last_coeffs.CM = ssys_var_regist->getVar<CFreal>("CM");


  // store latest the coefficients
  m_past_coeffs.push_back(last_coeffs);

  // keep only the lates m_nb_iter iterations
  if (m_past_coeffs.size() > m_nb_iter)
    m_past_coeffs.pop_front();

  // skip check if still not enough information
  if (m_past_coeffs.size() < m_nb_iter)
    return false;

  //If the number of iterations is bigger than MaxIterations, then should stop
  if (m_max_iter == SubSystemStatusStack::getActive()->getNbIter())
  { CFout << "!!! Max number iterations [" << m_max_iter << "] reached !!!\n";

    return true;
}

  // check if we have convergence
  bool return_value = true;
  std::deque<StoreCoefs>::const_iterator this_value = m_past_coeffs.begin();
  std::deque<StoreCoefs>::const_iterator past_value = this_value;
  for ( ++this_value; this_value != m_past_coeffs.end(); ++this_value)
  {
    cf_assert( this_value != past_value );

    // check CL if threshold is valid
    if ( m_check_cl )
      return_value &= (std::abs(past_value->CL - this_value->CL) < m_conv_cl);
    // check CD if threshold is valid
    if ( m_check_cd )
      return_value &= (std::abs(past_value->CD - this_value->CD) < m_conv_cd);
    // check CM if threshold is valid
    if ( m_check_cm )
      return_value &= (std::abs(past_value->CM - this_value->CM) < m_conv_cm);

//     if (return_value){
//         CFout << " Threshold reached in   [" << SubSystemStatusStack::getActive()->getNbIter() << "] iteration\n";
//         CFout << " error CL  =" <<  std::abs(past_value->CL - this_value->CL)<< "\n";
//         CFout << " error CD  =" <<  std::abs(past_value->CD - this_value->CD)<< "\n";
//         CFout << " error CM  =" <<  std::abs(past_value->CM - this_value->CM)<< "\n";
//     }

    past_value = this_value;

  }

  return return_value;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
