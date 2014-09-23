#include "RKRD/RungeKuttaRD.hh"

#include "Common/BadValueException.hh"

#include "RKRD/RKRDData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<RKRDData>, RKRDData, RKRDModule> nullRKRDComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void RKRDData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint > ("Order","Order of the Runge-Kutta scheme");
}

//////////////////////////////////////////////////////////////////////////////

RKRDData::RKRDData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_order(2),
    m_alpha(),
    m_beta()
{
  addConfigOptionsTo(this);

  setParameter("Order",&m_order);
}

//////////////////////////////////////////////////////////////////////////////

RKRDData::~RKRDData()
{
}

//////////////////////////////////////////////////////////////////////////////

void RKRDData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  switch( m_order )
  {
  case 1:

    m_alpha.resize(1,1);
    m_beta.resize(1,1);

    m_alpha(0,0) = 1.;
    m_beta (0,0) = 0.;

    break;

  case 2:

    m_alpha.resize(2,2);
    m_beta.resize(2,2);

    m_alpha(0,0) =  0.0;
    m_alpha(0,1) =  0.5;
    m_alpha(1,0) =  1.0;
    m_alpha(1,1) =  0.5;

    m_beta (0,0) =  0.0;
    m_beta (0,1) = -1.0;
    m_beta (1,0) =  0.0;
    m_beta (1,1) =  1.0;

    break;

  case 3:

    m_alpha.resize(3,3);
    m_beta.resize(3,3);

    m_alpha(0,0) = 0. ;
    m_alpha(0,1) = 0.25 ;
    m_alpha(0,2) = 1./6. ;

    m_alpha(1,0) = 0. ;
    m_alpha(1,1) = 0. ;
    m_alpha(1,2) = 1./6. ;

    m_alpha(2,0) = 1. ;
    m_alpha(2,1) = 0.25 ;
    m_alpha(2,2) = 2./3. ;

    m_beta(0,0) = 0.  ;
    m_beta(0,1) = -0.5 ;
    m_beta(0,2) = -2. ;

    m_beta(1,0) = 0. ;
    m_beta(1,1) = 0. ;
    m_beta(1,2) = 0. ;

    m_beta(2,0) = 0. ;
    m_beta(2,1) = 0.5 ;
    m_beta(2,2) = 2. ;

    break ;

  default:
    throw Common::BadValueException(FromHere(), "Runge-Kutta order " + StringOps::to_str(m_order) + " not supported\n");
  }

  std::cout << "Alpha " << " = " << m_alpha << std::endl;
  std::cout << "Beta "  << " = " << m_beta << std::endl;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

