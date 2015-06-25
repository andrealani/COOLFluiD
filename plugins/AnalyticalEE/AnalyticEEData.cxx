#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "AnalyticalEE/AnalyticalEE.hh"
#include "AnalyticalEE/AnalyticEEData.hh"

#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<AnalyticEEData>, AnalyticEEData, AnalyticalEEModule> nullAnalyticEEComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void AnalyticEEData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("updateVar","update variables set");
}

//////////////////////////////////////////////////////////////////////////////

AnalyticEEData::AnalyticEEData(Common::SafePtr<Method> owner): ErrorEstimatorData(owner),
  m_updateVarSet()
{ 
  addConfigOptionsTo(this);
  m_updateVarStr = "Null";
  setParameter("updateVar",&m_updateVarStr);
}

//////////////////////////////////////////////////////////////////////////////

AnalyticEEData::~AnalyticEEData()
{
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticEEData::configure ( Config::ConfigArgs& args )
{
  MethodData::configure(args);
  
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = 
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string provider = "Null";
  if (m_updateVarStr != "Null") {
    provider = (physModel->getConvectiveName() != "Null") ?
      (physModel->getConvectiveName() + m_updateVarStr ) :
      (physModel->getDiffusiveName() + m_updateVarStr ) ;
  }

  CFLog(VERBOSE, "TecWriterData::UpdateVarStr = " << provider << "\n");

  m_updateVarSet.reset(Environment::Factory<ConvectiveVarSet>::getInstance().
                      getProvider(provider)->create(physModel->getImplementor()->getConvectiveTerm()));

  cf_assert(m_updateVarSet.isNotNull());


  // configure the expression for the analytical solution

/*  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  */
//   if ( m_vFunction.getNbFuncs() != physModel->getNbEq() )
//   {
//     std::ostringstream msg;
//     msg << "Wrong number of functions configured in Analytical Error Estimator.\n"
//         << "Should be the same has the number of equations, "
//         << physModel->getNbEq() << ", but found "
//         << m_vFunction.getNbFuncs() << ".\n";
//     throw BadValueException (FromHere(), msg.str() );
//   }
// 
//   if ( m_vFunction.getNbVars() != physModel->getDim() )
//   {
//     std::ostringstream msg;
//     msg << "Wrong number of variables configured in Analytical Error Estimator.\n"
//         << "Should be the same has the number of dimensions, "
//         << physModel->getDim() << ", but found "
//         << m_vFunction.getNbVars() << ".\n";
//     throw BadValueException (FromHere(), msg.str() );
//   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

