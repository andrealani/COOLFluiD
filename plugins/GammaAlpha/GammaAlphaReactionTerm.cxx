#include "GammaAlphaReactionTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GammaAlpha {

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaReactionTerm::defineConfigOptions(Config::OptionList& options)
{
//   options.addConfigOption< std::string >("Model","Type of kOmega model to use (kOmega/BSL/SST).");
}

//////////////////////////////////////////////////////////////////////////////

GammaAlphaReactionTerm::GammaAlphaReactionTerm(const std::string& name) :
  BaseTerm(name)
{
   addConfigOptionsTo(this);

//    _modelStr = "kOmega";
//    setParameter("Model",&_modelStr);

}

//////////////////////////////////////////////////////////////////////////////

GammaAlphaReactionTerm::~GammaAlphaReactionTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaReactionTerm::configure ( Config::ConfigArgs& args )
{
  BaseTerm::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaReactionTerm::resizePhysicalData(RealVector& physicalData)
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  physicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

void GammaAlphaReactionTerm::setupPhysicalData()
{
  cf_assert(getDataSize() > 0);

  // set the size of each physical data in the StatesData
  
  resizePhysicalData(m_physicalData);
  resizePhysicalData(m_refPhysicalData);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GammaAlpha

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
