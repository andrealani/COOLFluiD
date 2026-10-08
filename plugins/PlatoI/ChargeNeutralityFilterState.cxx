#include "PlatoI/ChargeNeutralityFilterState.hh"
#include "PlatoI/PlatoLibrary.hh"
#include "PlatoI/Plato.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Plato {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ChargeNeutralityFilterState,
                            FilterState,
                            PlatoModule,
                            1>
chargeNeutralityFilterStateProvider("ChargeNeutrality");

//////////////////////////////////////////////////////////////////////////////

void ChargeNeutralityFilterState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >
    ("DensityVariables","Log if the state stores ln(rho_i) (update variables LogRhoivLogTTv), Linear if it stores rho_i (e.g. RhoivtTv); no default, it must match the update variables.");
}

//////////////////////////////////////////////////////////////////////////////

ChargeNeutralityFilterState::ChargeNeutralityFilterState(const std::string& name) :
  FilterState(name),
  m_logDensities(false),
  m_coef(),
  m_isSet(false)
{
  addConfigOptionsTo(this);
  m_densityVariables = "";
  Config::ConfigObject::setParameter("DensityVariables",&m_densityVariables);
}

//////////////////////////////////////////////////////////////////////////////

ChargeNeutralityFilterState::~ChargeNeutralityFilterState()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChargeNeutralityFilterState::configure(Config::ConfigArgs& args)
{
  FilterState::configure(args);

  if (m_densityVariables == "Log") {
    m_logDensities = true;
  }
  else if (m_densityVariables == "Linear") {
    m_logDensities = false;
  }
  else {
    throw BadValueException(FromHere(), "ChargeNeutralityFilterState: set DensityVariables = Log (state stores ln rho_i) or Linear (state stores rho_i) to match the update variables");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ChargeNeutralityFilterState::setCoefficients() const
{
  SafePtr<PlatoLibrary> library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PlatoLibrary>();
  if (library.isNull()) {
    throw BadValueException(FromHere(), "ChargeNeutralityFilterState: the physical property library is not PLATO");
  }
  if (!library->hasChargeNeutrality()) {
    throw BadValueException(FromHere(), "ChargeNeutralityFilterState: needs an ionized mixture with Plato.ChargeNeutrality = true");
  }
  m_coef.resize(library->getChargeNeutralityCoefficients().size());
  m_coef = library->getChargeNeutralityCoefficients();
  m_isSet = true;
}

//////////////////////////////////////////////////////////////////////////////

void ChargeNeutralityFilterState::filter(RealVector& state) const
{
  if (!m_isSet) setCoefficients();

  const CFuint nbSpecies = m_coef.size();
  cf_assert(state.size() > nbSpecies);

  // rho_e = sum over the ions of c_k rho_k
  CFreal rhoE = 0.;
  if (m_logDensities) {
    for (CFuint k = 1; k < nbSpecies; ++k) {
      if (m_coef[k] > 0.) rhoE += m_coef[k]*std::exp(state[k]);
    }
    // keep the stored value where no ion is present (ln 0 undefined)
    if (rhoE > 0.) state[0] = std::log(rhoE);
  }
  else {
    for (CFuint k = 1; k < nbSpecies; ++k) {
      rhoE += m_coef[k]*state[k];
    }
    state[0] = rhoE;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Plato

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
