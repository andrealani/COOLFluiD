#include <cmath>

#include "NEQ.hh"
#include "Euler2DNEQLogRhoivLogTTv.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQLogRhoivLogTTv, ConvectiveVarSet, NEQModule, 1>
euler2DNEQLogRhoivLogTTvProvider("Euler2DNEQLogRhoivLogTTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLogRhoivLogTTv::Euler2DNEQLogRhoivLogTTv(Common::SafePtr<BaseTerm> term) :
  Euler2DNEQRhoivtTv(term),
  m_rvtState(),
  m_parentExtra()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);

  vector<std::string> names(nbSpecies + 3 + nbTv);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "lnrho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "u";
  names[nbSpecies + 1] = "v";
  names[nbSpecies + 2] = "lnT";
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names[nbSpecies + 3 + ie] = "lnTv" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQLogRhoivLogTTv::~Euler2DNEQLogRhoivLogTTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::setup()
{
  Euler2DNEQRhoivtTv::setup();

  m_rvtState.reset(new State());
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::setRhoivtTvState(const State& state)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  State& rvt = *m_rvtState;

  // exp of ln rho_i, u and v copied, exp of ln T and ln Tv
  for (CFuint i = 0; i < state.size(); ++i) {
    rvt[i] = (i < nbSpecies || i >= TID) ? std::exp(state[i]) : state[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::computePhysicalData(const State& state, RealVector& data)
{
  setRhoivtTvState(state);
  Euler2DNEQRhoivtTv::computePhysicalData(*m_rvtState, data);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::computePerturbedPhysicalData(const State& state,
                                                            const RealVector& pdataBkp,
                                                            RealVector& pdata,
                                                            CFuint iVar)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = getTempID(nbSpecies);

  if (iVar < nbSpecies || iVar >= TID) {
    computePhysicalData(state, pdata);
  }
  else {
    // the parent u, v branch only reads u and v from the state, which are the
    // same in both variable sets; the copy must not be passed here, the
    // parent would convert it a second time through computePhysicalData
    Euler2DNEQRhoivtTv::computePerturbedPhysicalData(state, pdataBkp, pdata, iVar);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::setDimensionalValues(const State& state,
                                                    RealVector& result)
{
  Euler2DNEQRhoivtTv::setDimensionalValues(state, result);

  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal lnRhoRef = std::log(refData[EulerTerm::RHO]);
  const CFreal lnTRef = std::log(refData[EulerTerm::T]);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie] + lnRhoRef;
  }
  for (CFuint i = nbSpecies + 2; i < state.size(); ++i) {
    result[i] = state[i] + lnTRef;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::setAdimensionalValues(const State& state,
                                                     RealVector& result)
{
  Euler2DNEQRhoivtTv::setAdimensionalValues(state, result);

  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal lnRhoRef = std::log(refData[EulerTerm::RHO]);
  const CFreal lnTRef = std::log(refData[EulerTerm::T]);
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[ie] = state[ie] - lnRhoRef;
  }
  for (CFuint i = nbSpecies + 2; i < state.size(); ++i) {
    result[i] = state[i] - lnTRef;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::setDimensionalValuesPlusExtraValues(const State& state,
                                                                   RealVector& result,
                                                                   RealVector& extra)
{
  // parent extra values from the copy; its result is overwritten below
  setRhoivtTvState(state);
  Euler2DNEQRhoivtTv::setDimensionalValuesPlusExtraValues(*m_rvtState, result, m_parentExtra);

  Euler2DNEQLogRhoivLogTTv::setDimensionalValues(state, result);

  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint TID = nbSpecies + 2;
  const CFuint nbTemps = state.size() - TID;
  const CFuint nbParentExtra = m_parentExtra.size();

  // parent extras, then rho_i, T and the Tv's
  extra.resize(nbParentExtra + nbSpecies + nbTemps);
  for (CFuint i = 0; i < nbParentExtra; ++i) {
    extra[i] = m_parentExtra[i];
  }
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    extra[nbParentExtra + ie] = (*m_rvtState)[ie]*refData[EulerTerm::RHO];
  }
  for (CFuint i = 0; i < nbTemps; ++i) {
    extra[nbParentExtra + nbSpecies + i] = (*m_rvtState)[TID + i]*refData[EulerTerm::T];
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQLogRhoivLogTTv::getExtraVarNames() const
{
  vector<std::string> names = Euler2DNEQRhoivtTv::getExtraVarNames();

  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names.push_back("rho" + StringOps::to_str(ie));
  }

  const CFuint nbTv = getModel()->getNbScalarVars(1);
  names.push_back("T");
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names.push_back("Tv" + StringOps::to_str(ie));
  }

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQLogRhoivLogTTv::computePressureDerivatives(const State& state,
                                                          RealVector& dp)
{
  setRhoivtTvState(state);
  Euler2DNEQRhoivtTv::computePressureDerivatives(*m_rvtState, dp);

  // dp/d(ln x) = x dp/dx for the partial densities, T and the Tv's
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  for (CFuint i = 0; i < state.size(); ++i) {
    if (i < nbSpecies || i >= nbSpecies + 2) {
      dp[i] *= (*m_rvtState)[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool Euler2DNEQLogRhoivLogTTv::isValid(const RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  for (CFuint i = 0; i < nbSpecies + 3 + nbTv; ++i) {
    if ((i < nbSpecies || i >= nbSpecies + 2) && !std::isfinite(data[i])) {
      CFLog(VERBOSE, "Euler2DNEQLogRhoivLogTTv::isValid() => logarithm " << i << " = " << data[i] << "\n");
      return false;
    }
  }

  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
