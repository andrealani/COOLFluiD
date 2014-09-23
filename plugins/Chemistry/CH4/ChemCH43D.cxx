#include "Chemistry/CH4/CH4.hh"
#include "Framework/State.hh"
#include "ChemCH43D.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ChemCH43D, PhysicalModelImpl, CH4Module, 1> chemCH43DProvider("ChemCH43D");

//////////////////////////////////////////////////////////////////////////////

ChemCH43D::ChemCH43D(const std::string& name) : ChemCH4PhysicalModel(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ChemCH43D::~ChemCH43D()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43D::configure ( Config::ConfigArgs& args )
{
  ChemCH4PhysicalModel::configure(args);

  /// @note Should be here, but is not yet being used.
  _jacobians = std::vector<RealMatrix>(getDimension());
  for (CFuint i = 0; i < getDimension(); ++i) {
    _jacobians[i] = RealMatrix(getNbEquations(), getNbEquations());
    _jacobians[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint ChemCH43D::getDimension() const
{
  return 3;
}

//////////////////////////////////////////////////////////////////////////////

CFuint ChemCH43D::getNbEquations() const
{
  return 4;
}

//////////////////////////////////////////////////////////////////////////////

bool ChemCH43D::validate(const State& state) const
{
  return (state[0] >= 0.0 && state[0] <= 1.0) &&
         (state[1] >= 0.0 && state[1] <= 1.0) &&
         (state[2] >= 0.0 && state[2] <= 1.0) &&
         (state[3] >= 0.0 && state[3] <= 1.0);
}

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
