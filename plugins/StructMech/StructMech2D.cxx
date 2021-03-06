#include "StructMech/StructMech.hh"
#include "Framework/State.hh"
#include "StructMech2D.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech2D,
               PhysicalModelImpl,
	       StructMechModule,
               1>
structMech2DProvider("StructMech2D");

//////////////////////////////////////////////////////////////////////////////

StructMech2D::StructMech2D(const std::string& name)
  : StructMechPhysicalModel(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StructMech2D::~StructMech2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech2D::configure ( Config::ConfigArgs& args )
{
  StructMechPhysicalModel::configure(args);

  /// @note Should be here, but is not yet being used.
  _jacobians = std::vector<RealMatrix>(getDimension());
  for (CFuint i = 0; i < getDimension(); ++i) {
    _jacobians[i] = RealMatrix(getNbEquations(), getNbEquations());
    _jacobians[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

CFuint StructMech2D::getDimension() const
{
  return 2;
}

//////////////////////////////////////////////////////////////////////////////

CFuint StructMech2D::getNbEquations() const
{
  return 2;
}

//////////////////////////////////////////////////////////////////////////////

bool StructMech2D::validate(const State& state) const
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
