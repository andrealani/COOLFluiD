#include "StructMech/StructMech.hh"
#include <numeric>

#include "StructMech3DDisp.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech3DDisp, ConvectiveVarSet, StructMechModule, 1> structMech3DDispProvider("StructMech3DDisp");

//////////////////////////////////////////////////////////////////////////////

StructMech3DDisp::StructMech3DDisp(Common::SafePtr<Framework::BaseTerm> term) :
  ConvectiveVarSet(term),
    _model(CFNULL)
{
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

StructMech3DDisp::~StructMech3DDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech3DDisp::setup()
{
  //Framework::ConvectiveVarSet::setup();

  _model = Framework::PhysicalModelStack::getActive()->
    getImplementor().d_castTo<StructMechPhysicalModel>();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
