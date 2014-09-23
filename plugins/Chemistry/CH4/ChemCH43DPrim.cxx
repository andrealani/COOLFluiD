#include "Chemistry/CH4/CH4.hh"
#include <numeric>

#include "ChemCH43DPrim.hh"
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

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ChemCH43DPrim, ConvectiveVarSet, CH4Module, 1> ChemCH43DPrimProvider("ChemCH43DPrim");

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrim::ChemCH43DPrim(Common::SafePtr<Framework::BaseTerm> term) :
  ConvectiveVarSet(term),
    _model(
      Framework::PhysicalModelStack::getActive()->
      getImplementor().d_castTo<ChemCH4PhysicalModel>())
{
  vector<std::string> names(4);
  names[XCH4] = "XCH4";
  names[XO2 ] = "XO2";
  names[XCO2] = "XCO2";
  names[XH2O] = "XH2O";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrim::~ChemCH43DPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrim::computeFlux (const RealVector& vars,
                                 const RealVector& normals)
{
 /// @todo broken after release 2009.3
 throw Common::NotImplementedException (FromHere(),"ChemCH43DPrim::computeFlux()");
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrim::computeFlux (const RealVector& vars)
{
 /// @todo broken after release 2009.3
 throw Common::NotImplementedException (FromHere(),"ChemCH43DPrim::computeFlux()");
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrim::computeStateFlux(const RealVector& vars)
{
 /// @todo broken after release 2009.3
 throw Common::NotImplementedException (FromHere(),"ChemCH43DPrim::computeStateFlux()");
}

/////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrim::computePhysicalData(const State& state, RealVector& data)
{
   /// @todo broken after release 2009.3
 throw Common::NotImplementedException (FromHere(),"ChemCH43DPrim::computePhysicalData()");
}


//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
