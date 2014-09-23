#include "Chemistry/CH4/CH4.hh"
#include "ChemCH43DPrimDiffusive.hh"
#include "ChemCH43DDiffusiveVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Chemistry {

      namespace CH4 {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ChemCH43DPrimDiffusive, DiffusiveVarSet, CH4Module, 2> ChemCH43DPrimDiffusiveProvider("ChemCH43DPrimDiffusive");

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimDiffusive::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<CFreal> >("DiffusiveCoeff","Coefficients for the diffusive term");
}

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrimDiffusive::ChemCH43DPrimDiffusive(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  ChemCH43DDiffusiveVarSet(name, model)
{
   addConfigOptionsTo(this);
  vector<std::string> names(4);
  names[XCH4] = "XCH4";
  names[XO2 ] = "XO2";
  names[XCO2] = "XCO2";
  names[XH2O] = "XH2O";
  setVarNames(names);

  const CFuint size = getModel()->getNbEquations();

  _vecD.resize(size);
  _vecD[XCH4] = 1.0;
  _vecD[XO2 ] = 1.0;
  _vecD[XCO2] = 1.0;
  _vecD[XH2O] = 1.0;
   setParameter("DiffusiveCoeff",&_vecD);


}

//////////////////////////////////////////////////////////////////////////////

ChemCH43DPrimDiffusive::~ChemCH43DPrimDiffusive()
{
}

//////////////////////////////////////////////////////////////////////////////

void ChemCH43DPrimDiffusive::configure ( Config::ConfigArgs& args )
{
  ChemCH43DDiffusiveVarSet::configure(args);

  // add here configuration, specific of this class

  if(_vecD.size() != getModel()->getNbEquations()) {
    throw Common::BadValueException (FromHere(),"ChemCH43DPrimDiffusive::configure() : Wrong number of diffusive parameters.");
  }
}

//////////////////////////////////////////////////////////////////////////////

      } // namespace CH4

    } // namespace Chemistry

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
