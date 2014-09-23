#include "MutationI/OnlySigma.hh"
#include "MutationI/Mutation.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Common/OSystem.hh"
#include "Common/BadValueException.hh"
#include "Environment/DirPaths.hh"
#include "Common/Stopwatch.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Mutation {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<OnlySigma,
			    PhysicalPropertyLibrary,
			    MutationModule,
			    1>
onlySigmaLibraryProvider("OnlySigma");

//////////////////////////////////////////////////////////////////////////////

void OnlySigma::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("mixtureName","Name of the mixture.");
}

//////////////////////////////////////////////////////////////////////////////

OnlySigma::OnlySigma(const std::string& name)
  : PhysicalChemicalLibrary(name
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);
 }

//////////////////////////////////////////////////////////////////////////////

OnlySigma::~OnlySigma()
{
}

//////////////////////////////////////////////////////////////////////////////

void OnlySigma::configure ( Config::ConfigArgs& args )
{
  PhysicalChemicalLibrary::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

CFdouble OnlySigma::sigma(CFdouble& temp,   //electrical conductivity
			  CFdouble& pressure,
			  CFreal* tVec)
{
  //composition must be called before!
  CFdouble sigma = 0.0;
  // here implement
  return sigma;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Mutation

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

