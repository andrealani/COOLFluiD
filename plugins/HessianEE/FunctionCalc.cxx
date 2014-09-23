#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ProxyDofIterator.hh"

#include "HessianEE/HessianEE.hh"
#include "HessianEE/FunctionCalc.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FunctionCalc, HessEEData, HessianEEModule> FunctionCalcProvider("FunctionCalc");

//////////////////////////////////////////////////////////////////////////////

void FunctionCalc::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("RefH","Reference H.");
}

//////////////////////////////////////////////////////////////////////////////

FunctionCalc::FunctionCalc(const std::string& name) :
  HessEECom(name),
  socket_states("states"),
  socket_adapt_func("adapt_func"),
  socket_nstatesProxy("nstatesProxy")
{
   addConfigOptionsTo(this);
  _refH = 10.0;
   setParameter("RefH",&_refH);
}

//////////////////////////////////////////////////////////////////////////////

FunctionCalc::~FunctionCalc()
{
}

//////////////////////////////////////////////////////////////////////////////

void FunctionCalc::setup()
{
  CFAUTOTRACE;

  // first call parent method
  HessEECom::setup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
FunctionCalc::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_adapt_func);
  result.push_back(&socket_nstatesProxy);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FunctionCalc::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> adapt_func = socket_adapt_func.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy = socket_nstatesProxy.getDataHandle();

  // this isa sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  cf_assert( adapt_func.size() == nodalStates.getSize());

  for ( CFuint i = 0; i < nodalStates.getSize(); ++i)
  {
    const RealVector& currState = *nodalStates.getState(i);
    adapt_func[i] = currState[0];

    //cout << currState[0] << endl;

    //adapt_func[i] = (*states[i])[0];
  }

#if 0
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  for(CFuint iTR = 0; iTR < trs->getNbTRs(); ++iTR)
  {
    TopologicalRegion* tr = trs->getTopologicalRegion(iTR);
    for(CFuint iGeoEnt = 0; iGeoEnt < tr->getNbGeomEnt(); ++iGeoEnt)
    {
      GeometricEntity* cell = tr->getGeomEnt(iGeoEnt);
      vector<State*>* cellStates = cell->getStates();

      CFLogDebugMax("Cell " << iGeoEnt << "\n");

      /// @todo Implement here!!!!

      CFuint nbStatesInCell = cell->getNbNodesSolutionShapeFunction();
      for ( CFuint iState = 0; iState < nbStatesInCell; ++iState )
      {
      }
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
