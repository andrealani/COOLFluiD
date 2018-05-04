#include "FiniteVolume/FiniteVolume.hh"
#include "LeastSquareP1Setup.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/CFMultiMap.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteVolume/ComputeStencil.hh"
#include "Framework/ComputeNormals.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LeastSquareP1Setup, CellCenterFVMData, FiniteVolumeModule> 
leastSquareP1SetupProvider("LeastSquareP1Setup");

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1Setup::LeastSquareP1Setup(const std::string& name) :
  StdSetup(name),
  socket_stencil("stencil"),
  socket_weights("weights"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_uZ("uZ")
{
}

//////////////////////////////////////////////////////////////////////////////

LeastSquareP1Setup::~LeastSquareP1Setup()
{
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1Setup::configure ( Config::ConfigArgs& args )
{
  StdSetup::configure(args);

  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<MeshData> meshData = MeshDataStack::getInstance().getEntryByNamespace(nsp);
  
  // set the number of overlap layers to two for parallel runs
  if (PE::GetPE().IsParallel()) {
    meshData->setNbOverlapLayers(2);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LeastSquareP1Setup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = StdSetup::providesSockets();

  result.push_back(&socket_stencil);
  result.push_back(&socket_weights);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1Setup::execute()
{
  CFAUTOTRACE;

  StdSetup::execute();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  stencil.resize(nbStates);

  // ghost states are FOR THE MOMENT excluded from the stencil
  computeStencil();

  const CFuint nbEdges = countEdges();
  CFLog(VERBOSE,"LeasSquareP1 : Detected " << nbEdges << " edges\n");

  DataHandle<CFreal> weights = socket_weights.getDataHandle();
  weights.resize(nbEdges);
  
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  uX.resize(nbStates*nbEqs);
  
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  uY.resize(nbStates*nbEqs);

  CFLog(VERBOSE, "LeastSquareP1Setup::execute() => (ux.size(), nbEqs, nbState) = (" <<
	uX.size() << ", " << nbEqs  << ", " << nbStates << ")\n");
  
  if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
    DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
    uZ.resize(nbStates*nbEqs);
  }
  else {
    DataHandle<CFreal> uZ = socket_uZ.getDataHandle();
    uZ.resize(1); // AL: minimum size > 0 to avoid problems with NULL pointers
  }
  
  // resize the limiter storage
  // by default would have size == 0
  const CFuint sizeLimiter = nbStates*nbEqs;
  DataHandle<CFreal> limiter = socket_limiter.getDataHandle();
  limiter.resize(sizeLimiter);
  
  CFLog(VERBOSE, "_limiterSocketName " << _limiterSocketName << "\n");
  
  if (m_dynamicSockets.sinkSocketExists(_limiterSocketName)) {
    DataHandle<CFreal> initLimiter = m_dynamicSockets.
      getSocketSink<CFreal>(_limiterSocketName)->getDataHandle();
    
    cf_assert(initLimiter.size() == limiter.size());
    for (CFuint i = 0; i < sizeLimiter; ++i) {
      limiter[i] = initLimiter[i];
    }
  } 
}

//////////////////////////////////////////////////////////////////////////////

void LeastSquareP1Setup::computeStencil()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  // table that stores the connectivity state-neighbor states
  vector< vector<State*> > neighborStates(nbStates);
  
  SelfRegistPtr<ComputeStencil> computeStencil
    (FACTORY_GET_PROVIDER(getFactoryRegistry(), ComputeStencil, _stencilType)->create(_stencilType));
  configureNested ( computeStencil.getPtr(), m_stored_args );
  
  computeStencil->setDataSocketSinks(socket_states, socket_nodes, socket_stencil, socket_gstates);
  (*computeStencil)();
}

//////////////////////////////////////////////////////////////////////////////

CFuint LeastSquareP1Setup::countEdges()
{
  CFAUTOTRACE;

  // storage for the stencil (pointers to neighbors)
  DataHandle<vector<State*> > stencil = socket_stencil.getDataHandle();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  const CFuint nbStates = states.size();
  CFuint nbEdges = 0;
  
  for(CFuint iState = 0; iState < nbStates; ++iState) {
    const State* const first = states[iState];
    const CFuint stencilSize = stencil[iState].size();
    // loop over the neighbor cells belonging to the chosen stencil
    for(CFuint in = 0; in < stencilSize; ++in) {
      const State* const last = stencil[iState][in];
      const CFuint firstID = first->getLocalID();
      if (!last->isGhost()) {
	const CFuint lastID = last->getLocalID();
	cf_assert(firstID != lastID);
	
	if (lastID > firstID) {
	  ++nbEdges;
	}
      }
      else {
	// you cannot ask the localID to a ghost state
	++nbEdges;
      }
    }
  }
  
  return nbEdges;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
