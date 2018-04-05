#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/NeumannBC.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBC,
                      FluctuationSplitData,
                      FluctSplitModule>
neumannBCProvider("NeumannBC");

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::defineConfigOptions(Config::OptionList& options)
{
  //  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::NeumannBC(const std::string& name) :
  SuperInlet(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell")
{
  addConfigOptionsTo(this);
  
  // m_functions = std::vector<std::string>();
  // setParameter("Def",&m_functions);
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::~NeumannBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setup()
{
  CFAUTOTRACE;

  SuperInlet::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::unsetup()
{
  CFAUTOTRACE;
  
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "NeumannBC::execute() called for TRS: " << trs->getName() << "\n");
  
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle< State*, GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector variables(dim);
  bool applyBC = true;
  
  SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
    CFLogDebugMed( "NeumannBC::executeOnTrs() => iFace = " << iFace << "\n");
    
    // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity& currFace = *geoBuilder->buildGE();
    const CFreal faceAreaOvDim = getFaceArea(currFace.getID())/(CFreal)dim;
    cf_assert(faceAreaOvDim > 0.);
    const CFuint nbStatesInFace = currFace.nbStates();
    for (CFuint iState = 0 ; iState < nbStatesInFace; ++iState) {
      const CFuint stateID = currFace.getState(iState)->getLocalID();
      State *const state = states[stateID];
      
      // find if we should apply the BC to this node
      if (m_checkCondition) {
	const CFreal applyBCvalue = m_condition.Eval(&state->getCoordinates()[0]);
	applyBC = ((!isUpdated[stateID]) && (applyBCvalue > 0.0));
      }
      else {
	applyBC = (!isUpdated[stateID]);
      }
      
      // apply the BC to this node
      if (applyBC) {
	// CFLog(INFO, "NeumannBC::executeOnTrs() => applying BC\n");
	//Set the values of the variables xyz + time
	for (CFuint iDim = 0; iDim < dim; ++iDim){
	  variables[iDim] = state->getCoordinates()[iDim];
	}
	
	//Evaluate the function
	m_vFunction.evaluate(variables, *_input);
	
	for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	  // this sign has to be checked
	  rhs[stateID] += (*_input)[iEq]*faceAreaOvDim;
	}
	
	// isUpdated[stateID] = true; 
      }
    }
    
    //Release the geometric entity
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal NeumannBC::getFaceArea(const CFuint faceID)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Trio<CFuint, CFuint, SafePtr<TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);
  return normals[cellLocalID]->getAreaFace(iFaceLocal); 
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > NeumannBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = SuperInlet::needsSockets();
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit
  
} // namespace COOLFluiD
