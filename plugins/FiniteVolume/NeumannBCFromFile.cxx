#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/NeumannBCFromFile.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBCFromFile, CellCenterFVMData, FiniteVolumeModule> 
neumannBCFromFileFVMCCProvider("NeumannBCFromFileFVMCC");
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBCFromFile::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VarIDs","IDs of the variables from which values are read by file");
}
      
//////////////////////////////////////////////////////////////////////////////

NeumannBCFromFile::NeumannBCFromFile(const std::string& name) :
  SuperInlet(name)
{
  addConfigOptionsTo(this);
  
  m_varIDs = vector<CFuint>();
  setParameter("VarIDs",&m_varIDs);
}
      
//////////////////////////////////////////////////////////////////////////////

NeumannBCFromFile::~NeumannBCFromFile()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCFromFile::setup()
{
  SuperInlet::setup();
  
  if (m_varIDs.size() == 0) {m_varIDs.resize(1, 0);}
  
  m_mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCFromFile::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBCFromFile::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  cf_assert(getMethodData().getNodalStatesExtrapolator()->getMapTrs2NodalValues()->size() > 0);
  
  vector<Node*>& nodesInFace = *face->getNodes();
  const CFuint nbNodesInFace = nodesInFace.size();
  SafePtr<TopologicalRegionSet> trs = m_mapGeoToTrs->getTrs(face->getID());
  cf_assert(trs.isNotNull());
  
  // build the mapTrs2NodalValues storage
  SafePtr<NodalStatesExtrapolator<CellCenterFVMData> > nse =
    this->getMethodData().getNodalStatesExtrapolator();
  SafePtr<vector<NodalStatesExtrapolator<CellCenterFVMData>::MapTrs2NodalValues*> >
    mapTrs2NodalValues = nse->getMapTrs2NodalValues();
  
  cf_assert(m_varIDs.size() > 0);
  
  for (CFuint iVar = 0; iVar < m_varIDs.size(); ++iVar) {
    const CFuint varID = m_varIDs[iVar];
    cf_assert(varID < mapTrs2NodalValues->size());

    RealVector& bArray = *(*mapTrs2NodalValues)[varID]->find(&*trs);
    CFMap<CFuint,CFuint>* mapNodeIDs = nse->getMapTrs2NodeIDs()->find(&*trs);
    
    // compute the temperature in the face centroid as the average of the nodal values
    cf_assert(varID < _input->size());
    
    (*_input)[varID] = 0.0;
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      const CFuint localNodeID = nodesInFace[iNode]->getLocalID();
      (*_input)[varID] += bArray[mapNodeIDs->find(localNodeID)];
    }
    (*_input)[varID] /= nbNodesInFace;
  }
  
  // we assume that (U_i - U_g)/dr = f(x,y,z) => U_g = U_i - f*dr
  const CFreal dr = MathTools::MathFunctions::getDistance
    (ghostState->getCoordinates(), innerState->getCoordinates());
  
  CFLog(DEBUG_MED, "NeumannBCFromFile::setGhostState() => (*_input) = " << *_input << ", dr = " << dr  << "\n");
  *ghostState = *innerState - (*_input)*dr;
  
  CFLog(DEBUG_MED, "NeumannBCFromFile::setGhostState() => ghostState = " << *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
