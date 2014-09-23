#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
 
    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

FVMCC_BC::FVMCC_BC(const std::string& name) :
  CellCenterFVMCom(name),
  _putGhostsOnFace(false),
  socket_normals("normals"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nodes("nodes"),
  m_velocityIDs(),
  m_computeVars(),
  m_t(0.), 
  m_drXiXw(0.), 
  m_drXiXg(0.),
  m_factor(0.),
  m_innerNode(CFNULL), 
  m_faceMidPoint(),
  m_XgXm(),
  m_tempNode(), 
  m_midNode(), 
  m_tempGhostNode(), 
  m_tempGhostNodeBkp(), 
  m_faceNormal(),
  m_fullLoop(true)
{
  addConfigOptionsTo(this);
  
  _zeroGradient = vector<bool>();
  setParameter("ZeroGradientFlags",&_zeroGradient); 
  
  m_coeff = 2.0; 
  setParameter("CoeffMove",&m_coeff); 
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_BC::~FVMCC_BC()
{
}

//////////////////////////////////////////////////////////////////////////////
      
      void FVMCC_BC::defineConfigOptions(Config::OptionList& options)
      {
	options.addConfigOption< vector<bool> >
	  ("ZeroGradientFlags",
	   "Flag specifying constantly extrapolated variables");
  
	options.addConfigOption< CFreal > 
	  ("CoeffMove", "coefficient for the ghost node movement"); 
      }
      
      //////////////////////////////////////////////////////////////////////////////

      void FVMCC_BC::setup()
      {
	// set the IDs corresponding to the velocity components
	getMethodData().getUpdateVar()->setStateVelocityIDs(m_velocityIDs);
  
	const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
	if(_zeroGradient.size() > 0 && _zeroGradient.size() != nbEqs) {
	  CFout << "WARNING: FVMCC_BC " << getName() << " has zero gradient flags not properly set\n";
	  CFout << "WARNING: FVMCC_BC " << getName() << " sets all zero gradient flags=false\n";
	}
  
	if(_zeroGradient.size() == 0) {
	  _zeroGradient.resize(nbEqs);
	  _zeroGradient.assign(nbEqs, false);
	}
  
	CFout << "TRS " << getName() << " has zero gradient flags = ";
	for (CFuint i = 0; i < _zeroGradient.size(); ++i) {
	  CFout << _zeroGradient[i] << " ";
	}
  CFout << "\n";
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  m_faceMidPoint.resize(dim);
  m_XgXm.resize(dim);
  m_tempNode.resize(dim); 
  m_midNode.resize(dim); 
  m_tempGhostNode.resize(dim); 
  m_tempGhostNodeBkp.resize(dim); 
  m_faceNormal.resize(dim); 
  
  m_computeVars.resize(PhysicalModelStack::getActive()->getNbEq(), true);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BC::executeOnTrs()
{
  m_fullLoop = true;
  
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "FVMCC_BC::execute() called for TRS: "
	<< trs->getName() << "\n");

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  geoData.trs = trs;

  preProcess();

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    CFLogDebugMed( "iFace = " << iFace << "\n");

    // build the GeometricEntity
    geoData.idx = iFace;

    GeometricEntity *const face = geoBuilder->buildGE();
    getMethodData().getCurrentFace() = face;
    
    setGhostState(face);

    // release the GeometricEntity
    geoBuilder->releaseGE();
  } 
  
  m_fullLoop = false;
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_BC::computeFlux(RealVector& result)
{
  getMethodData().getFluxSplitter()->computeFlux(result);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
FVMCC_BC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_BC::preProcess()
{
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_BC::computeGhostPosition(Framework::GeometricEntity *const face) 
{  
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  const CFuint faceID = face->getID(); 
  const CFuint startID = faceID*dim; 
  
  DataHandle<CFreal> normals = socket_normals.getDataHandle(); 
  
  // set the current normal 
  for (CFuint i = 0; i < dim; ++i) { 
    m_faceNormal[i] = normals[startID + i]; 
  } 
  
  // compute the original position of the ghost state @see ComputeDummyState 
  const CFuint nbNodesInFace = face->nbNodes();	
  m_faceMidPoint = 0.;
  for (CFuint i =0; i < nbNodesInFace; ++i) {
    m_faceMidPoint += *face->getNode(i);
  }
  m_faceMidPoint /= (CFreal)nbNodesInFace;

  CFreal k = 0.0;
  //if (face->nbNodes() < 4) { 
  // k = - MathFunctions::innerProd(m_faceNormal, *face->getNode(0)); 
  k = - MathFunctions::innerProd(m_faceNormal, m_faceMidPoint);
  // } 
  // else if (face->nbNodes() == 4) { 
  // if the face has 4 nodes, they could nopt lie all on the same plane 
  // so we take an average 
  //  const CFreal p0 = - MathFunctions::innerProd(m_faceNormal, *face->getNode(0)); 
  // const CFreal p1 = - MathFunctions::innerProd(m_faceNormal, *face->getNode(1)); 
  // const CFreal p2 = - MathFunctions::innerProd(m_faceNormal, *face->getNode(2)); 
  // const CFreal p3 = - MathFunctions::innerProd(m_faceNormal, *face->getNode(3)); 
  // k = 0.25*(p0 + p1 + p2 + p3); 
  //  } 
  
  const CFreal n2 = MathFunctions::innerProd(m_faceNormal, m_faceNormal);
  cf_assert (std::abs(n2) > 0.0); 
  State *const innerState = face->getState(0); 
  m_innerNode = &innerState->getCoordinates(); 
  m_t = (MathFunctions::innerProd(m_faceNormal,*m_innerNode) + k)/n2; 
  
  m_tempGhostNode = (*m_innerNode) - 2.*m_t*m_faceNormal; 
  
  /* if (face->nbNodes() == 4) {  
  // we check if the candidate ghost point G is on the other side with respect to the face mid point F  
  // and the inner state I   
  // consider the two vectors (G-F) and (I-F): if their dot product is > 0 they are on the same   
  // side of the face   
  
  static RealVector ghostF(3);  
  static RealVector innerF(3);  
  static RealVector midF(3);  
  
  midF = 0.25*((*face->getNode(0)) + (*face->getNode(1)) + (*face->getNode(2)) + (*face->getNode(3))); 
  ghostF = m_tempGhostNode - midF ;  
  innerF = *m_innerNode - midF;  
  if (MathFunctions::innerProd(ghostF,innerF) > 0.) {   
  //      cout << "NoSlipWallIsothermalNSvt => ghost wrong !! "<< endl; abort();  
  
  // xG = xI + 2 * (xF - xI) = xI - 2 * (xI - xF)   
  m_tempGhostNode = (*m_innerNode) - 2.0*innerF;  
  }  
  } */ 
  
  // this middle node is by construction on the boundary face 
  m_midNode = 0.5*(*m_innerNode + m_tempGhostNode); 
  
  // first calculate the "unmodified distances" inner-wall, inner-ghost 
  m_drXiXw = MathTools::MathFunctions::getDistance(*m_innerNode,m_midNode); 
  m_drXiXg = MathTools::MathFunctions::getDistance(*m_innerNode, m_tempGhostNode); 

} 
      
//////////////////////////////////////////////////////////////////////////////
      
} // namespace FiniteVolume
      
} // namespace Numerics
  
} // namespace COOLFluiD
