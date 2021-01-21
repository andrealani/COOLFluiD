#include "FVMCC_PolyRec.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption<CFuint, Config::DynamicOption<> >
    ("StopLimiting","Stop applying the limiter.");
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_PolyRec::FVMCC_PolyRec(const std::string& name) :
  PolyReconstructor<CellCenterFVMData>(name),
  socket_cellFlag("cellFlag"),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_normals("normals"),
  _isLimiterNull(false),
  _quadPointCoord(),
  _tmpLimiter(), 
  _gradientCoeff(),
  _vFunction()
{
  addConfigOptionsTo(this);
  
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);
  
  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
  
  _stopLimiting = 0;
  setParameter("StopLimiting",&_stopLimiting);
  
  // fix high default value   
  _limitIter = 1000000000;
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_PolyRec::~FVMCC_PolyRec()
{
  for (CFuint iFace = 0; iFace < _quadPointCoord.size(); ++iFace) {
    deletePtr(_quadPointCoord[iFace][0]);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::updateWeights()
{
  //nothing to do here
}


//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::setup()
{
  PolyReconstructor<CellCenterFVMData>::setup();
  
  _isLimiterNull = getMethodData().getLimiter()->isNull();
  _tmpLimiter.resize(PhysicalModelStack::getActive()->getNbEq());
  _gradientCoeff.resize(PhysicalModelStack::getActive()->getNbEq());
  
  // if the limiter is null set all its entries to 1.0
  if (_isLimiterNull) { 
    socket_limiter.getDataHandle() = 1.0;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::unsetup()
{
  for (CFuint iFace = 0; iFace < _quadPointCoord.size(); ++iFace) {
    deletePtr(_quadPointCoord[iFace][0]);
  }
  SwapEmpty(_quadPointCoord);
  
  PolyReconstructor<CellCenterFVMData>::unsetup();
}

//////////////////////////////////////////////////////////////////////////////
      
void FVMCC_PolyRec::prepareReconstruction()
{
  // special treatment if input vars is "i"  (= iteration number)
  if (_vars.size() == 1 && _vars[0] == "i") {
    RealVector input(1); input[0] = SubSystemStatusStack::getActive()->getNbIter();
    _vFunction.evaluate(input, _gradientCoeff);
    _gradientFactor = _gradientCoeff[0];
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::computeFaceLimiter(GeometricEntity* const face)
{
  if (!_isLimiterNull) { 
    DataHandle<CFreal> newLimiter = socket_limiter.getDataHandle();
    SafePtr<Limiter<CellCenterFVMData> > limiter = getMethodData().getLimiter();
    DataHandle<bool> cellFlag = socket_cellFlag.getDataHandle();
    
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
    const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
 
    for (CFuint iCell = 0; iCell < 2; ++iCell) {
      if (face->getState(iCell)->isGhost()) continue;
      
      const CFuint cellID = face->getState(iCell)->getLocalID();
      
      if (!cellFlag[cellID]) {
	GeometricEntity *const currCell = face->getNeighborGeo(iCell);
	const GeomEntList *const faces = currCell->getNeighborGeos();
	const CFuint nbFaces = faces->size();
	for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	  Node& faceMidCoord = *_quadPointCoord[iFace][0];	
	  computeMidPoint(*(*faces)[iFace]->getNodes(), faceMidCoord);
	}

	if (residual > _limitRes && (_limitIter > 0 && iter < _limitIter)) {	
	  const CFuint stateID = currCell->getState(0)->getLocalID();
	  const CFuint startID = stateID*nbEqs;
	  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
	    newLimiter[startID + iEq] = 1.0;
	  }
	  limiter->limit(_quadPointCoord, currCell, &newLimiter[startID]);
	}
	else {
	  if (!_freezeLimiter) {
	    // historical modification of the limiter
	    limiter->limit(_quadPointCoord, currCell, &_tmpLimiter[0]);
	    const CFuint stateID = currCell->getState(0)->getLocalID();
	    CFuint currID = stateID*nbEqs;
	    for (CFuint iVar = 0; iVar < nbEqs; ++iVar, ++currID) {
	      newLimiter[currID] = min(_tmpLimiter[iVar],newLimiter[currID]);
	    }
	  }
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::baseExtrapolateImpl(GeometricEntity* const face)
{ 
  //this is only usable if one quadrature point per face is needed
  //compute the position of the face centroid (quadrature point)
  
  Node& faceMidCoord = *_extrapCoord[0];	
  computeMidPoint(*face->getNodes(), faceMidCoord);
  
  getValues(0).setSpaceCoordinates(&faceMidCoord);
  getValues(1).setSpaceCoordinates(&faceMidCoord);
  
  cf_assert(!getValues(0).getCoordinates().isOwnedByState());
  cf_assert(!getValues(1).getCoordinates().isOwnedByState());
  
  getValues(0).setLocalID(face->getState(LEFT)->getLocalID());
  getValues(1).setLocalID(face->getState(RIGHT)->getLocalID());
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::allocateReconstructionData()
{
  // this data allocation is ok for 1st and 2nd order FV method 
  // one quadrature point position 
  _extrapCoord.resize(1);
  _extrapCoord[0] = new Node(true);
    
  _extrapValues.resize(2);
  for (CFuint i = 0; i < _extrapValues.size(); ++i) {
    _extrapValues[i] = new State();
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _backupValues.resize(2);
  for (CFuint i = 0; i < _backupValues.size(); ++i) {
    _backupValues[i].resize(nbEqs);
  }
  
  _extrapPhysData.resize(2); 
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm(); 
  for (CFuint i = 0; i < _extrapPhysData.size(); ++i) {
    convTerm->resizePhysicalData(_extrapPhysData[i]);
  } 
  
  _backupPhysData.resize(2); 
  for (CFuint i = 0; i < _backupPhysData.size(); ++i) {
    convTerm->resizePhysicalData(_backupPhysData[i]);
  } 
  
  // one quadrature point per face /// @TODO this fails for curved boundaries ...
  _quadPointCoord.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint iFace = 0; iFace < _quadPointCoord.size(); ++iFace) {
    _quadPointCoord[iFace].resize(1);
    _quadPointCoord[iFace][0] = new Node();
  }
  
  if (_functions.size() > 0) {
    _vFunction.setFunctions(_functions);
    _vFunction.setVariables(_vars);
    try {
      _vFunction.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
FVMCC_PolyRec::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
    Framework::PolyReconstructor<CellCenterFVMData>::needsSockets();
  
  result.push_back(&socket_cellFlag);
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_normals);
  
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_PolyRec::configure ( Config::ConfigArgs& args )
{
  Framework::PolyReconstructor<CellCenterFVMData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////
 
const RealVector& FVMCC_PolyRec::getCurrentGeoShapeFunction(GeometricEntity *const face)
{
  Common::SafePtr<Framework::VolumeIntegrator> volInt = getMethodData().getVolumeIntegrator();
  const vector<RealVector>& shapeFunctions = volInt->getGeometryIntegrator(face)->computeShapeFunctionsAtQuadraturePoints();
  cf_assert(shapeFunctions.size() >=  1);
  return shapeFunctions[0];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
