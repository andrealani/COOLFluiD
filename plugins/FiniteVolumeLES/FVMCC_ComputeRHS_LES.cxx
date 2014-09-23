#include "FVMCC_ComputeRHS_LES.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FVMCC_BC.hh"
#include "LES/LESVarSet.hh"
#include "FiniteVolumeLES/FiniteVolumeLES.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FVMCC_ComputeRHS_LES, CellCenterFVMData, FiniteVolumeLESModule>
FVMCC_ComputeRHS_LESProvider("FVMCC_LES");

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHS_LES::FVMCC_ComputeRHS_LES(const std::string& name) :
  FVMCC_ComputeRHS(name)
{
  CFLog(INFO, "FVMCC_ComputeRHS_LES() \n");
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

FVMCC_ComputeRHS_LES::~FVMCC_ComputeRHS_LES()
{
  
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS_LES::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS_LES::configure ( Config::ConfigArgs& args )
{
  CFLog(INFO, "FVMCC_ComputeRHS_LES::configure() \n");
  FVMCC_ComputeRHS::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS_LES::setup()
{
  CFAUTOTRACE;
  FVMCC_ComputeRHS::setup();
  
  m_lesVarSet = _diffVar.d_castTo<LES::LESVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS_LES::unsetup()
{
  FVMCC_ComputeRHS::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_ComputeRHS_LES::execute()
{
  CFTRACEBEGIN;
 
  // get the volumes
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
 
 
  initializeComputationRHS();
  
  // set the list of faces
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();

  _faceIdx = 0;
  
  // no variable perturbation is needed in explicit residual computation
  getMethodData().setIsPerturb(false);
  
  // prepare the building of the faces
  Common::SafePtr<GeometricEntityPool<FaceCellTrsGeoBuilder> > geoBuilder = getMethodData().getFaceCellTrsGeoBuilder();
  geoBuilder->getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceCellTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  vector<bool> zeroGrad(PhysicalModelStack::getActive()->getNbEq(), false);
  const bool hasSourceTerm = (getMethodData().isAxisymmetric() || getMethodData().hasSourceTerm());
  
  // this could be set during set up with no guarantee that it will be effective:
  // a MethodStrategy could set it to a different value afterwards, before entering here
  geoData.allCells = getMethodData().getBuildAllCells();
 
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
    SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];

    CFLog(VERBOSE, "TRS name = " << currTrs->getName() << "\n");

    // the faces on the boundary of the partition don't have to
    // be processed (their fluxes could give NaN)
    if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {
      if (currTrs->hasTag("writable")) {
        _currBC = _bcMap.find(iTRS);

        CFLog(VERBOSE, "BC name = " << _currBC->getName() << "\n");

        geoData.isBFace = true;

	// set the flags specifying the variables for which the boundary condition
	// imposes constant extrapolation (zero gradient)
	_polyRec->setZeroGradient(_currBC->getZeroGradientsFlags());
      }
      else {
        geoData.isBFace = false;
	_polyRec->setZeroGradient(&zeroGrad);
      }

      // set the current TRS in the geoData
      geoData.faces = currTrs;

        
      const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace, ++_faceIdx) {
        CFLogDebugMed( "iFace = " << iFace << "\n");
	

    	// reset the equation subsystem descriptor
	PhysicalModelStack::getActive()->resetEquationSubSysDescriptor();

        // build the GeometricEntity
        geoData.idx = iFace;
        _currFace = geoBuilder->buildGE();
	
	//	if (_currFace->getState(0)->isParUpdatable() || _currFace->getState(1)->isParUpdatable()) {
	  
	  // set the data for the FaceIntegrator
	  setFaceIntegratorData();
	  
	  // extrapolate (and LIMIT, if the reconstruction is linear or more)
	  // the solution in the quadrature points
	  _polyRec->extrapolate(_currFace);
	
	  // compute the physical data for each left and right reconstructed
	  // state and in the left and right cell centers
	  computePhysicalData();
	  
	  // a jacobian free method requires to re-compute the update coefficient every time the 
	  // residual is calculated to get F*v from the finite difference formula
	  // in particular the time dependent part of the residual depend on a updateCoeff
	  // that has to be up-to-date
	  getMethodData().setIsPerturb(false);
	  
	  const bool isBFace = _currFace->getState(1)->isGhost();
	  if (!isBFace) {
	    _fluxSplitter->computeFlux(_flux);
	  }
	  else {
	    _currBC->computeFlux(_flux);
	  }
 
	  computeInterConvDiff();
	  
	  if (_hasDiffusiveTerm) {
	    // reset to false the flag telling to freeze the diffusive coefficients
	    _diffVar->setFreezeCoeff(false);
	    
      // set the volume in the LESVarSet as the average of the cells 
      // on both sides of the face
      CFreal volume;
      volume  = volumes[_currFace->getState(0)->getLocalID()];
      volume += volumes[_currFace->getState(1)->getLocalID()];
      volume *= 0.5;
      m_lesVarSet->setVolume(volume);
          
	    //_nodalExtrapolator->extrapolateInNodes(*_currFace->getNodes());
	    _diffusiveFlux->computeFlux(_dFlux);
	    
	    _flux -= _dFlux;
	  }
	  
	  CFLogDebugMed("flux = " <<  _flux  << "\n");

	  // compute the source term
	  if (hasSourceTerm) {			
	    computeSourceTerm();
	  }
	  
	  // compute the contribution to the RHS
	  updateRHS();
	  // source term jacobians are only computed while processing internal faces 
	  computeRHSJacobian();
	  //	}
	
	geoBuilder->releaseGE();
      }
    }
  }
  
  finalizeComputationRHS();
 
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////
      
vector<SafePtr<BaseDataSocketSink> > FVMCC_ComputeRHS_LES::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = FVMCC_ComputeRHS::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////
 
 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
