#include <fstream>
#include <iostream>

#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/CUDA/CudaDeviceManager.hh"
#include "Common/CUDA/CFVec.hh"
#include "Common/CUDA/CudaTimer.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/SocketBundleSetter.hh"

#include "FiniteVolume/CellCenterFVM.hh"

#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "RadiativeTransfer/Solvers/FiniteVolumeDOM/RadiativeTransferFVDOMCUDA.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferFVDOMCUDA, 
		      DataProcessingData, 
		      RadiativeTransferModule>
radiativeTransferFVDOMCUDAProvider("RadiativeTransferFVDOMCUDA");

//////////////////////////////////////////////////////////////////////////////

__global__ void getFieldOpacitiesKernel(const bool useExponentialMethod,
					const CFuint TID, 
					const CFuint PID,
					const CFuint nbTemp,
					const CFuint nbPress,
					const CFuint nbBins,
					const CFuint ib,
					const CFuint nbEqs,
					const CFuint nbCells,
					const CFreal* states,
					const CFreal* volumes,
					const CFreal* Ttable,
					const CFreal* Ptable,
					const CFreal* opacities,
					const CFreal* radSource,
					CFreal* fieldSource,
					CFreal* fieldAbsor,
					CFreal* fieldAbSrcV,
					CFreal* fieldAbV)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (cellID < nbCells) {
    fieldSource[cellID] = 0.;
    if(useExponentialMethod) {
      fieldAbsor[cellID]  = 0.;
    }
    else {
      fieldAbSrcV[cellID] = 0.;
      fieldAbV[cellID]    = 0.;
    }
    
    //Get the field pressure and T commented because now we impose a temperature profile
    const CFuint sIdx = cellID*nbEqs; 
    const CFreal p = states[sIdx + PID];
    const CFreal T = states[sIdx + TID];
    const CFreal patm = p/101325.; //converting from Pa to atm
    
    CFreal val1 = 0;
    CFreal val2 = 0;
    
    RadiativeTransferFVDOM::Interpolator interp;
    interp.tableInterpolate(nbBins, nbTemp, nbPress, Ttable, Ptable,
			    opacities, radSource, T, patm, ib, val1, val2); 
    
    if(useExponentialMethod){
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	fieldSource[cellID] = 1e-30;
	fieldAbsor[cellID]  = 1e-30;
      }
      else {
	fieldSource[cellID] = val2/val1;
	fieldAbsor[cellID]  = val1;
      }
    }
    else{
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	fieldSource[cellID] = 1e-30;
	fieldAbV[cellID]    = 1e-30*volumes[cellID]; // Volumen converted from m^3 into cm^3
      }
      else {
	fieldSource[cellID] = val2/val1;
	fieldAbV[cellID]    = val1*volumes[cellID];
      }      
      fieldAbSrcV[cellID]   = fieldSource[cellID]*fieldAbV[cellID];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelExponential(const CFuint nbCells,
					  const CFreal weightIn,
					  const CFuint* cellFaces,
					  const CFint* isOutward,
					  const CFint* advanceOrder,
					  const CFreal* volumes,
					  const CFreal* fieldSource,
					  const CFreal* fieldAbsor,
					  const CFreal* dotProdInFace,
					  CFreal* In,
					  CFreal* II,
					  CFreal* divQ)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  __shared__ CFreal weight;
  weight = weightIn;
      
  if (cellID < nbCells) {
    // allocate the cell entity
    const CFuint iCell   = abs(advanceOrder[cellID]);
    CFreal Ic            = 0.;
    CFreal inDirDotnANeg = 0.;
    CFreal dirDotnANeg   = 0;
    CFreal Lc            = 0;
    CFreal halfExp       = 0;
    
    const CFuint nbFaces = 5; //// cellFaces->nbCols(iCell);  /////
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) { 
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      
      if(dirDotNA < 0.) {
	dirDotnANeg += dirDotNA;
        // const CFint fcellID = faceNeighborID[faceID*2]; 
        // const CFint neighborCellID = (fcellID == iCell) ? faceNeighborID[faceID*2+1] : fcellID;
	// const CFreal source = (neighborCellID >=0) ? In[neighborCellID] : fieldSource[iCell];
        // inDirDotnANeg +=source*dirDotNA;
	
	
	/*const bool isBFace = mapGeoToTrs->isBGeo(faceID); /////
	  if (!isBFace){
	  const CFuint neighborCellID = getNeighborCellID(faceID, iCell); /////
	  inDirDotnANeg += In[neighborCellID]*dirDotNA;
	  }
	  else {
	  const CFreal boundarySource = fieldSource[iCell];
	  inDirDotnANeg += boundarySource*dirDotNA;
	  }*/
      }
    } 
    Lc          = volumes[iCell]/(- dirDotnANeg); 
    halfExp     = std::exp(-0.5*Lc*fieldAbsor[iCell]);
    In[iCell] = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*fieldSource[iCell];
    Ic          = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*fieldSource[iCell];
    
    // m_q(iCell,XX) += Ic*m_dirs(d,0)*weight; /////
    // m_q(iCell,YY) += Ic*m_dirs(d,1)*weight; /////
    // m_q(iCell,ZZ) += Ic*m_dirs(d,2)*weight; /////
    
    CFreal inDirDotnA = inDirDotnANeg;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = cellFaces[iFace*nbCells + iCell];
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      if (dirDotNA > 0.) {
	inDirDotnA += In[iCell]*dirDotNA;
      }
    }
    
    divQ[iCell] += inDirDotnA*weight;
    II[iCell]   += Ic*weight;
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

__global__ void computeQKernelNoExponential(const CFuint nbCells,
					    const CFreal weight,
					    const CFuint* cellFaces,
					    const CFint* isOutward,
					    const CFint* advanceOrder,
					    const CFreal* volumes,
					    const CFreal* fieldSource,
					    const CFreal* fieldAbSrcV,
					    const CFreal* fieldAbV,
					    const CFreal* dotProdInFace,
					    CFreal* In,
					    CFreal* II,
					    CFreal* divQ)
{    
  // each thread takes care of computing the gradient for one single cell
  const int cellID = threadIdx.x + blockIdx.x*blockDim.x;
  
  // __shared__ typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE> s_dcol[32];
  // typename LIMITER::BASE::template DeviceConfigOptions<NOTYPE>* dcol = &s_dcol[threadIdx.x];
  // dcol->init(gdcol);
  
  // __shared__ typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE> s_dcor[32];
  // typename POLYREC::BASE::template DeviceConfigOptions<NOTYPE>* dcor = &s_dcor[threadIdx.x];
  // dcor->init(gdcor);
        
  if (cellID < nbCells) {
    /*  CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    
    // allocate the cell entity
    const CFuint iCell = std::abs(m_advanceOrder[cellID]);
      const CFuint nbFaces = cellFaces->nbCols(iCell);
    
    
    
      CFreal dirDotnAPos = 0;
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = (*cellFaces)(iCell, iFace);
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFreal dirDotNA = dotProdInFace[faceID]*factor;
	
	if (dirDotNA >= 0.){
	  dirDotnAPos += dirDotNA;
	}
	else {
	  const bool isBFace = m_mapGeoToTrs->isBGeo(faceID);
	  if (!isBFace){
	    const CFuint neighborCellID = getNeighborCellID(faceID, iCell);
	    inDirDotnANeg += m_In[neighborCellID]*dirDotNA;
	  }
	  else {
	    const CFreal boundarySource = m_fieldSource[iCell];
	    inDirDotnANeg += boundarySource*dirDotNA;
	  }
	}
      } 
      m_In[iCell] = (m_fieldAbSrcV[iCell] - inDirDotnANeg)/(m_fieldAbV[iCell] + dirDotnAPos);
      Ic = m_In[iCell];
    
    
    m_q(iCell,XX) += Ic*m_dirs(d,0)*m_weight[d];
    m_q(iCell,YY) += Ic*m_dirs(d,1)*m_weight[d];
    m_q(iCell,ZZ) += Ic*m_dirs(d,2)*m_weight[d];
    
    CFreal inDirDotnA = inDirDotnANeg;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint faceID = (*cellFaces)(iCell, iFace);
      const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
      const CFreal dirDotNA = dotProdInFace[faceID]*factor;
      if (dirDotNA > 0.) {
	inDirDotnA += m_In[iCell]*dirDotNA;
      }
    }
    
    divQ[iCell] += inDirDotnA*m_weight[d];
    m_II[iCell]   += Ic*m_weight[d];*/
  }  
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOMCUDA::RadiativeTransferFVDOMCUDA(const std::string& name) :
  RadiativeTransferFVDOM(name)
{
  addConfigOptionsTo(this);
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferFVDOMCUDA::~RadiativeTransferFVDOMCUDA()
{
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::defineConfigOptions(Config::OptionList& options)
{  
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::setup()
{
  CFAUTOTRACE;
  
  RadiativeTransferFVDOM::setup();
  
  // store invariant data on GPU
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  cellFaces->getPtr()->put();
  socket_isOutward.getDataHandle().getLocalArray()->put(); 
  socket_volumes.getDataHandle().getLocalArray()->put(); 
  socket_divq.getDataHandle().getLocalArray()->put(); 
  
  m_fieldSource.put();
  m_fieldAbsor.put();
  m_fieldAbSrcV.put();
  m_fieldAbV.put();  
  m_In.put();
  m_II.put();
  m_opacities.put();
  m_radSource.put();
  m_Ttable.put();
  m_Ptable.put();
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferFVDOMCUDA::unsetup()
{
  CFAUTOTRACE;
  
  RadiativeTransferFVDOM::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////
 
void RadiativeTransferFVDOMCUDA::loopOverDirs(const CFuint startBin, 
					      const CFuint endBin, 
					      const CFuint startDir,
					      const CFuint endDir)
{
  SafePtr<ConnectivityTable<CFuint> > cellFaces = 
    MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbCells = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  states.getGlobalArray()->put();
  
  const CFuint blocksPerGrid = 
    CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells);
  const CFuint nThreads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
  
  for (CFuint d = startDir; d < endDir; ++d) {
    CFLog(INFO, "( dir: " << d << " ), ( bin: ");
    const CFuint bStart = (d != startDir) ? 0 : startBin;
    const CFuint bEnd   = (d != m_startEndDir.second) ? m_nbBins : endBin;
    
    m_advanceOrder[d].put();
    
    // precompute dot products for all faces and directions (a part from the sign)
    computeDotProdInFace(d, m_dotProdInFace);
    m_dotProdInFace.put();
    
    for (CFuint ib = startBin; ib < endBin; ++ib) {
      CFLog(INFO, "[dir, bin] = [" << d << ", " << ib << "]\n");
      
      // precompute the radiation properties for all cells
      getFieldOpacitiesKernel<<<blocksPerGrid,nThreads>>>
	(m_useExponentialMethod, 
	 m_TID, m_PID, m_nbTemp, m_nbPress, m_nbBins,
	 ib, nbEqs, nbCells, m_Ttable.ptrDev(), m_Ptable.ptrDev(), 
	 states.getGlobalArray()->ptrDev(),
	 volumes.getLocalArray()->ptrDev(),
	 m_opacities.ptrDev(),
	 m_radSource.ptrDev(),
	 m_fieldSource.ptrDev(),
	 m_fieldAbsor.ptrDev(),
	 m_fieldAbSrcV.ptrDev(),
	 m_fieldAbV.ptrDev());  
      
      // m_fieldSource.get();
      //  m_fieldAbsor.get();
      // m_fieldAbSrcV.get();
      // m_fieldAbV.get();
      
      // RadiativeTransferFVDOM::computeQ(ib,d);
      
      // compute the radiative heat flux
      if (m_useExponentialMethod) {
	computeQKernelExponential<<<blocksPerGrid,nThreads>>> 
	  (nbCells, m_weight[d],
	   cellFaces->getPtr()->ptrDev(),
	   isOutward.getLocalArray()->ptrDev(),
	   m_advanceOrder[d].ptrDev(),
	   volumes.getLocalArray()->ptrDev(),
	   m_fieldSource.ptrDev(),
	   m_fieldAbsor.ptrDev(),
	   m_dotProdInFace.ptrDev(),
	   m_In.ptrDev(), m_II.ptrDev(), divQ.getLocalArray()->ptrDev());
      }
      else {
	computeQKernelNoExponential<<<blocksPerGrid,nThreads>>> 
	  (nbCells, m_weight[d],
	   cellFaces->getPtr()->ptrDev(),
	   isOutward.getLocalArray()->ptrDev(),
	   m_advanceOrder[d].ptrDev(),
	   volumes.getLocalArray()->ptrDev(),
	   m_fieldSource.ptrDev(),
	   m_fieldAbSrcV.ptrDev(),
	   m_fieldAbV.ptrDev(),
	   m_dotProdInFace.ptrDev(),
	   m_In.ptrDev(), m_II.ptrDev(), divQ.getLocalArray()->ptrDev());
	
      }
      CFLog(INFO, ")\n");
    }
  }
  socket_divq.getDataHandle().getLocalArray()->get(); 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

