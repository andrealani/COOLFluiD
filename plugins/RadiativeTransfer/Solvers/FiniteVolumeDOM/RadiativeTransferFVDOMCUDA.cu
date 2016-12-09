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

__device__ void tableInterpolate
(const CFuint nbBins, const CFuint nbTemp, const CFuint nbPress, 
 const CFreal* Ttable, const CFreal* Ptable, const CFreal* opacities, 
 const CFreal* radSource, CFreal T, CFreal p, CFuint ib, CFreal& val1, CFreal& val2)
{
  //Find the lower bound fo the temperature and the pressure ranges
  //we assume that the temperature and pressure always fall in the bounds.
  //If they don't then the value are still interpolated from the nearest
  //two points in the temperature or pressure list
  CFuint it = nbTemp - 2;
  for (CFuint i = 1; i < (nbTemp - 2); i++){
    if(Ttable[i] > T) { it = i - 1; break;}
  }
  
  CFuint ip = nbPress - 2;
  for (CFuint i = 1; i < (nbPress - 2); i++){
    if(Ptable[i] > p) { ip = i - 1; break;}
  }
  
  //Linear interpolation for the pressure
  
  const CFuint iPiBiT           = it + ib*nbTemp + ip*nbBins*nbTemp;
  const CFuint iPplus1iBiT      = it + ib*nbTemp + (ip + 1)*nbBins*nbTemp;
  const CFuint iPiBiTplus1      = (it + 1) + ib*nbTemp + ip*nbBins*nbTemp;
  const CFuint iPplus1iBiTplus1 = (it + 1) + ib*nbTemp + (ip + 1)*nbBins*nbTemp;
  
  // Linear interpolation for the pressure
  // Interpolation of the opacities
  const CFreal bt1op = (opacities[iPplus1iBiT] - opacities[iPiBiT])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + opacities[iPiBiT];
  
  const CFreal bt2op = (opacities[iPplus1iBiTplus1] - opacities[iPiBiTplus1])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + opacities[iPiBiTplus1];
  
  // Interpolation of the source
  const CFreal bt1so = (radSource[iPplus1iBiT] - radSource[iPiBiT])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + radSource[iPiBiT];
  
  const CFreal bt2so = (radSource[iPplus1iBiTplus1] - radSource[iPiBiTplus1])*
		    (p - Ptable[ip])/(Ptable[ip + 1] - Ptable[ip]) + radSource[iPiBiTplus1];    
  
  // Logarithmic interpolation for the temperature
  // Protect against log(0) and x/0 by switching to linear interpolation if either
  // bt1 or bt2 == 0.  (Note we can't allow log of negative numbers either)
  // Interpolation of the opacities   
  if(bt1op <= 0 || bt2op <= 0){
    val1 = (bt2op - bt1op)*(T - Ttable[it])/(Ttable[it + 1] - Ttable[it]) + bt1op;
//    cout <<"\nOption1 \n";
//    cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//    cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  else {
    val1 = std::exp((T - Ttable[it])/(Ttable[it + 1] - Ttable[it])*std::log(bt2op/bt1op))*bt1op;
//     cout <<"\nOption2 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  // Interpolation of the source
  if(bt1so <= 0 || bt2so <= 0){
    val2 = (bt2so - bt1so)*(T - Ttable[it])/(Ttable[it + 1] - Ttable[it]) + bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val1 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  else {
    val2 = std::exp((T - Ttable[it])/(Ttable[it + 1] - Ttable[it])*std::log(bt2so/bt1so))*bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<Ttable[it + 1]<<"\tTi = "<<Ttable[it] <<"\n";
//     cout <<"val2 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  
  //cf_assert(ib == 0);
}
      
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
    
    tableInterpolate(nbBins, nbTemp, nbPress, Ttable, Ptable,
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

__global__ void computeQKernel(const CFuint ib, 
			       const CFuint nbCells,
			       const CFuint* cellFaces,
			       const CFint* isOutward,
			       const CFint* advanceOrder,
			       const CFreal* states,
			       const CFreal* volumes,
			       const CFreal* fieldSource,
			       const CFreal* fieldAbsor,
			       const CFreal* fieldAbSrcV,
			       const CFreal* fieldAbV)
  
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

/*void RadiativeTransferFVDOMCUDA::computeQ(const CFuint ib, const CFuint d)
{      
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::computeQ() in (bin, dir) = ("
	<< ib << ", " << d << ") => start\n");
  
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  SafePtr<TopologicalRegionSet> cells = geoData.trs;
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  cf_assert(m_advanceOrder[d].size() == nbCells);
  
  // precompute the dot products for all faces and directions (a part from the sign)
  RealVector dotProdInFace;
  computeDotProdInFace(d, dotProdInFace);
  SafePtr<ConnectivityTable<CFuint> > cellFaces = MeshDataStack::getActive()->getConnectivity("cellFaces");
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  
  for (CFuint m = 0; m < nbCells; m++) {
    CFreal inDirDotnANeg = 0.;
    CFreal Ic            = 0.;
    
    // allocate the cell entity
    const CFuint iCell = std::abs(m_advanceOrder[d][m]);
    
    // new algorithm (more parallelizable): opacities are computed cell by cell
    // for a given bin
    if (!m_oldAlgo) {getFieldOpacities(ib, iCell);} 
    
    const CFuint nbFaces = cellFaces->nbCols(iCell);
    
    if(m_useExponentialMethod){
      inDirDotnANeg = 0.;
      CFreal dirDotnANeg = 0;
      CFreal Lc      = 0;
      CFreal halfExp = 0;
      
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	const CFuint faceID = (*cellFaces)(iCell, iFace);
	const CFreal factor = ((CFuint)(isOutward[faceID]) != iCell) ? -1. : 1.;
	const CFreal dirDotNA = dotProdInFace[faceID]*factor;
	
	if(dirDotNA < 0.) {
	  dirDotnANeg += dirDotNA;
	  
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
      Lc          = volumes[iCell]/(- dirDotnANeg); 
      halfExp     = std::exp(-0.5*Lc*m_fieldAbsor[iCell]);
      m_In[iCell] = (inDirDotnANeg/dirDotnANeg)*halfExp*halfExp + (1. - halfExp*halfExp)*m_fieldSource[iCell];
      Ic          = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*m_fieldSource[iCell];
    }
    else{
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
    }
    
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
    
    m_divq[iCell] += inDirDotnA*m_weight[d];
    m_II[iCell]   += Ic*m_weight[d];
  }  
  
  CFLog(VERBOSE, "RadiativeTransferFVDOMCUDA::computeQ() in (bin, dir) = ("
	<< ib << ", " << d << ") => end\n");
}*/
      
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
    
    for(CFuint ib = startBin; ib < endBin; ++ib) {
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
      
      // compute the radiative heat flux
      computeQKernel<<<blocksPerGrid,nThreads>>> 
	(ib, nbCells,
	 cellFaces->getPtr()->ptrDev(),
	 isOutward.getLocalArray()->ptrDev(),
	 m_advanceOrder[d].ptrDev(),
	 states.getGlobalArray()->ptrDev(),
	 volumes.getLocalArray()->ptrDev(),
	 m_fieldSource.ptrDev(),
	 m_fieldAbsor.ptrDev(),
	 m_fieldAbSrcV.ptrDev(),
	 m_fieldAbV.ptrDev());
    }    
    
    CFLog(INFO, ")\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

