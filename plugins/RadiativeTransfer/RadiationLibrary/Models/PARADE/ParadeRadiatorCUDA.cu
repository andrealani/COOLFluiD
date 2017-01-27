#include <fstream>

#include "RadiativeTransfer/RadiationLibrary/Models/PARADE/ParadeRadiatorCUDA.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

#include "Common/CFLog.hh"
#include "Common/DebugFunctions.hh"
#include "Common/CFPrintContainer.hh"
#include "Framework/CudaTimer.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ParadeRadiatorCUDA,
			    Radiator,
			    RadiativeTransferModule,
			    1>
paradeRadiatorCUDAProvider("ParadeRadiatorCUDA");

//////////////////////////////////////////////////////////////////////////////

__global__ void computeCellBinsKernel(const CFuint nbPoints, 
				      const CFuint nbCells, 
				      const CFuint nbBinsre,
				      const CFuint testID, 
				      const CFreal dWavIn,
				      const CFreal* vctBins,
				      const CFreal* data,
				      CFreal* alpha_bin, 
				      CFreal* emission_bin,
				      CFreal* B_binCurr)
{
  const CFuint idx = threadIdx.x + blockIdx.x*blockDim.x;
  const CFuint cellID = idx%nbCells;
  const CFuint pointID  = idx/nbPoints;
  const CFuint sizeLoop = nbCells*nbPoints;
  
  // fill the alpha and emission arrays for each local cell in parallel simulations
  // each row gives all the bins for a cell
  if (idx < sizeLoop) {
    ParadeRadiator::DeviceFunc df;
    df.computeCellBins<GPU>(nbPoints, pointID, cellID, nbBinsre, testID, dWavIn, vctBins, 
			    data, alpha_bin, emission_bin, B_binCurr);
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void ParadeRadiatorCUDA::defineConfigOptions(Config::OptionList& options)
{
}
  
//////////////////////////////////////////////////////////////////////////////

ParadeRadiatorCUDA::ParadeRadiatorCUDA(const std::string& name) :
  ParadeRadiator(name)
{
  addConfigOptionsTo(this);
}
  
//////////////////////////////////////////////////////////////////////////////

ParadeRadiatorCUDA::~ParadeRadiatorCUDA()
{
}

//////////////////////////////////////////////////////////////////////////////
      
void ParadeRadiatorCUDA::setup()
{
  CFLog(VERBOSE, "ParadeRadiatorCUDA::setup() => START\n");
  
  ParadeRadiator::setup();
    
  CFLog(VERBOSE, "ParadeRadiatorCUDA::setup() => END\n");
}
  
//////////////////////////////////////////////////////////////////////////////
            
void ParadeRadiatorCUDA::unsetup()
{
  CFLog(VERBOSE, "ParadeRadiatorCUDA::unsetup() => START\n");
  
  ParadeRadiator::unsetup();
  
  CFLog(VERBOSE, "ParadeRadiatorCUDA::unsetup() => END\n");
}
  
/////////////////////////////////////////////////////////////////////////////

void ParadeRadiatorCUDA::computeBinning()
{    
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBinning() => START\n");

  CudaEnv::CudaTimer& timer = CudaEnv::CudaTimer::getInstance();
  
  CFuint totalNbCells = m_pstates->getSize();
  const CFuint nbCols = m_nbPoints*3;
  CFuint nbCells = m_data.size()/nbCols;
  cf_assert(m_saveMemory && nbCells <= totalNbCells);
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
     MPI_Allreduce(&nbCells, &totalNbCells, 1, MPIStructDef::getMPIType(&nbCells), 
		   MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
  cf_assert(totalNbCells == m_pstates->getSize());
  
  //Function to calculate the binning algorithm with absorption & emission coefficients
  const CFuint nbBins = m_nbBins;
  CFreal num_alphatot = 0.0;
  CFreal den_alphatot = 0.0;
  m_alphaav.resize(m_nbPoints, 0.);
  LocalArray<CFreal>::TYPE* num_alphatotVec = new LocalArray<CFreal>::TYPE(m_nbPoints, 0.);
  LocalArray<CFreal>::TYPE* den_alphatotVec = new LocalArray<CFreal>::TYPE(m_nbPoints, 0.);
  
  // offset for the ID from which the cells need to be counted 
  CFuint offsetStateID = 0;
  const CFuint nbCellsPerProc = totalNbCells/m_nbProc;
  for (CFuint rank = 0; rank < m_rank; ++rank) {
    offsetStateID += nbCellsPerProc;
  }
  
  timer.start();
  
  for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint s=0;s<nbCells;++s) {
      const CFuint start = s*nbCols + i*3;
      const CFreal alpha = m_data[start+2];
      const CFreal epsilon = m_data[start+1];
      const CFreal Bs = epsilon/alpha;
      
      // Planck function if we are not in equilibrium
      //
      // CFreal h = 6.626070040e-34;     //SI units Js
      // CFuint c = 3e08; //SI units m/s
      // CFreal k_b = 1.3806485279e-23;  //SI units J/K
      //
      // T needs to be defined from the values extracted
      // const  B = (2*h*pow(c,2)/pow(m_data(j,i*3),5))*(1/(exp(h*c/(m_data(j,i*3)*k_b*T))-1));
      const CFuint stateID = s+offsetStateID;
      const CFreal volume = getCellVolume(m_pstates->getStateLocalID(stateID));
      const CFreal BsVolume = Bs*volume;
      const CFreal num_alpha_vol = alpha*BsVolume;
      const CFreal den_alpha_vol = BsVolume;
      num_alphatot += num_alpha_vol;
      den_alphatot += den_alpha_vol;
    }
    
    (*num_alphatotVec)[i] = num_alphatot;
    (*den_alphatotVec)[i] = den_alphatot;
  } 
  
  CFLog(INFO, "ParadeRadiatorCUDA::computeBinning() => num/den alpha took  " << timer.elapsed() << "s \n"); 
  
  /*computeNumDenAlphaKernel<<<blocks, threads>>>
    (nbPoints, nbCells, offsetStateID, 
    m_data.ptrDev(),
    num_alphatotVec.ptrDev(), 
    den_alphatotVec.ptrDev(), 
    m_radPhysicsHandlerPtr->getDataSockets()->volumes.getLocalArray()->ptrDev());*/
  
  // compute the total numerator and denominator per spectral point across all ranks
  vector<CFreal> num_alphatotGlobal(m_nbPoints, 0.);
  vector<CFreal> den_alphatotGlobal(m_nbPoints, 0.);
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
     MPI_Allreduce(&(*num_alphatotVec)[0], &num_alphatotGlobal[0], m_nbPoints, 
		   MPIStructDef::getMPIType(&(*num_alphatotVec)[0]), 
		   MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
     MPI_Allreduce(&(*den_alphatotVec)[0], &den_alphatotGlobal[0], m_nbPoints, 
		   MPIStructDef::getMPIType(&(*den_alphatotVec)[0]), 
		   MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
  
  cf_assert(m_alphaav.size() == m_nbPoints);
  for (CFuint i = 0; i < m_alphaav.size(); ++i) {
    cf_assert(std::abs(den_alphatotGlobal[i]) > 0.);
    m_alphaav[i] = num_alphatotGlobal[i] / den_alphatotGlobal[i];
    // CFLog(INFO, "m_alphaav[" << i << "] = " <<  m_alphaav[i] << "\n");
  }
  
  // free memory
  deletePtr(num_alphatotVec);
  deletePtr(den_alphatotVec);
  
  CFreal alphamin = m_alphaav[0];
  CFreal alphamax = m_alphaav[0];
  for(CFuint r=0;r<m_nbPoints;++r) {
    if(m_alphaav[r]<alphamin) {
      alphamin = m_alphaav[r];
    }
    if(m_alphaav[r]>alphamax) {
      alphamax = m_alphaav[r];
    }
  }
  
  CFLog(INFO,"ParadeLibrary::computeBinning () => [alphamin, alphamax] = [" 
	<< alphamin <<", " << alphamax << "]\n");
  
  // Logarithmic spacing
  //
  const CFreal alpha_minlog = std::log(alphamin);
  const CFreal alpha_maxlog = std::log(alphamax);
  
  CFLog(INFO,"ParadeLibrary::computeBinning () => [alpha_minlog, alpha_maxlog] = [" 
	<< alpha_minlog << ", " << alpha_maxlog << "]\n");
  
  const CFreal dy = (alpha_maxlog-alpha_minlog) / (nbBins-1);
  CFLog(VERBOSE,"ParadeLibrary::computeBinning () => dy = " << dy <<"\n");
  
  LocalArray<CFreal>::TYPE vctBins(0., nbBins);
  for(CFuint i = 0; i<nbBins; ++i) {
    vctBins[i] = std::exp(alpha_minlog + (dy * i));
    CFLog(VERBOSE,"ParadeLibrary::computeBinning () => vctBins(" << i << ") = " << vctBins[i] <<"\n");
  }
  
  // To search for minimum and maximum alpha
  CFreal alphamin_tot = m_alphaav[0];
  CFreal alphamax_tot = m_alphaav[0];
  for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint s=0;s<nbCells;s++) {
      const CFreal alphaA = m_data[s*nbCols + i*3+2];
      if(alphaA < alphamin_tot) {alphamin_tot = alphaA;}
      if(alphaA > alphamax_tot) {alphamax_tot = alphaA;}
    }
  }
  
  if (m_saveMemory) {
    CFreal alphaMinTotLocal = alphamin_tot;
    CFreal alphaMaxTotLocal = alphamax_tot;
    CFreal alphaMinTotGlobal = 0.;
    CFreal alphaMaxTotGlobal = 0.;
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
       MPI_Allreduce(&alphaMinTotLocal, &alphaMinTotGlobal, 1, 
		     MPIStructDef::getMPIType(&alphaMinTotLocal), 
		     MPI_MIN, PE::GetPE().GetCommunicator(m_namespace)));
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
       MPI_Allreduce(&alphaMaxTotLocal, &alphaMaxTotGlobal, 1, 
		     MPIStructDef::getMPIType(&alphaMaxTotLocal), 
		     MPI_MAX, PE::GetPE().GetCommunicator(m_namespace)));
    
    alphamin_tot = alphaMinTotGlobal;
    alphamax_tot = alphaMaxTotGlobal;
  }
  
  vctBins[0] = alphamin_tot;
  vctBins[nbBins-1] = alphamax_tot;
  
  CFLog(VERBOSE,"ParadeLibrary::computeBinning () => [alphamin_tot, alphamax_tot] = [" 
	<< alphamin_tot << ", "  << alphamax_tot <<"\n");
  
  computeAveragedBins(nbBins,2, vctBins);
  
  CFLog(INFO, "ParadeRadiatorCUDA::computeBinning() took " << timer.elapsed() << "s\n"); 
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBinning() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////
  
void ParadeRadiatorCUDA::computeBanding() 
{  
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBanding() => START\n");

  Stopwatch<WallTime> stp;
  stp.start();
  
  const CFuint totalNbCells = m_pstates->getSize();
  const CFuint nbCols = m_nbPoints*3;
  CFuint nbCells = m_data.size()/nbCols;
  cf_assert((!m_saveMemory && nbCells == totalNbCells) || 
	    ( m_saveMemory && nbCells <= totalNbCells));
  
  if (m_saveMemory) {
    CFuint totalNbCells = 0;
    MPIError::getInstance().check
      ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBinning()",
       MPI_Allreduce(&nbCells, &totalNbCells, 1, MPIStructDef::getMPIType(&nbCells), 
		     MPI_SUM, PE::GetPE().GetCommunicator(m_namespace)));
    cf_assert(totalNbCells == m_pstates->getSize());
  }

  // Banding
  CFreal alphamin_tot = m_data[0];
  CFreal alphamax_tot = m_data[0];
  
  for(CFuint i=0;i<m_nbPoints;++i) {
    for(CFuint s=0;s<nbCells;s++) {
      const CFreal wavelength = m_data[s*nbCols + i*3];
      if(wavelength < alphamin_tot) {
	alphamin_tot = wavelength;
      }
      if(wavelength > alphamax_tot) {
	alphamax_tot = wavelength;
      }
    }
  }
  
  CFreal alphamax_totGlobal = 0.;
  CFreal alphamin_totGlobal = 1e10;
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBanding()",
     MPI_Allreduce(&alphamin_tot, &alphamin_totGlobal, 1, 
		   MPIStructDef::getMPIType(&alphamin_tot), 
		   MPI_MIN, PE::GetPE().GetCommunicator(m_namespace)));
  
  MPIError::getInstance().check
    ("MPI_Allreduce", "ParadeRadiatorCUDA::computeBanding()",
     MPI_Allreduce(&alphamax_tot, &alphamax_totGlobal, 1, 
		   MPIStructDef::getMPIType(&alphamax_tot), 
		   MPI_MAX, PE::GetPE().GetCommunicator(m_namespace)));
  
  const CFreal alpha_minlog = std::log(alphamin_totGlobal);
  const CFreal alpha_maxlog = std::log(alphamax_totGlobal);

  CFLog(INFO,"ParadeLibrary::computeBanding () => [alpha_minlog, alpha_maxlog] = ["
	<< alpha_minlog << ", " << alpha_maxlog << "]\n");
  
  const CFreal dy = (alpha_maxlog-alpha_minlog) / (m_nbBands-1);
  
  CFLog(INFO,"ParadeLibrary::computeBanding() => dy = " << dy <<"\n");
  
  LocalArray<CFreal>::TYPE vctBins(0., m_nbBands);
  
  vctBins = 0.;
  for(int i = 0; i<m_nbBands; ++i) {
    vctBins[i] = std::exp(alpha_minlog + (dy * i));
    CFLog(VERBOSE,"ParadeLibrary::computeBanding() => vctBins(" << i << ") = " << vctBins[i] <<"\n");
  }
  
  // AL: I stop here 
  //To search for minimum alpha and maximum
  computeAveragedBins(m_nbBands, 0, vctBins);
  
  CFLog(INFO, "ParadeRadiatorCUDA::computeBanding() took " << stp.read() << "s\n");
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBanding() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ParadeRadiatorCUDA::computeBinningBanding() 
{ 
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBinningBanding() => START\n");
  Stopwatch<WallTime> stp;
  stp.start();
  CFLog(INFO, "ParadeRadiatorCUDA::computeBinningBanding() took "<< stp.read() << "s\n");  
  CFLog(VERBOSE, "ParadeRadiatorCUDA::computeBinningBanding() => END\n");
}
    
//////////////////////////////////////////////////////////////////////////////


  
//////////////////////////////////////////////////////////////////////////////

void ParadeRadiatorCUDA::computeAveragedBins(const CFuint nbBinsre, 
					     const CFuint testID,
					     LocalArray<CFreal>::TYPE& vctBins)
{
  const CFuint totalNbCells = m_pstates->getSize();
  const CFuint nbCols = m_nbPoints*3;
  CFuint nbCells = m_data.size()/nbCols;
  
  //alpha_avbin is average value for absorptivity for each bin
  SafePtr<SocketBundle> sockets = m_radPhysicsHandlerPtr->getDataSockets();
  DataHandle<CFreal> alpha_avbin = sockets->alpha_avbin; // array w GLOBAL cell size
  DataHandle<CFreal> B_bin = sockets->B_bin;             // array w GLOBAL cell size
  
  if (alpha_avbin.size() != nbBinsre*totalNbCells) {
    CFLog(ERROR, "ParadeRadiatorCUDA::computeAveragedBins() => alpha_avbin.size() != nbBinsre*totalNbCells => " << alpha_avbin.size() << " != " <<  nbBinsre*totalNbCells << "\n");
    cf_assert(alpha_avbin.size() == nbBinsre*totalNbCells);
  }
  if (B_bin.size() != nbBinsre*totalNbCells) {
    CFLog(ERROR, "ParadeRadiatorCUDA::computeAveragedBins() => B_bin.size() != nbBinsre*totalNbCells => " << B_bin.size() << " != " <<  nbBinsre*totalNbCells << "\n");
    cf_assert(B_bin.size() == nbBinsre*totalNbCells);
  }
  
  // from now only local arrays (=global if m_saveMemory==false) are used
  m_alpha_bin.resize(nbBinsre*nbCells);     // array w LOCAL cell size
  m_emission_bin.resize(nbBinsre*nbCells);  // array w LOCAL cell size
  
  cf_assert((nbCells < totalNbCells && m_nbProc > 1) || 
	    (nbCells == totalNbCells && m_nbProc == 1));
  
  LocalArray<CFreal>::TYPE alpha_avbinCurr(0., nbBinsre*nbCells);
  LocalArray<CFreal>::TYPE B_binCurr(0., nbBinsre*nbCells);
  
  const CFuint threads = CudaEnv::CudaDeviceManager::getInstance().getNThreads();
  const CFuint blocks = 
    CudaEnv::CudaDeviceManager::getInstance().getBlocksPerGrid(nbCells*m_nbPoints);
  
  CFLog(INFO, "ParadeRadiatorCUDA::computeAveragedBins() => [blocks, threads] = [" 
	<< blocks << ", " << threads << "]\n");
  
  // copy input data from CPU to GPU
  vctBins.put(); 
  m_data.put();
  
  // call a kernel that computes bins 
  computeCellBinsKernel<<<blocks, threads>>>
    (m_nbPoints, nbCells, nbBinsre, testID, m_dWav, vctBins.ptrDev(), m_data.ptrDev(),
     m_alpha_bin.ptrDev(), m_emission_bin.ptrDev(), B_binCurr.ptrDev());
  
  // copy output data from GPU to CPU
  m_alpha_bin.get();
  m_emission_bin.get();
  B_binCurr.get();
  
  /*for(CFuint k=1;k<nbBinsre;++k) {
    for(CFuint j=0;j<nbCells;++j) {
    CFLog(DEBUG_MAX,"ParadeLibrary::computeproperties () => m_alpha_bin(" << k << "," << j << ") = " << m_alpha_bin[nbBinsre*j+k] <<"\n");
    CFLog(DEBUG_MAX,"ParadeLibrary::computeproperties () => B_binCurr(" << k << "," << j << ") = " << B_binCurr[nbBinsre*j+k] <<"\n");
    }
    }*/
  
  for(CFuint j=0;j<nbCells;++j) {
    for(CFuint k=1;k<nbBinsre;++k) {
      const CFuint idx0 = nbBinsre*j;
      alpha_avbinCurr[idx0] = 0.;
      const CFuint idx = k + idx0;
      // AL: is this fix needed to mask an error or it is supposed to be like this?
      if(B_binCurr[idx] != 0.) {
	alpha_avbinCurr[idx] = m_alpha_bin[idx] / B_binCurr[idx];
	CFLog(DEBUG_MED,"ParadeRadiatorCUDA::computeAverageBins() => alpha_avbinCurr(" << k << "," << j << ") = "<< alpha_avbinCurr[idx] <<"\n");
      }
      else {
	alpha_avbinCurr[idx] = 0.;
	CFLog(DEBUG_MED,"ParadeRadiatorCUDA::computeAverageBins() => alpha_avbinCurr(" << k << ","<< j << ") = "<< alpha_avbinCurr[idx] <<"\n");
      }
    }
  }
  
  if (m_saveMemory) {
    // here we need to gather all the entries for alpha_avbin and B_bin from all processes, so that
    // every process keeps a global storage of them
    CFuint minSizeToSend = 0;
    CFuint maxSizeToSend = 0;
    vector<int> recvCounts(m_nbProc, 0);
    vector<int> displs(m_nbProc, 0);
    computeRecvCountsDispls(totalNbCells, nbBinsre, minSizeToSend, maxSizeToSend, recvCounts, displs);
    const CFuint sendSize = nbCells*nbBinsre;
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    cf_assert(sendSize <= maxSizeToSend);
    cf_assert(sendSize >= minSizeToSend);
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiatorCUDA::computeAverageBins() => alpha_avbin",
       MPI_Allgatherv(&alpha_avbinCurr[0], sendSize, MPIStructDef::getMPIType(&alpha_avbinCurr[0]),
		      &alpha_avbin[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&alpha_avbin[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
    
    MPIError::getInstance().check
      ("MPI_Allgatherv", "ParadeRadiatorCUDA::computeAverageBins() => B_bin",
       MPI_Allgatherv(&B_binCurr[0], sendSize, MPIStructDef::getMPIType(&B_binCurr[0]),
		      &B_bin[0], &recvCounts[0], &displs[0],  MPIStructDef::getMPIType(&B_bin[0]),
		      PE::GetPE().GetCommunicator(m_namespace)));
  }
  
  // To be commented out after verification
  if(PE::GetPE().GetRank(m_namespace) == 0) {
   
    /*for(CFuint j=0;j<totalNbCells;++j) {
	for(CFuint k=1;k<nbBinsre;++k) {
	CFLog(DEBUG_MAX,"alpha (" << k << "," << j << ") = " << alpha_avbin[k + nbBinsre*j] << "\n");
	}
	}*/
    ofstream fout1("alpha.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<nbBinsre;++k) {
	fout1 << "alpha (" << k << "," << j << ") = " << 
	  alpha_avbin[k + nbBinsre*j] << "\n";
      }
    }
    fout1.close();
    
    ofstream fout2("beta.txt");
    for(CFuint j=0;j<totalNbCells;++j) {
      for(CFuint k=1;k<nbBinsre;++k) {
	fout2 << "beta (" << k << "," << j << ") = " << 
	  B_bin[k + nbBinsre*j] << "\n";
      }
    }
    fout2.close();
  }
  
  PE::GetPE().setBarrier(m_namespace);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace RadiativeTransfer

} // namespace COOLFluiD

