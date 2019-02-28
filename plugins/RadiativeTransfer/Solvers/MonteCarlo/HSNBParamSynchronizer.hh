#ifndef HSNBPARAMSYNCHRONIZER_H
#define HSNBPARAMSYNCHRONIZER_H

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "Config/ConfigObject.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DofDataHandleIterator.hh"

#include "Common/COOLFluiD.hh"
#include <vector>
#include "RadiativeTransfer/RadiationLibrary/Models/HSNB/HSNBRadiator.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/PhotonData.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonData.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBLocalParameterSet.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/HSNBPhotonTrace.hh"
#include "LagrangianSolver/LagrangianSolver.hh"

#include "Common/MPI/MPIStructDef.hh"
#include "Common/CFPrintContainer.hh"
#include "Common/MPI/MPIError.hh"
#include "Common/PE.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "Common/Stopwatch.hh"
#include "Common/SafePtr.hh"

const COOLFluiD::CFuint nbAtomicParameters = 1;
const COOLFluiD::CFuint nbNonThickParameters = 1;
const COOLFluiD::CFuint nbThickParameters = 3;
const COOLFluiD::CFuint nbCO2Parameters = 1;
const COOLFluiD::CFuint nbSynchHeaderParameters=2;

namespace COOLFluiD {

    namespace RadiativeTransfer {

    typedef LagrangianSolver::Particle<HSNBPhotonData> Photon;
      
//////////////////////////////////////////////////////////////////////////////

    //!
    //! \brief Class for the synchronization of HSNBPhotonTrace structures between processes
    //!
    //! For the correct computation of HSNB absorption a specific set of information is needed for
    //! all the cells crossed by the photon along its trace. This set is stored in the HSNBPhotonTrace
    //! struct and has to be passed between processes if the computational domain is crossed.
    //! The HSNBParamSynchronizer is responsible of handling the serialization of all Traces committed to
    //! it using bufferCommitParticle(). The serial data is pushed into a single buffer of real values
    //! m_sendAllBuffer and distributed among processes using the MPI AllToAllV routine (in exactly the same
    //! fashion as in the LagrangianSolver module). The receive buffer then has to be parsed "manually" / in
    //! a predetermined manner to reconstruct the HSNBPhotonTraces and push them to the traceStack.
    //!
    //! Be really careful when changing something with the datastructure of PhotonTraces because the synchronization
    //! might easily break if you change the ordering of real values as they are sent and received.
    //!
class HSNBParamSynchronizer {
public:
    HSNBParamSynchronizer();
    void setup(CFuint mRank, CFuint nbProcs, CFuint nbThickDiatomics, CFuint nbNonThickDiatomics, CFuint nbContinua, CFuint nbAtoms, CFuint nbSpecies, bool co2Exists=false);
    void reset();

    /// Commit a new particle
    void bufferCommitParticle(HSNBDataContainer& newParticle);

    /// Synchronize the serialized sendbuffer using MPI (the pattern is the same as in the
    /// LagrangianSolver, it might be easier to understand by looking at the sincronize
    /// routine in there since the dataset to be synced is more simple there.
    void synchronize(std::vector<Photon> &photonStack);

    const vector<HSNBDataContainer>& getTraceStack() const;

    void clearTraceStackMemory();

    HSNBDataContainer getTraceStackAt(CFuint i) const;

    void debugSendToZero();
    void debugOrderRecvTraces();


    int totalSendBufferByteSize() const;

    void printSendTraces(bool printAll=false);
    void printRecvTraces(bool printAll=false);

    void printAll();
    void printPerformanceMetrics() const;

    /**
     * @brief Checks whether all processes together reserve buffers with more than "byte" bytes of memory usage.
     * Requires a synchronization operation between all processes.
     * @param byte
     * @return true if the total buffer size exceeds "byte"
     */
    bool totalMemoryUsageExceeds(CFreal byteSize);


private:
    bool m_co2Exists=false;

    //Telemetry data
    CFint m_nbSyncSteps;
    CFuint m_avgSendBufByteSize;
    CFuint m_avgRecvBufByteSize;
    CFuint m_maxRecvBufByteSize;
    CFuint m_maxSendBufByteSize;
    CFuint m_maxTotalBufByteSize;
    CFuint m_maxTotalBufByteSizeGlobal;

    CFuint m_avgTotalBufByteSize;
    CFuint m_avgTotalBufByteSizeGlobal;
    CFreal m_avgSyncDuration;
    CFreal m_avgSyncDurationGlobal;

    /// Total size of the continous send buffer in byte
    CFuint m_sendBufferByteSize;
    CFuint m_recvBufferByteSize;


    CFuint m_rank;
    CFuint m_nbProcs;
    CFuint m_nbSpecies;

    CFuint m_nbNonThickDiatomics;
    CFuint m_nbThickDiatomics;
    CFuint m_nbDiatomics;
    CFuint m_nbContinua;
    CFuint m_nbAtoms;

    MPI_Aint m_startAddress;
    MPI_Aint m_address;
    MPI_Comm m_comm;

    //TODO: Have vector of External States. Associate photons with states such
    //that no redundant states are sent


    Common::SafePtr< HSNBPhotonData > m_curHeader;
    Common::SafePtr< HSNBNonThickParameterSet> m_curNonThickParameters;
    Common::SafePtr< HSNBThickParameterSet>  m_curThickParameters;
    Common::SafePtr< HSNBAtomicParameterSet > m_curAtomicParameters;


    /// Unify storage for all processes // mechanisms
    /// where photonTraces(i,j) is the Trace j to be sent to
    /// process i
    std::vector< std::vector< TracePtr > > photonTraces;

    std::vector<CFuint> debugRecvProcs;
    std::vector<CFuint> debugRecvRealCount;
    std::vector<CFuint> debugRecvTraceCount;
    CFuint m_debugCurProc;
    CFuint m_debugCurTrace;

    CFuint m_nbTracesCommitted;

    /// stack of photon traces received from other processes
    /// The ordering corresponds to the photonStack that has been
    /// synchronized in the LagrangianSolver.
    std::vector<HSNBDataContainer> m_traceStack;

    //(Potentially very large) buffer for serialization of trace data
    std::vector< CFreal > m_sendAllBuffer;

    //Receive buffer for serial trace data. Has to be parsed manually in parseRecvData
    std::vector< CFreal > m_recvAllBuffer;

    std::vector< CFint > m_sendCountToProcess;
    std::vector< CFint > m_recvCountFromProcess;

    //Total number of CFreal values in the concat. send Buffer
    CFuint m_totalSendBufCount;
    CFuint m_totalRecvBufCount;

  Common::SafePtr<HSNBThickParameterSet> m_tempThickData;
  Common::SafePtr<HSNBNonThickParameterSet> m_tempNonThickData;
  Common::SafePtr<HSNBAtomicParameterSet> m_tempAtomicData;

    boost::shared_ptr<HSNBThickParameterSet> m_tempThickPtr;
    boost::shared_ptr<HSNBNonThickParameterSet> m_tempNonThickPtr;

private:
  
  void appendParameterSetToBuffer(Common::SafePtr<HSNBLocalParameterSet> paramSet);
  
    void appendThinParamSetToBuffer(HSNBNonThickParameterSet& paramSet);
    void appendThickParamSetToBuffer(HSNBThickParameterSet& paramSet);
    void appendAtomicParamSetToBuffer(HSNBAtomicParameterSet& paramSet);
    void appendCO2ParamSetToBuffer(HSNBCO2ParameterSet &paramSet);

    void appendTraceToBuffer(TracePtr newTrace);
    void appendVector(std::vector<CFreal> &v1, std::vector<CFreal> &v2);
    void appendVector2(std::vector<CFreal> &v1, std::vector<CFreal> &v2);
    void clearReceiveParameters();

    void setupReceiveContainer(std::vector<Photon> &photonStack);
    void setupSendContainer();
    void parseReceiveBuffer(std::vector<Photon> &photonStack);
  
    CFuint traceRealCount(CFuint nbCells) const;

   };

    }
}


#endif // HSNBPARAMSYNCHRONIZER_H
