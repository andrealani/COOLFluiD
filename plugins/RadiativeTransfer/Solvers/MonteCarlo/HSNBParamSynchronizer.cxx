#include "HSNBParamSynchronizer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {


HSNBParamSynchronizer::HSNBParamSynchronizer()
{
    m_nbSyncSteps=0;
    m_avgRecvBufByteSize=0;
    m_avgSendBufByteSize=0;
    m_maxRecvBufByteSize=0;
    m_maxSendBufByteSize=0;
    m_avgSyncDuration=0;
    m_avgSyncDurationGlobal=0;
    m_avgTotalBufByteSize=0;
    m_avgTotalBufByteSizeGlobal=0;
    m_maxTotalBufByteSize=0;
    m_maxTotalBufByteSizeGlobal=0;

}



void HSNBParamSynchronizer::setup(COOLFluiD::CFuint mRank, COOLFluiD::CFuint nbProcs, CFuint nbThickDiatomics, CFuint nbNonThickDiatomics, CFuint nbContinua, CFuint nbAtoms, CFuint nbSpecies)
{
m_rank=mRank;
m_nbProcs=nbProcs;
m_nbSpecies=nbSpecies;

m_nbThickDiatomics=nbThickDiatomics;
m_nbNonThickDiatomics=nbNonThickDiatomics;
m_nbAtoms=nbAtoms;
m_nbContinua=nbContinua;
m_nbDiatomics=nbThickDiatomics+nbNonThickDiatomics;

const std::string nsp = Framework::MeshDataStack::getActive()->getPrimaryNamespace();
m_comm = Common::PE::GetPE().GetCommunicator(nsp);

photonTraces.resize(m_nbProcs);

m_sendCountToProcess.resize(m_nbProcs);
m_recvCountFromProcess.resize(m_nbProcs);

debugRecvProcs.resize(m_nbProcs);
debugRecvTraceCount.resize(m_nbProcs);
debugRecvRealCount.resize(m_nbProcs);
m_nbTracesCommitted=0;

m_recvAllBuffer.clear();
m_sendAllBuffer.clear();
}

void HSNBParamSynchronizer::reset()
{

    CFLog(DEBUG_MAX, "HSNBParamSynchronizer::reset" << "\n");
    //std::cout << "P" << m_rank << ": HSNBParamSynchronizer::reset" << "\n";

    cf_assert(photonTraces.size()==m_nbProcs);

    for (int i=0; i<photonTraces.size(); i++) {
        photonTraces[i].clear();
    }

    photonTraces.clear();
    photonTraces.resize(m_nbProcs);

    m_nbTracesCommitted=0;
    m_totalRecvBufCount=0;
    m_totalSendBufCount=0;

    debugRecvProcs.clear();
    debugRecvProcs.resize(m_nbProcs);

    debugRecvTraceCount.clear();
    debugRecvTraceCount.resize(m_nbProcs);

    debugRecvRealCount.clear();
    debugRecvRealCount.resize(m_nbProcs);

    m_recvAllBuffer.clear();
    m_sendAllBuffer.clear();

    for(CFuint i=0; i<m_sendCountToProcess.size(); ++i){
      m_sendCountToProcess[i]=0;
      m_recvCountFromProcess[i]=0;
    }
}

void HSNBParamSynchronizer::bufferCommitParticle(HSNBDataContainer &newParticle)
{
    CFreal addCount;

    m_curHeader=newParticle.headerData;   


    photonTraces[m_curHeader->targetProcessId].push_back(boost::move(newParticle.trace));


    addCount=newParticle.trace->getBufferRealCount();
    m_sendCountToProcess[m_curHeader->targetProcessId]+=addCount;
    m_totalSendBufCount+=addCount;
    m_nbTracesCommitted++;

    //std::cout << "HSNBParamSynchronizer::bufferCommitParticle => newParticle.trace->m_nbCellsCrossed=" << newParticle.trace->m_nbCellsCrossed << ", m_curHeader->nbCrossedCells=" << m_curHeader->nbCrossedCells << std::endl;

    cf_assert(newParticle.trace->m_nbCellsCrossed==m_curHeader->nbCrossedCells);

}

void HSNBParamSynchronizer::synchronize(std::vector<Photon> &photonStack)
{
    m_nbSyncSteps++;

    CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => Start synchro step: \n");
    Common::Stopwatch<WallTime> stp;
    stp.start();

    setupSendContainer();

    CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => Total sendbuffer size is "<< m_totalSendBufCount*sizeof(CFreal)<< " bytes. / " << m_totalSendBufCount << " elements. \n");
    CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => PhotonStack.size()="<< photonStack.size() <<  "\n");

    setupReceiveContainer(photonStack);

    //Compute performance metrics
    if (m_sendBufferByteSize>m_maxSendBufByteSize) {
        m_maxSendBufByteSize=m_sendBufferByteSize;
    }
    if (m_recvBufferByteSize>m_maxRecvBufByteSize) {
        m_maxRecvBufByteSize=m_recvBufferByteSize;
    }
    if ((m_recvBufferByteSize+m_sendBufferByteSize)>m_maxTotalBufByteSize) {
        m_maxTotalBufByteSize=m_recvBufferByteSize+m_sendBufferByteSize;
    }
    m_avgSendBufByteSize+=(int(m_sendBufferByteSize)-int(m_avgSendBufByteSize))/int(m_nbSyncSteps);
    m_avgRecvBufByteSize+=(int(m_recvBufferByteSize)-int(m_avgRecvBufByteSize))/int(m_nbSyncSteps);
    //CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => m_avgRecvBufByteSize= "<< m_avgRecvBufByteSize<< ", m_recvBufferByteSize=" << m_recvBufferByteSize << ", m_nbSyncSteps= " << m_nbSyncSteps<< "\n");

    m_avgTotalBufByteSize+=(int((m_sendBufferByteSize+m_recvBufferByteSize))-int(m_avgTotalBufByteSize))/m_nbSyncSteps;

    CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => Total recvbuffer size is "<< m_totalRecvBufCount*sizeof(CFreal)<< " bytes. / " << m_totalSendBufCount << " elements. \n");

    CFuint nbPhotonsSend=0;
    std::vector<int> displacements(m_nbProcs);
    for(CFuint i=0; i< m_nbProcs ; ++i ){
      displacements[i] = nbPhotonsSend;
      nbPhotonsSend += m_sendCountToProcess[i];
    }

    //get the displacements for photons to receive
    std::vector<int> recvDisps(m_nbProcs);

    CFuint nbPhotonsRecv=0;
    for(CFuint i=0; i< m_nbProcs ; ++i ){
      recvDisps[i]   = nbPhotonsRecv;
      nbPhotonsRecv += m_recvCountFromProcess[i];
    }


    //Scatter continous buffer
    MPI_Alltoallv(&m_sendAllBuffer[0], &m_sendCountToProcess[0], &displacements[0],
          Common::MPIStructDef::getMPIType(&m_sendAllBuffer[0]), &m_recvAllBuffer[0], &m_recvCountFromProcess[0],
          &recvDisps[0],Common::MPIStructDef::getMPIType(&m_recvAllBuffer[0]), m_comm);


    //clear the sendbuffers
    m_sendAllBuffer.clear();

    for(CFuint i=0; i<m_sendCountToProcess.size(); ++i){
      m_sendCountToProcess[i]=0;
      m_recvCountFromProcess[i]=0;
    }

    parseReceiveBuffer(photonStack);

    stp.stop();
    m_avgSyncDuration+=(stp.read()-m_avgSyncDuration)/m_nbSyncSteps;
    CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => Sync finished in "<<stp.read()<< "\n");


//    for (int i=0; i<debugRecvTraceCount.size(); i++) {
//       CFLog(VERBOSE, "HSNBParamSynchronizer::synchronize => debugRecvTraceCount[" << i << "]=" << debugRecvTraceCount[i]<< ", RealCount: "<< debugRecvRealCount[i]<< " \n");
//   }

    MPI_Reduce(&m_avgTotalBufByteSize, &m_avgTotalBufByteSizeGlobal, 1,
          Common::MPIStructDef::getMPIType(&m_avgTotalBufByteSize), MPI_SUM, 0, m_comm);
    MPI_Reduce(&m_avgSyncDuration, &m_avgSyncDurationGlobal, 1,
          Common::MPIStructDef::getMPIType(&m_avgSyncDuration), MPI_SUM, 0, m_comm);
    MPI_Reduce(&m_maxTotalBufByteSize, &m_maxTotalBufByteSizeGlobal, 1,
          Common::MPIStructDef::getMPIType(&m_maxTotalBufByteSize), MPI_MAX, 0, m_comm);


    if (m_rank==0) {
        m_avgTotalBufByteSizeGlobal/=m_nbProcs;
        m_avgSyncDurationGlobal/=m_nbProcs;
    }

//    printPerformanceMetrics();


}

const vector<HSNBDataContainer> &HSNBParamSynchronizer::getTraceStack() const
{
    return m_traceStack;
}

void HSNBParamSynchronizer::clearTraceStackMemory()
{
    m_traceStack.clear();
}

HSNBDataContainer HSNBParamSynchronizer::getTraceStackAt(CFuint i) const
{
    return m_traceStack[i];
}

void HSNBParamSynchronizer::debugSendToZero()
{
    CFuint sendCount=0;
    CFuint addCount=0;

    std::cout << "P" <<m_rank <<" send to Zero " << photonTraces[0].size() << " / Realcount: " << sendCount << " // Addcount: " << addCount << ". m_sendCountToProcess[0]=" << m_sendCountToProcess[0] <<  " \n";

    for (int i=0; i<photonTraces[0].size(); i++) {
        sendCount+=photonTraces[0][i]->getBufferRealCount();
        addCount+=traceRealCount(photonTraces[0][i]->m_nbCellsCrossed);
    }

    for (int i=0; i<photonTraces[0].size(); i++) {
        photonTraces[0][i]->print(true);
    }


 }


int HSNBParamSynchronizer::totalSendBufferByteSize() const
{
    return m_totalSendBufCount*sizeof(CFreal);
}

void HSNBParamSynchronizer::printSendTraces(bool printAll)
{

    if (printAll) {
        std::cout << "HSNBParamSynchronizer::print " << m_rank << " \n";
        std::cout << "HSNBParamSynchronizer::print => " << m_rank << ". Number of traces committed=" << m_nbTracesCommitted << " \n";
        std::cout << "HSNBParamSynchronizer::print => "<< m_rank << ". BufferByteSize=" << totalSendBufferByteSize() << " \n";

        std::cout << "HSNBParamSynchronizer::print => Traces commited to process " << m_rank << "\n";
        for (int i=0; i<photonTraces.size(); i++) {
            std::cout << "\n HSNBParamSynchronizer::print => " << m_rank << " Traces to be sent to process "<< i << "\n ";
            for (int j=0; j<photonTraces[i].size(); j++) {
                std::cout << "------------------------------------ \n";
                std::cout << "\n HSNBParamSynchronizer::print " << m_rank << ", Trace[" << i  << "]["<<j << "]: \n ";
                photonTraces[i][j]->print();
            }
        }
    }
    else {
        CFLog(INFO, "HSNBParamSynchronizer::print \n");
        CFLog(INFO, "HSNBParamSynchronizer::print => Number of traces committed=" << m_nbTracesCommitted << " \n");
        CFLog(INFO, "HSNBParamSynchronizer::print => BufferByteSize=" << totalSendBufferByteSize() << " \n");



        //  std::cout << "BUFFER BYTE SIZE " << totalBufferByteSize() << std::endl;
        CFLog(INFO, "HSNBParamSynchronizer::print => Traces commited to process " << m_rank << "\n");
        for (int i=0; i<photonTraces.size(); i++) {
            CFLog(INFO, "\n HSNBParamSynchronizer::print => Traces to be sent to process "<< i << "\n ");
            for (int j=0; j<photonTraces[i].size(); j++) {
                CFLog(INFO, "------------------------------------ \n");
                CFLog(INFO, "\n HSNBParamSynchronizer::print, Trace[" << i  << "]["<<j << "]: \n ");
                photonTraces[i][j]->print();
            }
        }
    }
}

void HSNBParamSynchronizer::printRecvTraces(bool printAll)
{
    if (printAll) {
        std::cout << "HSNBParamSynchronizer::printRecvTraces => " << m_rank << ". Print Trace stack for process " << m_rank <<". Total Tracecount=" << m_traceStack.size() << "\n";
    }
    else {
        CFLog(INFO, "HSNBParamSynchronizer::printRecvTraces => Print Trace stack for process " << m_rank <<". Total Tracecount=" << m_traceStack.size() << "\n");
    }

    for (int i=0; i<m_traceStack.size(); i++) {
        if (printAll) {
            std::cout << "\n HSNBParamSynchronizer::printRecvTraces => " << m_rank << ". Print trace " << i << "\n";
        }
        else {
            CFLog(INFO, "\n HSNBParamSynchronizer::printRecvTraces => Print trace " << i << "\n");
        }

        m_traceStack[i].print(printAll);
    }


}



void HSNBParamSynchronizer::appendThinParamSetToBuffer(HSNBNonThickParameterSet &paramSet)
{
    //CFLog(VERBOSE, "HSNBParamSynchronizer::appendThinParamSetToBuffer => ADD NONTHICK PARAMSET \n");

    //Scalar for nonthick systems
    m_sendAllBuffer.push_back(paramSet.kappa);
}

void HSNBParamSynchronizer::appendThickParamSetToBuffer(HSNBThickParameterSet &paramSet)
{
    //CFLog(VERBOSE, "HSNBParamSynchronizer::appendThickParamSetToBuffer => ADD THICK PARAMSET \n");

    appendVector(m_sendAllBuffer,paramSet.kappa);
    appendVector(m_sendAllBuffer,paramSet.betaD);
    appendVector(m_sendAllBuffer,paramSet.betaL);
}

void HSNBParamSynchronizer::appendAtomicParamSetToBuffer(HSNBAtomicParameterSet &paramSet)
{
    //CFLog(VERBOSE, "HSNBParamSynchronizer::appendAtomicParamSetToBuffer => ADD ATOMIC PARAMSET \n");
    m_sendAllBuffer.push_back(paramSet.optThick);
}

void HSNBParamSynchronizer::appendTraceToBuffer(TracePtr newTrace)
{

   //Push message header first
   m_sendAllBuffer.push_back(newTrace->m_nbCellsCrossed);
   m_sendAllBuffer.push_back(newTrace->m_kappa0);
   m_sendAllBuffer.push_back(newTrace->m_mechLocalID);
   m_sendAllBuffer.push_back(newTrace->m_distance0);


   cf_assert(m_nbDiatomics==(newTrace->m_thinDiatomics.size()+newTrace->m_thickDiatomics.size()));
   cf_assert(m_nbContinua==newTrace->m_continua.size());
   cf_assert(m_nbAtoms==newTrace->m_atoms.size());
   cf_assert(m_nbNonThickDiatomics==newTrace->m_thinDiatomics.size());
   cf_assert(m_nbThickDiatomics==newTrace->m_thickDiatomics.size());


   for (int i=0; i<m_nbThickDiatomics; i++) {
       appendThickParamSetToBuffer(newTrace->m_thickDiatomics[i]);
   }

   for (int i=0; i<m_nbNonThickDiatomics; i++) {
       appendThinParamSetToBuffer(newTrace->m_thinDiatomics[i]);
   }

   for (int i=0; i<newTrace->m_nbContinua; i++) {
       appendThinParamSetToBuffer(newTrace->m_continua[i]);
   }

   for (int i=0; i<newTrace->m_nbAtoms; i++) {
       appendAtomicParamSetToBuffer(newTrace->m_atoms[i]);
   }

}


void HSNBParamSynchronizer::setupReceiveContainer(std::vector<Photon> &photonStack)
{
    //Set the size of the recvBuffer to fit all Data to be received for this process
    m_totalRecvBufCount=0;
    CFreal addCount;



//    for (int i=0; i<photonStack.size(); i++) {
//        debugRecvProcs[photonStack[i].userData.prevProcessId]++;
//    }

//    for (int i=0; i<m_nbProcs; i++) {
//        CFLog(VERBOSE, "HSNBParamSynchronizer::setupReceiveContainer => debugRecvProcs[" << i <<"]=" <<debugRecvProcs[i] << " \n");
//    }

    for (int i=0; i<photonStack.size(); i++) {
        addCount=traceRealCount(photonStack[i].userData.nbCrossedCells);
        //CFLog(VERBOSE, "HSNBParamSynchronizer::setupReceiveContainer => nbCrossedCells=" << photonStack[i].userData.nbCrossedCells << " \n");
        m_recvCountFromProcess[photonStack[i].userData.prevProcessId]+=addCount;
        m_totalRecvBufCount+=addCount;

        //ALL YOUR PHOTONS ARE BELONG TO US
        //photonStack[i].userData.prevProcessId=m_rank;
    }

    CFLog(VERBOSE, "HSNBParamSynchronizer::setupReceiveContainer => m_totalRecvBufCount=" << m_totalRecvBufCount << " \n");
    //Allocate continous memory to receive in

//    for (int i=0; i<m_recvCountFromProcess.size(); i++) {
//        CFLog(VERBOSE, "HSNBParamSynchronizer::setupReceiveContainer => m_recvCountFromProcess[" << i << "]=" << m_recvCountFromProcess[i] << " \n");
//    }

    m_recvAllBuffer.resize(m_totalRecvBufCount);
    m_recvBufferByteSize=m_totalRecvBufCount*sizeof(CFreal);
}

void HSNBParamSynchronizer::setupSendContainer()
{
    m_sendAllBuffer.reserve(m_totalSendBufCount);

    //std::cout << "PhotonTraces.size=" << photonTraces.size() <<"=!=" << m_nbProcs << "\n";

    for (int pi=0; pi<photonTraces.size(); pi++) {

        for (int ti=0; ti<photonTraces[pi].size(); ti++) {
            this->appendTraceToBuffer(photonTraces[pi][ti]);
        }

        //std::cout << "PROC " << m_rank << " FINISHED FIRST LOOP NO DELETE, PHOTONTRACES[" << pi<< "].size()=" <<  photonTraces[pi].size() << std::endl;

        //DONT FORGET TO FREE MEMORY ALWAYS (MEMORY IS LIFE)
        photonTraces[pi].clear();
    }


    //std::cout << "PROC " << m_rank << " FINISHED " << sendAllBuffer.size() << " " << m_totalSendBufCount<< std::endl;


    cf_assert(m_sendAllBuffer.size()==m_totalSendBufCount);
    m_sendBufferByteSize=totalSendBufferByteSize();
}

void HSNBParamSynchronizer::parseReceiveBuffer(std::vector<Photon> &photonStack)
{
    CFuint nbCells;
    CFuint emitMechID;
    CFreal kappa0;
    CFreal distance0;

    CFuint photonStackID=0;


    CFuint curPos=0;
    TracePtr curTrace;

    clearReceiveParameters();

    cf_assert(m_traceStack.size()==0);

    //while (recvAllBuffer.size()>TraceStaticRealCount){
    //We need to make sense of the continous datastream in the Buffer:
    while (m_recvAllBuffer.empty()==false){
        curPos=0;

        //std::cout << "P" << m_rank << " ENTERS LOOP, REMAINING BUFSIZE " << recvAllBuffer.size() << " \n";

        nbCells=CFuint(m_recvAllBuffer[0]);
        kappa0=m_recvAllBuffer[1];
        emitMechID=m_recvAllBuffer[2];
        distance0=m_recvAllBuffer[3];
        curPos=TraceStaticRealCount;

        //Setup new DataContainer
        m_traceStack.push_back(HSNBDataContainer());
        m_traceStack.back().trace=TracePtr(new HSNBPhotonTrace);
        m_traceStack.back().headerData=&photonStack[photonStackID].userData;


        debugRecvRealCount[photonStack[photonStackID].userData.prevProcessId]+=traceRealCount(nbCells);
        debugRecvTraceCount[photonStack[photonStackID].userData.prevProcessId]++;
        photonStack[photonStackID].userData.prevProcessId=m_rank;

        photonStackID++;

        //copy ptr to safe ptr temporarily
        curTrace= m_traceStack.back().trace;

        curTrace->m_nbDiatomics=m_nbDiatomics;
        curTrace->m_nbContinua=m_nbContinua;
        curTrace->m_nbAtoms=m_nbAtoms;
        curTrace->m_nbSpecies=m_nbSpecies;

        //Parse static variables
        curTrace->m_nbCellsCrossed=nbCells;
        curTrace->m_kappa0=kappa0;
        curTrace->m_mechLocalID=emitMechID;
        curTrace->m_distance0=distance0;

        //Check writing order
        for (int mi=0; mi<m_nbThickDiatomics; mi++) {
            curTrace->m_thickDiatomics.push_back(HSNBThickParameterSet());
            m_tempThickData=&(curTrace->m_thickDiatomics.back());
            m_tempThickData->nbCells=nbCells;

            m_tempThickData->kappa=std::vector<CFreal>(m_recvAllBuffer.begin()+curPos,m_recvAllBuffer.begin()+curPos+nbCells);
            curPos+=nbCells;
            m_tempThickData->betaD=std::vector<CFreal>(m_recvAllBuffer.begin()+curPos,m_recvAllBuffer.begin()+curPos+nbCells);
            curPos+=nbCells;
            m_tempThickData->betaL=std::vector<CFreal>(m_recvAllBuffer.begin()+curPos,m_recvAllBuffer.begin()+curPos+nbCells);
            curPos+=nbCells;
        }

        for (int mi=0; mi<m_nbNonThickDiatomics; mi++) {
            curTrace->m_thinDiatomics.push_back(HSNBNonThickParameterSet());
            m_tempNonThickData=&(curTrace->m_thinDiatomics.back());
            m_tempNonThickData->nbCells=nbCells;

            m_tempNonThickData->kappa=m_recvAllBuffer[curPos];
            curPos++;
        }

        for (int mi=0; mi<m_nbContinua; mi++) {
            curTrace->m_continua.push_back(HSNBNonThickParameterSet());
            m_tempNonThickData=&(curTrace->m_continua.back());
            m_tempNonThickData->nbCells=nbCells;

            m_tempNonThickData->kappa=m_recvAllBuffer[curPos];
            //m_tempNonThickData->kappa=std::vector<CFreal>(recvAllBuffer.begin()+curPos,recvAllBuffer.begin()+curPos+nbCells);
            curPos++;

        }

        for (int mi=0; mi<m_nbAtoms; mi++) {
            curTrace->m_atoms.push_back(HSNBAtomicParameterSet());
            m_tempAtomicData=&(curTrace->m_atoms.back());
            m_tempAtomicData->nbCells=nbCells;

            m_tempAtomicData->optThick=m_recvAllBuffer[curPos];
            curPos++;
        }


        cf_assert(curPos==traceRealCount(nbCells));

        //CHECK
        m_recvAllBuffer.erase(m_recvAllBuffer.begin(), m_recvAllBuffer.begin()+curPos);
    }

    cf_assert(photonStack.size()==m_traceStack.size());
}

void HSNBParamSynchronizer::printAll()
{

    for (int i=0; i<photonTraces.size(); i++) {
        std::cout << " TRACE COUNT " << photonTraces[i].size() << std::endl;
        for (int j=0; j<photonTraces[i].size(); j++) {
            std::cout << photonTraces[i][j]->name() << std::endl;
        }
    }
}

void HSNBParamSynchronizer::printPerformanceMetrics() const
{

    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => nbSyncSteps=" << m_nbSyncSteps << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => curSendBufByteSize=" << m_sendBufferByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgSendBufByteSize=" << m_avgSendBufByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => curRecvBufByteSize=" << m_recvBufferByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgRecvBufByteSize=" << m_avgRecvBufByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => maxRecvBufByteSize=" << m_maxRecvBufByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => maxSendBufByteSize=" << m_maxSendBufByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgTotalBufByteSize=" << m_avgTotalBufByteSize << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgTotalBufByteSizeGlobal=" << m_avgTotalBufByteSizeGlobal << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics =>  m_maxTotalBufByteSizeGlobal=" <<  m_maxTotalBufByteSizeGlobal << " bytes / " << m_maxTotalBufByteSizeGlobal/1e6<< " MB. \n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgSyncDuration=" << m_avgSyncDuration << "\n");
    CFLog(INFO, "HSNBParamSynchronizer::printPerformanceMetrics => avgSyncDurationGlobal=" << m_avgSyncDurationGlobal << "\n");
}


CFuint HSNBParamSynchronizer::traceRealCount(CFuint nbCells) const
{
    CFuint count=0;

    //Static variables in trace class (e.g. nbCrossedCells)
    count+=TraceStaticRealCount;

    //Add trace parameters for all mechanisms
    count+=nbThickParameters*nbCells*m_nbThickDiatomics;

    count+=nbAtomicParameters*m_nbAtoms;
    count+=nbNonThickParameters*m_nbNonThickDiatomics;
    count+=nbNonThickParameters*m_nbContinua;

    return count;
}


void HSNBParamSynchronizer::appendVector(std::vector<CFreal>& v1, std::vector<CFreal>& v2)
{
    v1.insert(v1.end(),v2.begin(),v2.end() );
}

void HSNBParamSynchronizer::appendVector2(std::vector<CFreal> &v1, std::vector<CFreal> &v2)
{
    for (int i=0; i<v2.size(); i++) {
        v1.push_back(2.0);
    }
}

void HSNBParamSynchronizer::clearReceiveParameters()
{
    m_traceStack.clear();
}




    }
}
