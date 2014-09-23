#include "LESDataProcessing/LESDataProcessing.hh"
#include "TimeAveraging.hh"
#include "TurbulenceFunction.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TimeAveraging, 
                      LESProcessingData, 
                      LESDataProcessingModule> 
LESDataProcessingTimeAveragingProvider("TimeAveraging");

// LESProcessingComProvider LESDataProcessingTimeAveragingProvider("TimeAveraging");

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool, Config::DynamicOption<> >("Average","Interactive flag to turn on/off the averaging. (default=true)");
  options.addConfigOption< bool >("Reset","Flag that tells the averaging to reset. (default=false)");
  options.addConfigOption< bool >("FirstTimeCreation","Flag that tells to create source socket instead of sink sockets. (default=false)");
  options.addConfigOption< bool >("Nodal","Flag that tells if the solution will be averaged from nodal values");
}

//////////////////////////////////////////////////////////////////////////////

TimeAveraging::TimeAveraging(const std::string& name) 
: LESProcessingCom(name),
  m_socketMap(6)
{
  addConfigOptionsTo(this);
  // by default the data processing is run once -> processRate=infinity
  m_averaging = true;
  setParameter("Average",&m_averaging);
  
  m_resetFlag = false;
  setParameter("Reset",&m_resetFlag);
  
  m_firstTimeCreation = false;
  setParameter("FirstTimeCreation",&m_firstTimeCreation);
  
  m_nodal = false;
  setParameter("Nodal",&m_nodal);
  
  
}

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::configure ( Config::ConfigArgs& args )
{
  LESProcessingCom::configure(args);
  if (m_nodal) {
    m_sockets.createSocketSink<RealVector>("nstates");
    m_globalSockets.createSocketSink<Framework::State*>("states");
  }
  else {
    m_globalSockets.createSocketSink<Framework::State*>("states");
  }

  makeSocketAvailable("averageSolution");
  makeSocketAvailable("averageVelocityProducts");  
  makeSocketAvailable("turbulenceAveragingSteps");
  makeSocketAvailable("averageTurbulenceFluctuations",SOURCE);

}

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::makeSocketAvailable(const std::string& socketName, SocketType socketType)
{
  switch (socketType) {
    case UNDEFINED:
    {
      // Make source socket if source socket doesn't exist
      if (m_firstTimeCreation) {
        m_sockets.createSocketSource<CFreal>(socketName);
        m_socketMap.insert(socketName,SOURCE);
        CFLog (INFO,"   +++ Created source socket for " << socketName << "\n");
      } 
      // Make sink socket if source socket already exists
      else {
        m_sockets.createSocketSink<CFreal>(socketName);
        m_socketMap.insert(socketName,SINK);
        CFLog (INFO,"   +++ Created sink socket for " << socketName << "\n");
      }
      break;
    }
    case SOURCE:
    {
      m_sockets.createSocketSource<CFreal>(socketName);
      m_socketMap.insert(socketName,SOURCE);
      CFLog (INFO,"   +++ Created source socket for " << socketName << "\n");
      break;
    }
    case SINK:
    {
      m_sockets.createSocketSink<CFreal>(socketName);
      m_socketMap.insert(socketName,SINK);
      CFLog (INFO,"   +++ Created source sink for " << socketName << "\n");
      break;
    }  
  }
}

//////////////////////////////////////////////////////////////////////////////

DataHandle<CFreal> TimeAveraging::getDataHandle(const std::string& socketName) 
{
  if (isSourceSocket(socketName)) {
    return m_sockets.getSocketSource<CFreal>(socketName)->getDataHandle();
  } else {
    return m_sockets.getSocketSink<CFreal>(socketName)->getDataHandle();
  }
}

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::setup()
{  
  CFLog(INFO, " +++ TimeAveraging::setup() \n");
  DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();


  if (m_nodal) {
    Framework::DataHandle<RealVector> nodalStates = m_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
    m_nbStates = nodalStates.size();
  }
  else {
    Framework::DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();
    m_nbStates = states.size();
  }

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  m_primState.resize(nbEqs);
  
  switch (dim) {
    case DIM_1D:
      m_stride = 1;
      break;
    case DIM_2D:
      m_stride = 3;
      break;
    case DIM_3D:
      m_stride = 6;
      break;
  }
  
  DataHandle<CFreal> avgSteps = getDataHandle("turbulenceAveragingSteps");   
  if (isSourceSocket("turbulenceAveragingSteps")) {
    avgSteps.resize(1);
    avgSteps = 0;
    m_resetFlag = true;
  }
  if (isSourceSocket("averageSolution")) {
    DataHandle<CFreal> avgSol = getDataHandle("averageSolution");
    avgSol.resize(m_nbStates*nbEqs);
    avgSol = 0.;
    m_resetFlag = true;
  }
  if (isSourceSocket("averageVelocityProducts")) {
    DataHandle<CFreal> avgVelProd =  getDataHandle("averageVelocityProducts");
    avgVelProd.resize(m_nbStates*m_stride);
    avgVelProd = 0.;
    m_resetFlag = true;
  }
  if (isSourceSocket("averageTurbulenceFluctuations")) {
    DataHandle<CFreal> avgTurbFluct = getDataHandle("averageTurbulenceFluctuations");
    avgTurbFluct.resize(m_nbStates*m_stride);
    avgTurbFluct = 0.;
  }
  
  
  if (!m_averaging) {
    m_resetFlag = true;
  }
    
  if (m_resetFlag) {
    // This will reset the averaging
    m_avgStepsCounter = 0;
  } else {
    m_avgStepsCounter = static_cast<CFreal>(avgSteps[0]);
  }
  
  m_gradients.resize(nbEqs);
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
    m_gradients[iEq].resize(dim);
  }
  
  m_turbulenceFunctions = getMethodData().getAverageTurbulenceFunctions();
  
}

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void TimeAveraging::execute()
{
  CFAUTOTRACE;
  
  if(m_averaging) {
    
    if(m_resetFlag) {
      m_avgStepsCounter = 0;
      m_resetFlag = false;
      CFLog(INFO, "\nTime Averaging (re)started \n\n");
    }
    

    Framework::DataHandle<Framework::State*,Framework::GLOBAL> states = m_globalSockets.getSocketSink<Framework::State*>("states")->getDataHandle();
    Framework::DataHandle<RealVector> nodalStates(CFNULL);
    if (m_nodal) {
      nodalStates = m_sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
    }
    
    DataHandle<CFreal> avgSol =       getDataHandle("averageSolution");
    DataHandle<CFreal> avgVelProd =   getDataHandle("averageVelocityProducts");
    DataHandle<CFreal> avgTurbFluct = getDataHandle("averageTurbulenceFluctuations");
    DataHandle<CFreal> avgSteps =     getDataHandle("turbulenceAveragingSteps");

    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbCells = states.size();       
    const CFuint nbTurbulenceFunctions = m_turbulenceFunctions.size();

    // factors for averaging
    const CFreal avgStepsCounterReal = static_cast<CFreal>(m_avgStepsCounter);
    const CFreal oldAvgSolWght = avgStepsCounterReal/(avgStepsCounterReal+1.0);
    const CFreal curInsSolWght = 1.0/(avgStepsCounterReal+1.0);

    // Compute the averages and store in the sockets
    for(CFuint iCell=0; iCell<m_nbStates; ++iCell) {

      // Calculate DIMENSIONAL primitive state
      if (m_nodal) {
        m_primState = getMethodData().transformToPrimDim(nodalStates[iCell]);
      }
      else {
        m_primState = getMethodData().transformToPrimDim((*states[iCell]));
      }

      CFreal u=0. , v=0. , w=0.;
      switch (dim) {
        case DIM_3D:
          w = m_primState[3];
        case DIM_2D:
          v = m_primState[2];
        case DIM_1D:
          u = m_primState[1];
      }

      for (CFuint iEq = 0; iEq<nbEqs; ++iEq) {
        avgSol(iCell,iEq,nbEqs) = avgSol(iCell,iEq,nbEqs)*oldAvgSolWght + m_primState[iEq]*curInsSolWght;
      }
      // Average
      switch (dim) {
        case DIM_1D:
          avgVelProd(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride)*oldAvgSolWght + u*u*curInsSolWght;
          avgTurbFluct(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,1,nbEqs);
          break;
        case DIM_2D:
          avgVelProd(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride)*oldAvgSolWght + u*u*curInsSolWght;
          avgVelProd(iCell,1,m_stride) = avgVelProd(iCell,1,m_stride)*oldAvgSolWght + v*v*curInsSolWght;
          avgVelProd(iCell,2,m_stride) = avgVelProd(iCell,2,m_stride)*oldAvgSolWght + u*v*curInsSolWght;
          avgTurbFluct(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,1,nbEqs);
          avgTurbFluct(iCell,1,m_stride) = avgVelProd(iCell,1,m_stride) - avgSol(iCell,2,nbEqs)*avgSol(iCell,2,nbEqs);
          avgTurbFluct(iCell,2,m_stride) = avgVelProd(iCell,2,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,2,nbEqs);
          break;
        case DIM_3D:
          avgVelProd(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride)*oldAvgSolWght + u*u*curInsSolWght;
          avgVelProd(iCell,1,m_stride) = avgVelProd(iCell,1,m_stride)*oldAvgSolWght + v*v*curInsSolWght;
          avgVelProd(iCell,2,m_stride) = avgVelProd(iCell,2,m_stride)*oldAvgSolWght + w*w*curInsSolWght;
          avgVelProd(iCell,3,m_stride) = avgVelProd(iCell,3,m_stride)*oldAvgSolWght + u*v*curInsSolWght;
          avgVelProd(iCell,4,m_stride) = avgVelProd(iCell,4,m_stride)*oldAvgSolWght + u*w*curInsSolWght;
          avgVelProd(iCell,5,m_stride) = avgVelProd(iCell,5,m_stride)*oldAvgSolWght + v*w*curInsSolWght;
          avgTurbFluct(iCell,0,m_stride) = avgVelProd(iCell,0,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,1,nbEqs);
          avgTurbFluct(iCell,1,m_stride) = avgVelProd(iCell,1,m_stride) - avgSol(iCell,2,nbEqs)*avgSol(iCell,2,nbEqs);
          avgTurbFluct(iCell,2,m_stride) = avgVelProd(iCell,2,m_stride) - avgSol(iCell,3,nbEqs)*avgSol(iCell,3,nbEqs);
          avgTurbFluct(iCell,3,m_stride) = avgVelProd(iCell,3,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,2,nbEqs);
          avgTurbFluct(iCell,4,m_stride) = avgVelProd(iCell,4,m_stride) - avgSol(iCell,1,nbEqs)*avgSol(iCell,3,nbEqs);
          avgTurbFluct(iCell,5,m_stride) = avgVelProd(iCell,5,m_stride) - avgSol(iCell,2,nbEqs)*avgSol(iCell,3,nbEqs);
          break;
      }
    }

    
    // Calculate gradients for every cell and timeAverage turbulence functions
    if (nbTurbulenceFunctions) { // Only if there are turbulence functions to be averaged
      for(CFuint iCell=0; iCell<nbCells; ++iCell) {
        
        // compute gradients
        getMethodData().getGradientComputer()->compute(m_gradients,iCell);
        
        // compute and average turbulence functions
        for(CFuint iFunc=0; iFunc<nbTurbulenceFunctions; ++iFunc) {
          m_turbulenceFunctions[iFunc]->computeAverage(oldAvgSolWght,curInsSolWght,(*states[iCell]),m_gradients,iCell);
        }
      }      
    }
        
    // update m_avgStepsCounter
    ++m_avgStepsCounter;

    avgSteps[0] = m_avgStepsCounter;
  }
  else {
    if (m_avgStepsCounter != 0) {
      CFLog(INFO, "\nTime Averaging stopped after " << m_avgStepsCounter << " time steps. \n\n");
      m_avgStepsCounter = 0;
    }
    m_resetFlag = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LESDataProcessing

  } // namespace Numerics

} // namespace COOLFluiD
