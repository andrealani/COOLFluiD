// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <sstream>

#include "Common/PE.hh"
#include "Common/ProcessInfo.hh"
#include "Common/OSystem.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/CFEnvVars.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/StopConditionController.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string > ("ConvergenceFile","Name of Convergence File where to write the convergence history.");
   options.addConfigOption< std::string > ("SpaceResidualFile","Name of SpaceResidualFile where to write the spatial residual history.");
   options.addConfigOption< bool > ("OutputSpaceResidual","Indicate if the Space Residual must be outputted seperately for implicit time marching");
   options.addConfigOption< CFuint >   ("ConvRate","Rate to save convergence data to convergence file.");
   options.addConfigOption< CFint >    ("ScreenOutputPrecision","Screen Output Precision.");
   options.addConfigOption< CFuint >   ("ShowRate","Rate to show convergence message to stdout.");
   options.addConfigOption< bool >     ("ConvergenceFileOnlyP0","Indictate if only the processor 0 should write to file.");
   options.addConfigOption< std::string > ("StopCondition","The stop condition to control the iteration procedure.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::build_dynamic_functions()
{
  DynamicFunctionCaller<ConvergenceMethod>::add_dynamic_function("takeStep", &ConvergenceMethod::takeStep);
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::registActionListeners()
{
  Method::registActionListeners();

  //   const std::string ssname = SubSystemStatusStack::getCurrentName();  
  //   event_handler->addListener(event_handler->key(ssname, "CF_ON_CONVERGENCE_TAKESTEP"),
  //                              this,&ConvergenceMethod::solve);
}
    
//////////////////////////////////////////////////////////////////////////////

ConvergenceMethod::ConvergenceMethod(const std::string& name)
  : Method(name),
    m_lss(),
    m_stopwatch()
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  m_stopConditionStr = "AbsoluteNormAndMaxIter";
  setParameter("StopCondition",&m_stopConditionStr);

  m_nameConvergenceFile = "convergence.plt";
  setParameter("ConvergenceFile",&m_nameConvergenceFile);
  
  m_outputSpaceResidual = false;
  setParameter("OutputSpaceResidual",&m_outputSpaceResidual);
  
  m_nameSpaceResidualFile = "spaceResidual.plt";
  setParameter("SpaceResidualFile",&m_nameSpaceResidualFile);

  m_convRate = 1;
  setParameter("ConvRate",&m_convRate);

  m_showRate = 1;
  setParameter("ShowRate",&m_showRate);

  m_precision = 8;
  setParameter("ScreenOutputPrecision",&m_precision);

  m_onlyP0 = true;
  setParameter("ConvergenceFileOnlyP0",&m_onlyP0);
}

//////////////////////////////////////////////////////////////////////////////

ConvergenceMethod::~ConvergenceMethod()
{
  CFLog(VERBOSE, "ConvergenceMethod::~ConvergenceMethod()\n");
}

//////////////////////////////////////////////////////////////////////////////

std::string ConvergenceMethod::getFileName(const std::string name) const
{
  return name + "." + getMethodData()->getNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Method::configure(args);

  // builds the stop condition
  Common::SafePtr<StopCondition::PROVIDER> stopCondProv =
    FACTORY_GET_PROVIDER(getFactoryRegistry(), StopCondition, m_stopConditionStr);
  SelfRegistPtr<StopCondition> stopCondition = stopCondProv->create(stopCondProv->getName());
  
  configureNested ( stopCondition.getPtr(), args );
  
  // Pass the stop condition to the stop condition controller
  m_stopCondControler.reset(StopConditionController::Create(stopCondition));
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::takeStep()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  if (m_stopwatch.isNotRunning()) { m_stopwatch.start(); }

  takeStepImpl();
  if ( hasToUpdateConv() ) updateConvergenceFile();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::CFL> ConvergenceMethod::getCFL()
{
  return getConvergenceMethodData()->getCFL();
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::setMethodImpl()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ConvergenceMethod::setMethodImpl() => start\n");

  if ( hasToUpdateConv() ) {
    prepareConvergenceFile();
  }

  m_stopwatch.reset();

  SubSystemStatusStack::getActive()->adimensionalizeTimeData();

  CFLog(VERBOSE, "ConvergenceMethod::setMethodImpl() => end\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::unsetMethodImpl()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::prepareConvergenceFile()
{
  path fpath = getFileName(m_nameConvergenceFile);
  fpath = Environment::DirPaths::getInstance().getResultsDir() /
    Framework::PathAppender::getInstance().appendParallel( fpath );
  
  SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  ofstream& convergenceFile = (*fhandle)->open(fpath);

  convergenceFile << "TITLE  =  Convergence of SubSystem" << "\n";
  convergenceFile << "VARIABLES = Iter";

  SafePtr<Framework::ComputeNorm> norm_computer = getConvergenceMethodData()->getNormComputer();
  for(CFuint iVar = 0; iVar < norm_computer->getComputedNormVarIDs().size(); iVar++)
  {

    convergenceFile << " Res[" << norm_computer->getComputedNormVarIDs()[iVar] << "]";
  }
  convergenceFile << " CFL PhysTime DT WallTime MemUsage" << "\n";

  (*fhandle)->close();
  delete fhandle; 
 
  if (m_outputSpaceResidual) {
    path fpath = getFileName(m_nameSpaceResidualFile);
    fpath = Environment::DirPaths::getInstance().getResultsDir() /
            Framework::PathAppender::getInstance().appendParallel( fpath );

    SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
    ofstream& spaceResidualFile = (*fhandle)->open(fpath);

    spaceResidualFile << "TITLE  =  Space Residual of SubSystem" << "\n";
    spaceResidualFile << "VARIABLES = Iter";
    
    SafePtr<Framework::ComputeNorm> norm_computer = getConvergenceMethodData()->getNormComputer();
    for(CFuint iVar = 0; iVar < norm_computer->getComputedNormVarIDs().size(); iVar++)
    {

      spaceResidualFile << " Res[" << norm_computer->getComputedNormVarIDs()[iVar] << "]";
    }
    spaceResidualFile << " CFL PhysTime DT WallTime MemUsage" << "\n";

    (*fhandle)->close();
    delete fhandle; 
 }
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::updateConvergenceFile()
{
  CFAUTOTRACE;

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  if(!(subSysStatus->getNbIter() % m_convRate))
  {
    path fpath = getFileName(m_nameConvergenceFile);
    fpath = Environment::DirPaths::getInstance().getResultsDir() /
            Framework::PathAppender::getInstance().appendParallel( fpath );

    SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
    ofstream& convergenceFile = (*fhandle)->open(fpath, ios_base::app);

    convergenceFile << subSysStatus->getNbIter()       << " "
                    << subSysStatus->getAllResiduals() << " "
                    << getConvergenceMethodData()->getCFL()->getCFLValue()    << " "
                    << subSysStatus->getCurrentTimeDim() << " "
                    << subSysStatus->getDTDim()          << " "
                    << m_stopwatch.read()             << " "
                    << OSystem::getInstance().getProcessInfo()->memoryUsageBytes()  << "\n";

    (*fhandle)->close();
    delete fhandle;
  }
  
  if (outputSpaceResidual()) {
    if(!(subSysStatus->getNbIter() % m_convRate))
    {
      path fpath = getFileName(m_nameSpaceResidualFile);
      fpath = Environment::DirPaths::getInstance().getResultsDir() /
              Framework::PathAppender::getInstance().appendParallel( fpath );

      SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
      ofstream& spaceResidualFile = (*fhandle)->open(fpath, ios_base::app);

      spaceResidualFile << subSysStatus->getNbIter()       << " "
                      << getConvergenceMethodData()->getSpaceResidual() << " "
                      << getConvergenceMethodData()->getCFL()->getCFLValue()    << " "
                      << subSysStatus->getCurrentTimeDim() << " "
                      << subSysStatus->getDTDim()          << " "
                      << m_stopwatch.read()             << " "
                      << OSystem::getInstance().getProcessInfo()->memoryUsageBytes()  << "\n";

      (*fhandle)->close();
      delete fhandle; 
   }
  }
}


//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::syncGlobalDataComputeResidual(const bool computeResidual)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  const bool isParallel = Common::PE::GetPE().IsParallel();
  Common::Stopwatch<Common::WallTime> syncTimer;
  
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(getNamespace());
  DataHandle<State*, GLOBAL> statedata = 
    MeshDataStack::getInstance().getEntryByNamespace(nsp)->getStateDataSocketSink().getDataHandle();
  DataHandle<Node*, GLOBAL> nodedata = 
    MeshDataStack::getInstance().getEntryByNamespace(nsp)->getNodeDataSocketSink().getDataHandle();
  
  if (CFEnv::getInstance().getVars()->SyncAlgo != "Old") {
    if (isParallel) {
      statedata.synchronize();
      nodedata.synchronize();
    }
    if (computeResidual) {
      getConvergenceMethodData()->updateResidual();
    }
  }
  else {
    // after each update the states have to be syncronized
    if (isParallel) {
      syncTimer.start();
      statedata.beginSync ();
    }
    
    if (computeResidual) {
      getConvergenceMethodData()->updateResidual();
    }
    
    if (isParallel) {
      statedata.endSync();
      syncTimer.stop();
    }
  }

  popNamespace();
  
  CFLog(VERBOSE, "ConvergenceMethod::syncGlobalDataComputeResidual()\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::syncAllAndComputeResidual(const bool computeResidual)
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  const bool isParallel = Common::PE::GetPE().IsParallel();
  Common::Stopwatch<Common::WallTime> syncTimer;
  
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(getNamespace());
  DataHandle<State*, GLOBAL> statedata = 
    MeshDataStack::getInstance().getEntryByNamespace(nsp)->getStateDataSocketSink().getDataHandle();
  DataHandle<Node*, GLOBAL> nodedata = 
    MeshDataStack::getInstance().getEntryByNamespace(nsp)->getNodeDataSocketSink().getDataHandle();
  
  // after each update the states have to be syncronized
  if (isParallel)
  {
    syncTimer.start();
    statedata.beginSync ();
    nodedata.beginSync ();
  }

  if (computeResidual)
  {
    getConvergenceMethodData()->updateResidual();
  }

  if (isParallel)
  {
    statedata.endSync();
    nodedata.endSync();
    syncTimer.stop();
  }

  popNamespace();
  
  CFLog(VERBOSE, "ConvergenceMethod::syncAllAndComputeResidual()\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvergenceMethod::writeOnScreen()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());

  pushNamespace();

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  const CFuint iter = subSysStatus->getNbIter();
  if( m_showRate != 0 && !(iter % m_showRate))
  {
    ostringstream out;
    out.precision(m_precision);
    out.setf( std::ios::left );

    //    out.setf(ios::scientific,ios::floatfield);

    // print subsystem status name
    if(subSysStatus->getName() != "SubSystemStatus") {
      out << "[ ";
      out.width(25);
      out <<  subSysStatus->getName() << "] ";
    }
    
    // Print number of iterations
    out << "Iter: ";
    out.width(6);
    out <<  subSysStatus->getNbIter() << " ";

    // Print Residual
    
    out.precision(m_precision);
    //out.setf( std::ios::fixed, std:: ios::floatfield );
    out.width(3 + m_precision);
    out << " Res: [" << subSysStatus->getAllResiduals() << "] ";

    out.unsetf(std::ios::fixed );
    // Print Current CFL
    out << " CFL: ";
    out.width( m_precision );
    out <<  getConvergenceMethodData()->getCFL()->getCFLValue() ;
    // Print Current Time (if Unsteady)
    if (subSysStatus->getDT() > 0.) {
      out << " PhysTime:";
      out.width( m_precision );
      out << subSysStatus->getCurrentTimeDim()*PhysicalModelStack::getActive()->getImplementor()->getTimeFactor();
      // Print Current Time
      out << " DT:";
      out.width( m_precision );
      out << subSysStatus->getDTDim()*PhysicalModelStack::getActive()->getImplementor()->getTimeFactor();
    }

    // Print CPU Time
    out << " CPUTime: ";
    out.width( m_precision );
    out << subSysStatus->readWatch() << " ";

    // Print Memory Usage
    out.precision( m_precision );
    out.setf( std::ios::fixed, std:: ios::floatfield );
    //out.setf( std::ios::left );
    out << "Mem: " << Common::OSystem::getInstance().getProcessInfo()->memoryUsage() << "\n";
    
    CFLog(INFO, out.str() << "\n");
    
    popNamespace();
  }
}

//////////////////////////////////////////////////////////////////////////////

bool ConvergenceMethod::hasToUpdateConv()
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  // AL: should we use "Default" here to avoid I/O troubles
  return !m_onlyP0 || ( m_onlyP0 && (Common::PE::GetPE().GetRank(nsp) == 0) );
}
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
