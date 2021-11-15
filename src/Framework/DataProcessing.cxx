// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/DataProcessing.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<DataProcessing,
               DataProcessingMethod,
               FrameworkLib,
               1>
dataProcessingProvider("DataProcessing");

//////////////////////////////////////////////////////////////////////////////

void DataProcessing::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Comds","Types of the Data Processing commands.");
   options.addConfigOption< std::vector<std::string> >("Names","Names for the configuration of the Data Processing commands.");
   options.addConfigOption< bool >("SkipFirstIteration","Skip the processing after first iteration.");
   options.addConfigOption< bool >("RunAtSetup","Run at setup time.");
   options.addConfigOption< bool >("RunAtSetupAndAfter","Run at setup time and after use process rate.");  
}

//////////////////////////////////////////////////////////////////////////////

DataProcessing::DataProcessing(const std::string& name)
  : DataProcessingMethod(name),
    m_data(),
    m_dataprocessing(0)
{
  addConfigOptionsTo(this);
  m_data.reset(new DataProcessingData(this));
  cf_assert(m_data.getPtr() != CFNULL);

  m_dataprocessTypeStr = vector<std::string>();
  setParameter("Comds",&m_dataprocessTypeStr);

  m_dataprocessNameStr = vector<std::string>();
  setParameter("Names",&m_dataprocessNameStr); 

  m_skipFirstIteration = false;
  setParameter("SkipFirstIteration",&m_skipFirstIteration);
  
  m_runAtSetup = false;
  setParameter("RunAtSetup",&m_runAtSetup);
  
  m_runAtSetupAndAfter = false;
  setParameter("RunAtSetupAndAfter",&m_runAtSetupAndAfter);
}

//////////////////////////////////////////////////////////////////////////////

DataProcessing::~DataProcessing()
{
  clearComs();
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessing::clearComs()
{
  vector<SelfRegistPtr<DataProcessingCom> >().swap(m_dataprocessing);
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<MethodData> DataProcessing::getMethodData () const
{
  return m_data.getPtr();
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessing::setMethodImpl()
{
  DataProcessingMethod::setMethodImpl();

  setupCommandsAndStrategies();
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessing::unsetMethodImpl()
{
  unsetupCommandsAndStrategies();

  DataProcessingMethod::unsetMethodImpl();
}

//////////////////////////////////////////////////////////////////////////////

void DataProcessing::processDataImpl()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "DataProcessing::processDataImpl() -> m_runAtSetup = " << m_runAtSetup << "\n");
  CFLog(VERBOSE, "DataProcessing::processDataImpl() -> SubSystemStatusStack::getActive()->isSetup() = " << SubSystemStatusStack::getActive()->isSetup() << "\n");
  
  // AL: this is needed at least when calculating the distance to the wall
  if ((m_runAtSetup || m_runAtSetupAndAfter) && SubSystemStatusStack::getActive()->isSetup()) {
    for(CFuint i = 0; i < m_dataprocessing.size(); ++i) {
      cf_assert(m_dataprocessing[i].isNotNull()); 
      CFLog(VERBOSE, "DataProcessing " << m_dataprocessing[i]->getClassName() << "->execute() during setup()\n");
      m_dataprocessing[i]->execute();
    }

    if (m_runAtSetup) return;
  }
 
  if ((!SubSystemStatusStack::getActive()->isSetup()) && (!m_runAtSetup)) {
    const CFuint rest = (SubSystemStatusStack::getActive()->getNbIter()-1) % m_processRate; 
    if ((SubSystemStatusStack::getActive()->getNbIter() == 0 && !m_skipFirstIteration) ||
	(rest == 0 && SubSystemStatusStack::getActive()->getNbIter() > 1 && m_skipFirstIteration) ||
	(rest == 0 && SubSystemStatusStack::getActive()->getNbIter() >= 1 && !m_skipFirstIteration))
      {
	//if((SubSystemStatusStack::getActive()->getNbIter() == 0) || ((SubSystemStatusStack::getActive()->getNbIter()-1) % m_processRate == 0))
	//if(!(SubSystemStatusStack::getActive()->getNbIter() % m_processRate == 0) || (SubSystemStatusStack::getActive()->getNbIter() == 0))
	
	CFLog(VERBOSE, "DataProcessing::processDataImpl() => SpaceMethod->extrapolateStatesToNodes()\n");
	// set the nodal states at first
	getMethodData()->getCollaborator<SpaceMethod>()->extrapolateStatesToNodes();
	
	for(CFuint i = 0; i < m_dataprocessing.size(); ++i) {
	  cf_assert(m_dataprocessing[i].isNotNull()); 
	  CFLog(VERBOSE, "DataProcessing " << m_dataprocessing[i]->getClassName() << "->execute()\n");
	  m_dataprocessing[i]->execute();
	}
      }
  } 
}
    
//////////////////////////////////////////////////////////////////////////////

void DataProcessing::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingMethod::configure(args);
  
  m_data->setFactoryRegistry(getFactoryRegistry());
  configureNested ( m_data.getPtr(), args );
  
  clearComs();
  cf_assert(m_dataprocessTypeStr.size() == m_dataprocessNameStr.size());

  m_dataprocessing.resize(m_dataprocessTypeStr.size());

  for(CFuint i = 0; i < m_dataprocessing.size(); ++i) {

    configureCommand<DataProcessingCom,
                     DataProcessingData,
                     DataProcessingComProvider>
      (args, m_dataprocessing[i], m_dataprocessTypeStr[i],m_dataprocessNameStr[i], m_data);

    cf_assert(m_dataprocessing[i].isNotNull());

    CFLog(VERBOSE, "DataProcess type : " << m_dataprocessing[i]->getClassName() << "\n");
    CFLog(VERBOSE, "DataProcess name : " << m_dataprocessing[i]->getName() << "\n");

  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
