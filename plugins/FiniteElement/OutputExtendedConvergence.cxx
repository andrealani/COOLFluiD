#include "Common/PE.hh"
#include "Common/SafePtr.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/OutputExtendedConvergence.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<OutputExtendedConvergence,
                      DataProcessingData,
                      FiniteElementModule>
outputExtendedConvergenceProvider("OutputExtendedConvergence");

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file");
   options.addConfigOption< std::string >("OutputFile","Name of Output File");

}

//////////////////////////////////////////////////////////////////////////////

OutputExtendedConvergence::OutputExtendedConvergence(const std::string& name) :
  DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  _maxValues(),
  _minValues()
{

  m_outFile = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

   addConfigOptionsTo(this);
  _nameOutputFile = "ExtraOutput.plt";
   setParameter("OutputFile",&_nameOutputFile);

  _saveRate = 1;
   setParameter("SaveRate",&_saveRate);

}

//////////////////////////////////////////////////////////////////////////////

OutputExtendedConvergence::~OutputExtendedConvergence()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
OutputExtendedConvergence::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::setup()
{
  CFAUTOTRACE;

  prepareOutputFile();
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::execute()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % _saveRate)) {
    computeExtraValues();
  }

}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::computeExtraValues()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _minValues.resize(nbEqs);
  _maxValues.resize(nbEqs);
  _stateValue.resize(nbEqs);
  _averageState.resize(nbEqs);

  _minValues = std::numeric_limits<CFreal>::max();
  _maxValues = -std::numeric_limits<CFreal>::min();

  //in case the state is not found
  _stateValue = std::numeric_limits<CFreal>::max();
  _averageState = 0.;

  for(CFuint iState=0; iState < nbStates; ++iState)
  {
    for(CFuint iEq=0; iEq < nbEqs; ++iEq)
    {
      if((*(states[iState]))[iEq] > _maxValues[iEq]) _maxValues[iEq] = (*(states[iState]))[iEq];
      if((*(states[iState]))[iEq] < _minValues[iEq]) _minValues[iEq] = (*(states[iState]))[iEq];
    }

    //if state is
    if((std::fabs(states[iState]->getCoordinates()[XX] - 25.) < MathTools::MathConsts::CFrealEps()) &&
       (std::fabs(states[iState]->getCoordinates()[YY] - 0.) < MathTools::MathConsts::CFrealEps()))
    {
      for(CFuint iEq=0; iEq < nbEqs; ++iEq)
      {
        _stateValue[iEq] = (*(states[iState]))[iEq];
      }
    }

    for(CFuint iEq=0; iEq < nbEqs; ++iEq)
    {
      _averageState[iEq] += (*(states[iState]))[iEq];
    }
  }
  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    _averageState[iEq] /= nbStates;
  }


  // Output to file the coefficients
  updateOutputFile();
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::prepareOutputFile()
{
  boost::filesystem::path fpath = _nameOutputFile;

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  convergenceFile <<"TITLE  =  Extended Convergence Info"  << "\n";
  convergenceFile << "VARIABLES = Iter PhysTime";

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " MaxValue[" << iEq <<"]";
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " MinValue[" << iEq <<"]";
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " StateValue[" << iEq <<"]";
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " AverageState[" << iEq <<"]";
  }

  convergenceFile << "\n";
  convergenceFile.close();

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergence::updateOutputFile()
{
  boost::filesystem::path fpath (_nameOutputFile);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  CFreal time = subSysStatus->getCurrentTimeDim();
  if(SubSystemStatusStack::getActive()->doingSubIterations()
     && !SubSystemStatusStack::getActive()->isSubIterationLastStep() )
  {
    time += SubSystemStatusStack::getActive()->getDTDim();
  }

  convergenceFile << subSysStatus->getNbIter()
    << " "
    << time;

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " " << _maxValues[iEq];
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " " << _minValues[iEq];
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " " << _stateValue[iEq];
  }

  for(CFuint iEq=0; iEq < nbEqs; ++iEq)
  {
    convergenceFile << " " << _averageState[iEq];
  }

  convergenceFile << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

