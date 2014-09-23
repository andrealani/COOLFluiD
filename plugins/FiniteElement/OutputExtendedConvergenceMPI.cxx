#include "Common/PE.hh"
#include "Common/MPI/MPIStructDef.hh"
#include "Common/SafePtr.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/OutputExtendedConvergenceMPI.hh"
#include "Framework/PhysicalModel.hh"

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

MethodCommandProvider<OutputExtendedConvergenceMPI,
                      DataProcessingData,
                      FiniteElementModule>
OutputExtendedConvergenceMPIProvider("OutputExtendedConvergenceMPI");

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("SaveRate","Rate for saving the output file");
   options.addConfigOption< std::string >("OutputFile","Name of Output File");

}

//////////////////////////////////////////////////////////////////////////////

OutputExtendedConvergenceMPI::OutputExtendedConvergenceMPI(const std::string& name) :
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

OutputExtendedConvergenceMPI::~OutputExtendedConvergenceMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
OutputExtendedConvergenceMPI::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::setup()
{
  CFAUTOTRACE;

  if (PE::GetPE().GetRank() == 0) {
    prepareOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::execute()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % _saveRate)) {
    computeExtraValues();
  }

}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::computeExtraValues()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  const CFuint nbStates = states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  _minValues.resize(nbEqs);
  _maxValues.resize(nbEqs);
  _averageState.resize(nbEqs);
  _minValuesGlobal.resize(nbEqs);
  _maxValuesGlobal.resize(nbEqs);
  _averageStateGlobal.resize(nbEqs);

  _minValues = std::numeric_limits<CFreal>::max();
  _maxValues = -std::numeric_limits<CFreal>::min();
  _minValuesGlobal = std::numeric_limits<CFreal>::max();
  _maxValuesGlobal = -std::numeric_limits<CFreal>::min();

  _averageState = 0.;
  _averageStateGlobal = 0.;

  for(CFuint iState=0; iState < nbStates; ++iState)
  {
    for(CFuint iEq=0; iEq < nbEqs; ++iEq)
    {
      if((*(states[iState]))[iEq] > _maxValues[iEq]) _maxValues[iEq] = (*(states[iState]))[iEq];
      if((*(states[iState]))[iEq] < _minValues[iEq]) _minValues[iEq] = (*(states[iState]))[iEq];
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

  MPI_Reduce(&_minValues[0], &_minValuesGlobal[0], nbEqs,
     MPI_DOUBLE, MPI_MIN, 0, PE::GetPE().GetCommunicator());

  MPI_Reduce(&_maxValues[0], &_maxValuesGlobal[0], nbEqs,
     MPI_DOUBLE, MPI_MAX, 0, PE::GetPE().GetCommunicator());

  MPI_Reduce(&_averageState[0], &_averageStateGlobal[0], nbEqs,
     MPI_DOUBLE, MPI_SUM, 0, PE::GetPE().GetCommunicator());


CFout <<"Out\n";
  if (PE::GetPE().GetRank() == 0) {
CFout <<"Averaging\n";
    for(CFuint iEq=0; iEq < nbEqs; ++iEq)
    {
      _averageStateGlobal[iEq] /= PE::GetPE().GetProcessorCount();
    }
  }

  // Output to file the coefficients
  if (PE::GetPE().GetRank() == 0) {
CFout <<"Writing to file\n";
    updateOutputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::prepareOutputFile()
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
    convergenceFile << " AverageState[" << iEq <<"]";
  }

  convergenceFile << "\n";
  convergenceFile.close();

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void OutputExtendedConvergenceMPI::updateOutputFile()
{
  boost::filesystem::path fpath (_nameOutputFile);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  convergenceFile << subSysStatus->getNbIter()
    << " "
    << subSysStatus->getCurrentTimeDim();

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

