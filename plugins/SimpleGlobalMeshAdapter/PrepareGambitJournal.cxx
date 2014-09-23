// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "PrepareGambitJournal.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/OSystem.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PrepareGambitJournal, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> PrepareGambitJournalProvider("PrepareGambitJournal");

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("JournalFile","Name of the journal file.");
  options.addConfigOption< std::string >("NewJournalFile","Name of the modified journal file.");

  options.addConfigOption< std::vector<std::string> >("GambitParameters","Definition of the Parameters in the journal file.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");

}

//////////////////////////////////////////////////////////////////////////////

PrepareGambitJournal::PrepareGambitJournal(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{
   addConfigOptionsTo(this);

   _journalFile = "";
   setParameter("OriginalJournalFile",&_journalFile);

   _newJournalFile = "";
   setParameter("NewJournalFile",&_newJournalFile);

   _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

   _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);

   _gambitParams = std::vector<std::string>();
   setParameter("GambitParameters",&_gambitParams);

}

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  SimpleMeshAdapterCom::configure(args);

  // parsing the functions that the user inputed
  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  cf_assert(_gambitParams.size() == _functions.size());
}

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::setup()
{

  SimpleMeshAdapterCom::setup();

  _variables.resize(1);
  _result.resize(_gambitParams.size());
}

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::execute()
{
  CFAUTOTRACE;

  const std::string otherNameSpace = getMethodData().getOtherNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(otherNameSpace);
  Common::SafePtr<SubSystemStatus> subsystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

  _variables[0] = subsystemStatus->getCurrentTimeDim();

  //Evaluate the value of each of the gambit parameters
  _vFunction.evaluate(_variables,_result);

  //Read old journal file
  //check if parameter is in the line
  //if it is, replace it by his value
  //Write the new journal file
  if (PE::GetPE().IsParallel()) {

    PE::GetPE().setBarrier();

    if (PE::GetPE().GetRank () == 0) {
      readWriteFile();
    }

    PE::GetPE().setBarrier();
  }
  else{
    readWriteFile();
  }

}

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::readWriteFile()
{
  CFAUTOTRACE;

  boost::filesystem::path oldFile =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(_journalFile);

  Common::SelfRegistPtr<Environment::FileHandlerInput> oldJournalHandle =
     Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& oldJournal = oldJournalHandle->open(oldFile);

  boost::filesystem::path newFile =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(_newJournalFile);

  SelfRegistPtr<Environment::FileHandlerOutput> newJournalHandle =
     Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = newJournalHandle->open(newFile);

  CFuint lineNb = 0;
  std::string line;
  vector<std::string> words;

  // read line
///@todo change this to read the full file
CFuint nbLines =10;
  for (CFuint i=0; i< nbLines; ++i)
  {
    getWordsFromLine(oldJournal,line,lineNb,words);

    const CFuint nbWords = words.size();
    for(CFuint iWord=0; iWord< nbWords; ++iWord)
    {
      const CFuint nbParameters = _gambitParams.size();
      for(CFuint iPar=0; iPar< nbParameters; ++iPar)
      {
        if(words[iWord] == _gambitParams[iPar]) words[iWord] = StringOps::to_str(_result[iPar]);
      }
    }
    //output the modified line into a new file
    for(CFuint iWord=0; iWord< nbWords; ++iWord)
    {
      fout << words[iWord] << " ";
    }
    fout << "\n";
  }

  oldJournalHandle->close();
  newJournalHandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void PrepareGambitJournal::getWordsFromLine(ifstream& fin,
                                            std::string& line,
                                            CFuint&  lineNb,
                                            vector<std::string>& words)
{
  getline(fin,line);
  ++lineNb;
  words = Common::StringOps::getWords(line);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
