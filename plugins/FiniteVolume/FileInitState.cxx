

#include "Common/CFLog.hh"
#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MeshData.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FileInitState.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<FileInitState, CellCenterFVMData, FiniteVolumeModule> FileInitStateProvider("FileInitState");

//////////////////////////////////////////////////////////////////////////////

void FileInitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("FileName","Name of file with states.");
  options.addConfigOption< bool >("AdimensionalValues","Flag to input adimensional values.");
}

//////////////////////////////////////////////////////////////////////////////

FileInitState::FileInitState(const std::string& name) :
  CellCenterFVMCom(name),
  _varSet(CFNULL)
{
   addConfigOptionsTo(this);
   setParameter("FileName",&_fileName);

  _inputAdimensionalValues = false;
   setParameter("AdimensionalValues",&_inputAdimensionalValues);

}

//////////////////////////////////////////////////////////////////////////////

FileInitState::~FileInitState()
{
}

//////////////////////////////////////////////////////////////////////////////

void FileInitState::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void FileInitState::setup()
{
  CellCenterFVMCom::setup();

  _varSet = getMethodData().getUpdateVar();
}

//////////////////////////////////////////////////////////////////////////////

void FileInitState::execute()
{
  boost::filesystem::path fname =
    Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(_fileName);

  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& file = (*fhandle)->open(fname);

  CFuint nbStates;
  CFuint nbEqs;
  file >> nbStates >> nbEqs;

  CFout << "Number of states = " << nbStates << "\n";
  CFout << "Length of state vector = " << nbEqs << "\n";

  SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  if (nbCells != nbStates)
  {
    CFout << "File has " << nbStates << " states while there are " << nbCells << " cells \n";
        throw Common::FilesystemException (FromHere(), std::string( "Bad number of states inside the file: ") + fname.string());
  }
  else
  {
    CFout << "OK\n";
  }

  GeometricEntityPool<TrsGeoWithNodesBuilder> geoBuilder;
  geoBuilder.setup();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder.getDataGE();
  geoData.trs = cells;

  State dimState;
  for (CFuint iCell=0; iCell<nbCells; ++iCell)
  {
    geoData.idx = iCell;
    GeometricEntity& cell = *geoBuilder.buildGE();

    State* pstate = cell.getState(0);
    cf_assert(pstate);

    if ( pstate->size() != nbEqs)
    {
      throw Common::FilesystemException (FromHere(), std::string( "Bad length of vector state inside the file: ") +  fname.string());
    }

    if(_inputAdimensionalValues){
      for ( CFuint k=0; k<nbEqs; ++k)
      {
        file >> (*pstate)[k];
      }
    }
    else
    {
      for ( CFuint k=0; k<nbEqs; ++k)
      {
        file >> dimState[k];
      }
      _varSet->setAdimensionalValues(dimState, *pstate);
    }

    geoBuilder.releaseGE();
  }

  (*fhandle)->close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////
//
// void FileInitState::executeOnTrs()
// {
//   cf_assert(0);
//
//   std::string workingDir = Environment::DirPaths::getInstance().getWorkingDir();
//
//   ifstream file( (workingDir +_fileName).c_str());
//
//   if ( !file )
//      throw Common::FilesystemException (FromHere(), workingDir +_fileName );
//
//   CFout << workingDir << _fileName << "\n";
//
//   CFLogDebugMin( "" << "::write() => opened file: " << _fileName << "\n");
//
//   CFuint nst, neq;
//   file >> nst >> neq;
//
//   CFout << "number of states = " << nst << "\n";
//   CFout << "length of state vector = " << neq << "\n";
//
//   SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
//   CFLogDebugMax( "FileInitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");
//
//   vector<State*>* trsStates = trs->getStatesInTrs();
//
//   if ( trsStates->size() != nst )
//   {
//     CFout << "file n " << nst << " cool n " << trsStates->size() << "\n";
//      throw Common::FilesystemException (FromHere(), std::string( "bad number of states inside the file: ") + workingDir +_fileName );
//   }
//   else
//      CFout << "OK\n";
//
//   vector<State>      tab;
//   tab.resize( nst);
//
//   for ( CFuint ii=0; ii<nst; ++ii)
//     file >> tab[ii][0] >> tab[ii][1] >> tab[ii][2] >> tab[ii][3] >> tab[ii][4];
//
//   CFout << "finished tmp tab\n";
//
//   std::vector<State*>::iterator itd;
//   CFuint idx = 0;
//   for (itd = trsStates->begin(); itd != trsStates->end(); ++itd, ++idx)
//   {
//     //CFout << (*itd)->getLocalID() << " \n ";
//
//      //if (idx != (*itd)->getLocalID())
//      //{
//              //CFout << idx << "  " << (*itd)->getLocalID() << "  ERROR\n";
//       //}
//
//      if ( (*(*itd)).size() != neq)
//       throw Common::FilesystemException (FromHere(), std::string( "bad length of vector state inside the file: ") + workingDir +_fileName );
//
//      (*(*itd)) = tab[ (*itd)->getLocalID() ];
//
// //           for ( CFuint k=0; k<neq; ++k)
// //   {
// //     file >> (*(*itd))[k];
// //   }
//
//      //CFout << " [ " << (*(*itd))[0] << ", "  << (*(*itd))[1] << ", "  << (*(*itd))[2] << ", "  << (*(*itd))[3] << ", "  << (*(*itd))[4] << "]\n" ;
//     //_vFunction.evaluate((*itd)->getCoordinates(),*(*itd));
//   }
//
//   CFout << "read OK\n";
//
//
//
//   file.close();
//   CFLogDebugMin( "" << "::write() => closed file: " << _fileName << "\n");
//   CFLogDebugMin( "" << "::write() => end" << "\n");
// }
//
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
