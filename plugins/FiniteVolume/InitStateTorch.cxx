#include "FiniteVolume/FiniteVolume.hh"


#include "InitStateTorch.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Common/FilesystemException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateTorch, CellCenterFVMData, FiniteVolumeModule>
initStateTorchProvider("InitStateTorch");

//////////////////////////////////////////////////////////////////////////////

void InitStateTorch::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Datafile","input file with u(y),v(y),T(y)");
}

//////////////////////////////////////////////////////////////////////////////

InitStateTorch::InitStateTorch(const std::string& name) :
  InitState(name),
  _yMax(0.),
  _lookupTableU(),
  _lookupTableV(),
  _lookupTableT()
{
   addConfigOptionsTo(this);

   _datafile = "";
   setParameter("Datafile",&_datafile);
}

//////////////////////////////////////////////////////////////////////////////

InitStateTorch::~InitStateTorch()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateTorch::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitStateTorch::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitStateTorch not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  State dimState;
  std::vector<CFuint>::iterator itd;
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
    State* const currState = states[(*itd)];
    _vFunction.evaluate(currState->getCoordinates(),
			dimState);

    const CFreal yCoord = currState->getCoordinates()[YY];
    if (yCoord <= _yMax) {
      //      (*currState)[1] = _lookupTableU.get(yCoord);
      //  (*currState)[2] = _lookupTableV.get(yCoord);
      dimState[3] = _lookupTableT.get(yCoord);
    }

    if(_inputAdimensionalValues)
    {
      *currState = dimState;
    }
    else
    {
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void InitStateTorch::setup()
{

  InitState::setup();

  // open the given file with the following format
  // nbPoints
  // y, u(y), v(y), T(y)

  ifstream file(_datafile.c_str());

  if (!file) {
    throw Common::FilesystemException
      (FromHere(), "InitStateTorch::setup() => trying to open file");
  }
  
  CFuint nbPoints = 0;
  file >> nbPoints;

  _lookupTableU.reserve(nbPoints);
  _lookupTableV.reserve(nbPoints);
  _lookupTableT.reserve(nbPoints);

  _yMax = 0.0;
  for (CFuint i = 0; i < nbPoints; ++i) {
    CFreal y = 0.0;
    CFreal u = 0.0;
    CFreal v = 0.0;
    CFreal T = 0.0;

    file >> y >> u >> v >> T;

    _lookupTableU.insert(y,u);
    _lookupTableV.insert(y,v);
    _lookupTableT.insert(y,T);

    _yMax = max(_yMax,y);
  }

  _lookupTableU.sortKeys();
  _lookupTableV.sortKeys();
  _lookupTableT.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
