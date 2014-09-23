#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "FiniteVolume/DistanceBasedExtrapolatorGMoveMultiTRS.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "Common/CFPrintContainer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<DistanceBasedExtrapolatorGMoveMultiTRS,
                                  CellCenterFVMData,
		                  NodalStatesExtrapolator<CellCenterFVMData>,
                                  FiniteVolumeModule>
distanceBasedExtrapolatorGMoveMultiTRSProvider("DistanceBasedGMoveMultiTRS");

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveMultiTRS::DistanceBasedExtrapolatorGMoveMultiTRS
(const std::string& name) :
  Framework::DistanceBasedExtrapolator<CellCenterFVMData>(name),
  _trsIDs()
{
  addConfigOptionsTo(this);

  _tv = vector<TrsValuesTuple*>();
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveMultiTRS::~DistanceBasedExtrapolatorGMoveMultiTRS()
{
  for (CFuint i = 0; i < _tv.size(); ++i) {
    deletePtr(_tv[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::extrapolate()
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  
  // reset to 0 the nodal states
  nodalStates[_currNodeID] = 0.;
  applyInner(); 
  
  const int trsID = _trsIDs[_currNodeID];
  if (trsID >= 0) {
    TrsValuesTuple& tvt = *_tv[trsID];
    for (CFuint i = 0; i < tvt._wValuesIdx.size(); ++i) {
      nodalStates[_currNodeID][tvt._wValuesIdx[i]] = tvt._wValues[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::setup()
{
  DistanceBasedExtrapolator<CellCenterFVMData>::setup();
  
  if ( _trsName.empty() ) {
    throw Config::ConfigOptionException (FromHere(),"DistanceBasedExtrapolatorGMoveMultiTRS has List of TRS empty");
  }
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  _trsIDs.resize(nodes.size());
  _trsIDs.assign(_trsIDs.size(), -1);
  const CFuint nbTRSs = _orderedTrsList.size();
  
  CFuint count = 0;
  for (CFuint iName = 0; iName < _trsName.size(); ++iName) {
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = _orderedTrsList[iTRS];
      
      if (currTrs->getName() == _trsName[iName]) {
	Common::SafePtr<std::vector<CFuint> > trsNodes = currTrs->getNodesInTrs();
	const CFuint nbTrsNodes = trsNodes->size();
	for (CFuint i = 0; i < nbTrsNodes; ++i) {
	  // don't override existing TRS IDs >= 0
	  if (_trsIDs[(*trsNodes)[i]] == -1) {
	    _trsIDs[(*trsNodes)[i]] = iName;
	  }
	}
	count++;
	break;
      }
    }
  }
  
  if (count != _trsName.size()) {
    throw BadValueException
      (FromHere(), "DistanceBasedExtrapolatorGMoveMultiTRS::setup() => TRS name not found");
  }

  // for (CFuint i = 0; i < _tv.size(); ++i) {
  //     cout << "Name TRS = " << _trsName[i] << endl;
  //     _tv[i]->sanityCheck();
  //   }
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::configure ( Config::ConfigArgs& args )
{
  DistanceBasedExtrapolator<CellCenterFVMData>::configure(args);

  for (CFuint i = 0; i < _trsName.size(); ++i) {
    TrsValuesTuple* tvt = new TrsValuesTuple(_trsName[i]);
    _tv.push_back(tvt);
    configureNested ( tvt, args );
  }
}

//////////////////////////////////////////////////////////////////////////////

DistanceBasedExtrapolatorGMoveMultiTRS::TrsValuesTuple::TrsValuesTuple
(const std::string& name) :
  Config::ConfigObject(name),
  _idx(0)
{
  addConfigOptionsTo(this);

  _wValuesIdx = vector<CFuint>();
  setParameter("ValuesIdx",&_wValuesIdx);

  _wValues = vector<CFreal>();
  setParameter("Values",&_wValues);
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::TrsValuesTuple::defineConfigOptions
(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("ValuesIdx","Indices of the prescribed values");
  options.addConfigOption< vector<CFreal> >
    ("Values","Prescribed values");
}

//////////////////////////////////////////////////////////////////////////////

void DistanceBasedExtrapolatorGMoveMultiTRS::TrsValuesTuple::sanityCheck()
{
  CFLogInfo(Common::CFPrintContainer<vector<CFuint> >("ValuesIdx  = ", &_wValuesIdx));
  CFLogInfo(Common::CFPrintContainer<vector<CFreal> >("Values  = ", &_wValues));
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
