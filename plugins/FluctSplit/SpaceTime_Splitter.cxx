#include "Framework/State.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

SpaceTime_Splitter::SpaceTime_Splitter(const std::string& name):
  Splitter(name),
  _cellVolume(),
  _pastCellVolume(),
  _nodeArea(),
  _consStates(),
  _cellSpeed(),
  _timeStep(),
  _updateCoeff(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

SpaceTime_Splitter::~SpaceTime_Splitter()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime_Splitter::setup()
{
  Splitter::setup();

  _cellSpeed.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime_Splitter::setConsStates(const vector<State*>& p1)
{
  cf_assert(p1.size()!=0);

  ///@todo this will have to be changed to make it safer
  _consStates.resize(p1.size());

  for (CFuint i=0;i < p1.size();++i){
    _consStates[i] = p1[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime_Splitter::setInterStates(const vector<State*>& p1)
{
  cf_assert(p1.size()!=0);

  ///@todo this will have to be changed to make it safer
  _interStates.resize(p1.size());

  for (CFuint i=0;i<p1.size();++i){
    _interStates[i] = p1[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void SpaceTime_Splitter::setInterConsStates(const vector<State*>& p1)
{
  cf_assert(p1.size()!=0);

  ///@todo this will have to be changed to make it safer
  _interConsStates.resize(p1.size());

  for (CFuint i=0;i<p1.size();++i){
    _interConsStates[i] = p1[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

