// #include "Environment/DirPaths.hh"

// #include "Framework/BadFormatException.hh"
#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/LUSGSSetup.hh"

#include "Common/ConnectivityTable.hh"
// #include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< LUSGSSetup,SpectralFDMethodData,SpectralFDModule >
  lusgsSetupProvider("LUSGSSetup");

//////////////////////////////////////////////////////////////////////////////

LUSGSSetup::LUSGSSetup(const std::string& name) :
  StdSetup(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_statesSetStateIDs("statesSetStateIDs"),
  socket_isStatesSetParUpdatable("isStatesSetParUpdatable"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

LUSGSSetup::~LUSGSSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  LUSGSSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > >
      result = StdSetup::needsSockets();

  result.push_back(&socket_states                 );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  LUSGSSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > >
      result = StdSetup::providesSockets();

  result.push_back(&socket_diagBlockJacobMatr     );
  result.push_back(&socket_rhsCurrStatesSet       );
  result.push_back(&socket_statesSetStateIDs      );
  result.push_back(&socket_isStatesSetParUpdatable);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSSetup::execute()
{
  CFAUTOTRACE;

  // call execute() of the parent class
  StdSetup::execute();

  CFLog(INFO,"\n\nLUSGSSetup::execute\n\n");

  // get cell-state connectivity
  SafePtr< ConnectivityTable<CFuint> > cellToStates =
      MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  // number of cells
  const CFuint nbrCells = cellToStates->nbRows();

  // number of equations
  const CFuint nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // maximum number of states in a cell
  CFuint maxNbrStatesInCell = 0;

  // create statesSetStateIDs (will be the same as the cell-state connectivity for SpectralFD),
  // resize diagBlockJacobMatr and set isStatesSetParUpdatable
  DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
  statesSetStateIDs.resize(nbrCells);
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
  diagBlockJacobMatr.resize(nbrCells);
  DataHandle< bool > isStatesSetParUpdatable = socket_isStatesSetParUpdatable.getDataHandle();
  isStatesSetParUpdatable.resize(nbrCells);

  // get state datahandle
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  for (CFuint iCell = 0; iCell < nbrCells; ++iCell)
  {
    const CFuint nbrStatesInCell = cellToStates->nbCols(iCell);
    maxNbrStatesInCell = maxNbrStatesInCell > nbrStatesInCell ?
                         maxNbrStatesInCell : nbrStatesInCell;
    for (CFuint iState = 0; iState < nbrStatesInCell; ++iState)
    {
      statesSetStateIDs[iCell].push_back((*cellToStates)(iCell,iState));
    }
    const CFuint nbrVarsInCell = nbrStatesInCell*nbrEqs;
    diagBlockJacobMatr[iCell].resize(nbrVarsInCell,nbrVarsInCell);

    // set isStatesSetParUpdatable
    const CFuint firstStateID = (*cellToStates)(iCell,0);
    isStatesSetParUpdatable[iCell] = states[firstStateID]->isParUpdatable();
  }

  // maximum number of variables in a cell
  const CFuint maxNbrVarsInCell = nbrEqs*maxNbrStatesInCell;

  // resize rhsCurrStatesSet and updateCoefCurrStatesSet
  DataHandle< CFreal > rhsCurrStatesSet        = socket_rhsCurrStatesSet.getDataHandle();
  rhsCurrStatesSet       .resize(maxNbrVarsInCell);
}

/////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD
