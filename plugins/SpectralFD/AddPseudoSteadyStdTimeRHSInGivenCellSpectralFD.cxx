#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD,
                       SpectralFDMethodData,
                       SpectralFDModule>
AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD("PseudoSteadyTimeRHSInGivenCell");

//////////////////////////////////////////////////////////////////////////////

AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_volumes("volumes",false),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_updateCoeff("updateCoeff"),
  socket_statesSetIdx("statesSetIdx"),
  socket_pastStates("pastStates"),
  m_cellBuilder(CFNULL),
  m_cell(CFNULL),
  m_nbrEqs(),
  m_cellStates(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_diagValues(),
  m_isUnsteady(),
  m_iElemType(),
  m_currIter()
{
}

//////////////////////////////////////////////////////////////////////////////

AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::~AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::setup()
{
  SpectralFDMethodCom::setup();

  // get number of equations in the physical model
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get maximum number of solution points in a cell
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  const CFuint nbrElemTypes = sdLocalData.size();
  CFuint maxNbrSolPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrSolPnts = sdLocalData[iElemType]->getNbrOfSolPnts();
    maxNbrSolPnts = maxNbrSolPnts > nbrSolPnts ? maxNbrSolPnts : nbrSolPnts;
  }

  // resize m_diagValues
  m_diagValues.resize(maxNbrSolPnts);

  // initialize to 10, to make sure that solution point local coordinates are set at first iteration
  m_currIter = 10;
}

//////////////////////////////////////////////////////////////////////////////

void AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::execute()
{
  CFAUTOTRACE;

  // update current iteration number
  const CFuint currIter = SubSystemStatusStack::getActive()->getNbIter();
  if (m_currIter != currIter)
  // if we are at a new time step, then the diagonal block Jacobian matrices could
  // have been computed, and the data in the volume and face term computers should be reset.
  {
    // update m_currIter
    m_currIter = currIter;

    // get the local spectral FD data
    vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

    // get solution point local coordinates
    m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
  }

  // get cfl number and time step
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal dt  = SubSystemStatusStack::getActive()->getDT();

  // check if computation is unsteady (time accurate)
  m_isUnsteady = dt > 0.;

  // get datahandle of volumes if necessary
  DataHandle<CFreal> volumes(CFNULL);
  if(m_isUnsteady)
  {
    cf_assert(socket_volumes.isConnected());
    volumes = socket_volumes.getDataHandle();
  }

  // get index of cell for which to compute the rhs
  DataHandle< CFint > statesSetIdx = socket_statesSetIdx.getDataHandle();
  const CFuint cellIdx = statesSetIdx[0];

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the cell builder and set the TRS and the cell index
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  geoData.idx = cellIdx;

  // build the cell
  m_cell = m_cellBuilder->buildGE();

  // get the states in this cell
  m_cellStates = m_cell->getStates();

  // if cell is parallel updatable, compute the rhs
  if ((*m_cellStates)[0]->isParUpdatable())
  {
    /// @note KVDA: for now, it is assumed that there is only one element type in the mesh.
    /// Could pass the element type from LUSGSMethod, in the same socket as the statesSetIdx.
    const CFuint currElemType = 0;
    if (m_iElemType != currElemType)
    {
      m_iElemType = currElemType;

      // get the local spectral FD data
      vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

      // get solution point local coordinates
      m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
    }

    // compute diagonal values (jacobDet/dt)
    const CFuint nbrSolPnts = m_cellStates->size();
    cf_assert(m_diagValues.size() >= nbrSolPnts);
    if (m_isUnsteady)
    {
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        m_diagValues[iSol] = volumes[(*m_cellStates)[iSol]->getLocalID()]/dt;
      }
    }
    else
    {
      // get the cell volume
      const CFreal invCellVolume = 1.0/m_cell->computeVolume();

      // get jacobian determinants at solution points
      const std::valarray<CFreal> jacobDet =
          m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

      DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
      const CFuint firstStateID = (*m_cellStates)[0]->getLocalID();
      const CFreal updateCoeffDivCFL = invCellVolume*updateCoeff[firstStateID]/cfl;// should be the same for all states
      for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
      {
        m_diagValues[iSol] = updateCoeffDivCFL*jacobDet[iSol];
      }
    }

    // add the contribution of the time residual to the current cell rhs
    addTimeResidual();
  }

  //release the GeometricEntity
  m_cellBuilder->releaseGE();
}

//////////////////////////////////////////////////////////////////////////////

void AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::addTimeResidual()
{
  // get factor for the residual
  const CFreal resFactor = getMethodData().getResFactor();

  // multiply diagonal values with the residual factor
  m_diagValues *= resFactor;

  // get current cell rhs and past states
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();
  DataHandle< State* > pastStates       = socket_pastStates      .getDataHandle();

  // add time residual contribution
  const CFuint nbrSolPnts = m_cellStates->size();
  CFuint resIdx = 0;
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
    const State& currState = *(*m_cellStates)[iSol];

    // get state ID
    const CFuint stateID = currState.getLocalID();

    // get past state
    const State& pastState = *pastStates[stateID];

    // add contribution to rhs
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resIdx)
    {
      rhsCurrStatesSet[resIdx] -= (currState[iEq] - pastState[iEq])*m_diagValues[iSol];
//       CF_DEBUG_OBJ(rhsCurrStatesSet[resIdx]);
    }
  }
//   CF_DEBUG_POINT;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
AddPseudoSteadyStdTimeRHSInGivenCellSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_volumes);
  result.push_back(&socket_rhsCurrStatesSet);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_statesSetIdx);
  result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
