#include "Framework/CFL.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/PseudoSteadyStdTimeDiagBlockJacob.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< PseudoSteadyStdTimeDiagBlockJacob,
                       SpectralFDMethodData,
                       SpectralFDModule>
BE_DiagBlockJacob("PseudoSteadyTimeDiagBlockJacob");

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyStdTimeDiagBlockJacob::PseudoSteadyStdTimeDiagBlockJacob(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_volumes("volumes",false),
  socket_updateCoeff("updateCoeff"),
//   socket_pastStates("pastStates"),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  m_cellBuilder(CFNULL),
  m_currDiagMatrix(),
  m_nbrEqs(),
  m_cellStates(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_diagValues(),
  m_isUnsteady(),
  m_elemIdx()
{
}

//////////////////////////////////////////////////////////////////////////////

PseudoSteadyStdTimeDiagBlockJacob::~PseudoSteadyStdTimeDiagBlockJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeDiagBlockJacob::setup()
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
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeDiagBlockJacob::execute()
{
  CFAUTOTRACE;

  // get cfl number and time step
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

  // check if computation is unsteady (time accurate)
  m_isUnsteady = dt > 0.;

  // get the diagonal block Jacobian matrix datahandle
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();

  // get datahandle of volumes if necessary
  DataHandle<CFreal> volumes(CFNULL);
  if(m_isUnsteady)
  {
    cf_assert(socket_volumes.isConnected());
    volumes = socket_volumes.getDataHandle();
  }

  // get update coefficients
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // get the local spectral FD data
    vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

    // get solution point local coordinates
    m_solPntsLocalCoords = sdLocalData[iElemType]->getSolPntsLocalCoords();

    // add the diagonal entries in the jacobian (updateCoeff/CFL)
    for (m_elemIdx = startIdx; m_elemIdx < endIdx; ++m_elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = m_elemIdx;
      GeometricEntity* cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = cell->getStates();

      // compute diagonal values
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // set current diagonal block Jacobian matrix
        m_currDiagMatrix = &diagBlockJacobMatr[m_elemIdx];

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
          const CFreal invCellVolume = 1.0/cell->computeVolume();

          // get jacobian determinants at solution points
          const std::valarray<CFreal> jacobDet =
              cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

          const CFuint firstStateID = (*m_cellStates)[0]->getLocalID();
          const CFreal updateCoeffDivCFL = invCellVolume*updateCoeff[firstStateID]/cfl;// should be the same for all states
          for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
          {
            m_diagValues[iSol] = updateCoeffDivCFL*jacobDet[iSol];
          }
        }

        // add the contribution of the time residual to the rhs and the jacobian
        addTimeResidual();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PseudoSteadyStdTimeDiagBlockJacob::addTimeResidual()
{
//   // get and past states
//   DataHandle<State*> pastStates = socket_pastStates.getDataHandle();

  // get factor for the residual (in the Jacobian)
  const CFreal resFactor = getMethodData().getResFactor();

  // multiply diagonal values with the residual factor
  m_diagValues *= resFactor;

  // add time residual contribution (e.g Backward Euler or Cranck-Nicholson, depending on resFactor)
  const CFuint nbrSolPnts = m_cellStates->size();
  CFuint resIdx = 0;
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state
//     const State& currState = *(*m_cellStates)[iSol];

    // get past state
//     const State& pastState = *pastStates[stateID];

    // add contribution to diagonal block jacobian
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resIdx)
    {
      (*m_currDiagMatrix)(resIdx,resIdx) += m_diagValues[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
PseudoSteadyStdTimeDiagBlockJacob::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_volumes);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_diagBlockJacobMatr);
//   result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
