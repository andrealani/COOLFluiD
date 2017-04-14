#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/TVBLimiter.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<TVBLimiter, FluxReconstructionSolverData, FluxReconstructionModule> TVBLimiterFluxReconstructionProvider("TVBLimiter");

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("TVBFactor","Factor used to determine whether to limit the solution or not.");
}

//////////////////////////////////////////////////////////////////////////////

TVBLimiter::TVBLimiter(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_nodeNghbCellMinAvgStates("nodeNghbCellMinAvgStates"),
  socket_nodeNghbCellMaxAvgStates("nodeNghbCellMaxAvgStates"),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_cellNodes(),
  m_cellAvgSolCoefs(),
  m_cellCenterDerivCoefs(),
  m_flxPntsRecCoefs(),
  m_allFlxPntIdxs(),
  m_solPntsLocalCoords(),
  m_cellAvgState(),
  m_cellCenterDerivVar(),
  m_minAvgState(),
  m_maxAvgState(),
  m_minAvgStateAll(),
  m_maxAvgStateAll(),
  m_nbrEqs(),
  m_nbrFlxPnts(),
  m_nbrSolPnts(),
  m_nbrCornerNodes(),
  m_applyLimiter(),
  m_applyLimiterToPhysVar(),
  m_tvbLimitFactor(),
  m_lengthScaleExp()
{
  addConfigOptionsTo(this);

  m_tvbLimitFactor = 0.0;
  setParameter( "TVBFactor", &m_tvbLimitFactor );
}

//////////////////////////////////////////////////////////////////////////////

TVBLimiter::~TVBLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // RESET THE MINIMUM AND MAXIMUM CELL AVERAGED SOLUTIONS (STORED IN THE NODES!!!)
  resetNodeNghbrCellAvgStates();

  // LOOP OVER ELEMENT TYPES, TO SET THE MINIMUM AND MAXIMUM CELL AVERAGED STATES IN THE NODES
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // set reconstruction data for averaged solutions
    setAvgReconstructionData(iElemType);

    // loop over cells, to put the minimum and maximum neighbouring cell averaged states in the nodes
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes ();

      // reconstruct cell averaged state
      reconstructCellAveragedState();

      // set the min and max node neighbour cell averaged states
      setMinMaxNodeNghbrCellAvgStates();

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  // LOOP OVER ELEMENT TYPES, TO LIMIT THE SOLUTIONS
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // set reconstruction data
    setAllReconstructionData(iElemType);

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes ();

      // set the min and max neighbouring cell averaged states
      setMinMaxNghbrCellAvgStates();

      // check if limiting is necessary
      setLimitBooleans();

      // apply limiter if necessary
      if (m_applyLimiter)
      {
        limitStates();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::resetNodeNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // get the number of nodes
  const CFuint nbrNodes = nodeNghbCellMinAvgStates.size();
  cf_assert(nbrNodes == nodeNghbCellMaxAvgStates.size());

  // reset the cell averaged states
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodeNghbCellMinAvgStates[iNode] = +MathTools::MathConsts::CFrealMax();
    nodeNghbCellMaxAvgStates[iNode] = -MathTools::MathConsts::CFrealMax();
  }
  m_minAvgStateAll = +MathTools::MathConsts::CFrealMax();
  m_maxAvgStateAll = -MathTools::MathConsts::CFrealMax();
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setAvgReconstructionData(CFuint iElemType)
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // number of cell corner nodes
  /// @note in the future, hanging nodes should be taken into account here
  m_nbrCornerNodes = frLocalData[iElemType]->getNbrCornerNodes();

    // get solution point local coordinates
  m_cellAvgSolCoefs = frLocalData[iElemType]->getCellAvgSolCoefs();

    // number of solution points
  m_nbrSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setAllReconstructionData(CFuint iElemType)
{
  // set reconstruction data for averaged solutions
  setAvgReconstructionData(iElemType);

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get coefficients for cell center derivatives
  m_cellCenterDerivCoefs = frLocalData[iElemType]->getCellCenterDerivCoefs();

  // get the mapped coordinates of the solution points
  m_solPntsLocalCoords = frLocalData[iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::reconstructCellAveragedState()
{
  m_cellAvgState = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0]);
  for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::reconstructCellAveragedVariable(const CFuint iEq)
{
  m_cellAvgState[iEq] = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0])[iEq];
  for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState[iEq] += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol])[iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::computeCellCenterDerivVariable(const CFuint iEq)
{
  for (CFuint iDeriv = 0; iDeriv < m_dim; ++iDeriv)
  {
    m_cellCenterDerivVar[iDeriv] = (*m_cellCenterDerivCoefs)[iDeriv][0]*(*(*m_cellStates)[0])[iEq];
    for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellCenterDerivVar[iDeriv] += (*m_cellCenterDerivCoefs)[iDeriv][iSol]*(*(*m_cellStates)[iSol])[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setMinMaxNodeNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // loop over corner nodes, to set the neighbouring min and max states
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    // get min and max states
    RealVector& minAvgState = nodeNghbCellMinAvgStates[nodeID];
    RealVector& maxAvgState = nodeNghbCellMaxAvgStates[nodeID];

    // loop over physical variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // minimum average state
      minAvgState[iEq] = m_cellAvgState[iEq] < minAvgState[iEq] ? m_cellAvgState[iEq] : minAvgState[iEq];

      // maximum average state
      maxAvgState[iEq] = m_cellAvgState[iEq] > maxAvgState[iEq] ? m_cellAvgState[iEq] : maxAvgState[iEq];
    }
  }

  // set min and max of all states on the mesh
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
      // minimum average state
    m_minAvgStateAll[iEq] = m_cellAvgState[iEq] < m_minAvgStateAll[iEq] ? m_cellAvgState[iEq] : m_minAvgStateAll[iEq];

      // maximum average state
    m_maxAvgStateAll[iEq] = m_cellAvgState[iEq] > m_maxAvgStateAll[iEq] ? m_cellAvgState[iEq] : m_maxAvgStateAll[iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setMinMaxNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // set the minimum and maximum average states in the node neighbours
  const CFuint nodeID = (*m_cellNodes)[0]->getLocalID();
  m_minAvgState = nodeNghbCellMinAvgStates[nodeID];
  m_minAvgState = nodeNghbCellMaxAvgStates[nodeID];
  for (CFuint iNode = 1; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    // get min and max states in nodes
    const RealVector& noEminAvgState = nodeNghbCellMinAvgStates[nodeID];
    const RealVector& nodeMaxAvgState = nodeNghbCellMaxAvgStates[nodeID];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_minAvgState[iEq] = m_minAvgState[iEq] < noEminAvgState[iEq] ? m_minAvgState[iEq] : noEminAvgState[iEq];
      m_maxAvgState[iEq] = m_maxAvgState[iEq] > nodeMaxAvgState[iEq] ? m_maxAvgState[iEq] : nodeMaxAvgState[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setLimitBooleans()
{
  // compute length scale factor
  const CFreal lengthScaleFactor = m_tvbLimitFactor*pow(m_cell->computeVolume(),m_lengthScaleExp);

  // reset m_applyLimiter
  m_applyLimiter = false;

  // check if limiting is necessary
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_applyLimiterToPhysVar[iEq] = false;

    // compute upper and lower boundaries
    const CFreal dVar = lengthScaleFactor*(m_maxAvgStateAll[iEq] - m_minAvgStateAll[iEq]);
    const CFreal lowerBnd = m_minAvgState[iEq] - dVar;
    const CFreal upperBnd = m_maxAvgState[iEq] + dVar;
    for (CFuint iSol = 0; iSol < m_nbrSolPnts && !m_applyLimiterToPhysVar[iEq]; ++iSol)
    {
      m_applyLimiterToPhysVar[iEq] = (*((*m_cellStates)[iSol]))[iEq] < lowerBnd ||
                                     (*((*m_cellStates)[iSol]))[iEq] > upperBnd;
    }

    if (m_applyLimiterToPhysVar[iEq])
    {
      m_applyLimiter = true;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::limitStates()
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
/*    if (m_applyLimiterToPhysVar[iEq])
    {*/
      // reconstruct the cell averaged variable and the derivatives in the cell center
      reconstructCellAveragedVariable(iEq);
      computeCellCenterDerivVariable (iEq);

      // limit the gradient to ensure that the limited solution is within the boundaries [UAvgMin,UAvgMax]
      // U = Uavg + dUdksi*ksi + dUdeta*eta + dUdzta*zta, with ksi, eta and zta in [-1,+1]
      // ==> Umax = Uavg + abs(dUdksi) + abs(dUdeta) + abs(dUdzta)
      // ==> Umin = Uavg - abs(dUdksi) - abs(dUdeta) - abs(dUdzta)
      CFreal dSolCellMax = std::abs(m_cellCenterDerivVar[KSI]);
      for (CFuint iCoor = 1; iCoor < m_dim; ++iCoor)
      {
        dSolCellMax += std::abs(m_cellCenterDerivVar[iCoor]);
      }
      if (dSolCellMax > MathTools::MathConsts::CFrealMin())
      {
        CFreal limitFact = 1.0;
        if (m_cellAvgState[iEq] + dSolCellMax > m_maxAvgState[iEq])
        {
          limitFact = (m_maxAvgState[iEq]-m_cellAvgState[iEq])/dSolCellMax;
        }
        if (m_cellAvgState[iEq] - dSolCellMax < m_minAvgState[iEq])
        {
          const CFreal limitFact2 = (m_cellAvgState[iEq]-m_minAvgState[iEq])/dSolCellMax;
          limitFact = limitFact2 < limitFact ? limitFact2 : limitFact;
        }
        m_cellCenterDerivVar *= limitFact;
      }

      // set the solution in the solution points
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        (*(*m_cellStates)[iSol])[iEq] = m_cellAvgState[iEq];
        for (CFuint iCoor = 0; iCoor < m_dim; ++iCoor)
        {
          (*(*m_cellStates)[iSol])[iEq] +=
              m_cellCenterDerivVar[iCoor]*(*m_solPntsLocalCoords)[iSol][iCoor];
        }
      }
//     }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_applyLimiterToPhysVar.resize(m_nbrEqs);
  m_cellAvgState.resize(m_nbrEqs);
  m_cellCenterDerivVar.resize(m_dim);
  m_minAvgState.resize(m_nbrEqs);
  m_maxAvgState.resize(m_nbrEqs);
  m_minAvgStateAll.resize(m_nbrEqs);
  m_maxAvgStateAll.resize(m_nbrEqs);

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);

  // get the maximum number of flux points
  CFuint maxNbrFlxPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    maxNbrFlxPnts = maxNbrFlxPnts > nbrFlxPnts ? maxNbrFlxPnts : nbrFlxPnts;
  }

  // get the number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();

  // resize the sockets for the minimum and maximum states in the node neighbouring cells
  /// @todo this socket is too big, only the first order nodes (in the cell corners, and later on, hanging nodes)
  /// have to be taken into account
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
  nodeNghbCellMinAvgStates.resize(nbrNodes);
  nodeNghbCellMaxAvgStates.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodeNghbCellMinAvgStates[iNode].resize(m_nbrEqs);
    nodeNghbCellMaxAvgStates[iNode].resize(m_nbrEqs);
  }

  m_lengthScaleExp = 2.0/static_cast<CFreal>(m_dim);
}

//////////////////////////////////////////////////////////////////////////////

void TVBLimiter::unsetup()
{
  CFAUTOTRACE;

  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
  nodeNghbCellMinAvgStates.resize(0);
  nodeNghbCellMaxAvgStates.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> >
TVBLimiter::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nodeNghbCellMinAvgStates           );
  result.push_back(&socket_nodeNghbCellMaxAvgStates           );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
