#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BasePhysicality.hh"
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

void BasePhysicality::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

BasePhysicality::BasePhysicality(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_solPntsLocalCoords(),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_order(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_cellStatesFlxPnt()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

BasePhysicality::~BasePhysicality()
{
}

//////////////////////////////////////////////////////////////////////////////

void BasePhysicality::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BasePhysicality::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  
  const CFuint nbrElemTypes = elemType->size();
  
  CFuint nbLimits = 0;
  
  // LOOP OVER ELEMENT TYPES, TO SET THE MINIMUM AND MAXIMUM CELL AVERAGED STATES IN THE NODES
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // loop over cells, to put the minimum and maximum neighbouring cell averaged states in the nodes
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      
      computeFlxPntStates(m_cellStatesFlxPnt);

      if (!checkPhysicality())
      {
	enforcePhysicality();
	
	++nbLimits;
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  CFLog(NOTICE, "Number of times physicality enforced: " << nbLimits << "\n");
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void BasePhysicality::computeFlxPntStates(std::vector< RealVector >& statesFlxPnt)
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_maxNbrFlxPnts; ++iFlxPnt)
  {        
    statesFlxPnt[iFlxPnt] = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      statesFlxPnt[iFlxPnt] += (*m_solPolyValsAtFlxPnts)[iFlxPnt][iSol]*(*((*m_cellStates)[iSol]));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePhysicality::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  
  // get nbr of sol pnts
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get the maximum number of flux points
  m_maxNbrFlxPnts = 0;
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    m_maxNbrFlxPnts = m_maxNbrFlxPnts > nbrFlxPnts ? m_maxNbrFlxPnts : nbrFlxPnts;
  }
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    RealVector temp(m_nbrEqs);
    m_cellStatesFlxPnt.push_back(temp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePhysicality::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
