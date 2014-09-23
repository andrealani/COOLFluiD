#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFD/DiffVolTermRHSSpectralFD.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DiffVolTermRHSSpectralFD, SpectralFDMethodData, SpectralFDModule> DiffVolTermRHSSpectralFDProvider("DiffVolTermRHS");

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSSpectralFD::DiffVolTermRHSSpectralFD(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_rhs("rhs"),
  socket_gradients("gradients"),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_cellGrads(CFNULL),
  m_resUpdates(),
  m_nbrEqs()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSSpectralFD::~DiffVolTermRHSSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::configure ( Config::ConfigArgs& args )
{
  SpectralFDMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // loop over element types
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    // set the volume term data for this element type
    setVolumeTermData();

    // resize the variables m_resUpdates and m_cellGrads
    resizeResUpdatesAndCellGrads();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the gradients in this cell (in variable m_cellGrads)
      setGradients();

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);

        // reconstruct the gradients in the flux points
        m_volTermComputer->reconstructGradients(m_cellGrads);

        // compute the volume term
        m_volTermComputer->computeCellDiffVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::setVolumeTermData()
{
  // set the volume term data in the volume term computer
  m_volTermComputer->setVolumeTermData(m_iElemType);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::resizeResUpdatesAndCellGrads()
{
  // get the local spectral FD data
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // resize m_resUpdates
  m_resUpdates.resize(nbrSolPnts*m_nbrEqs);

  // resize m_cellGrads
  m_cellGrads.resize(nbrSolPnts);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::setGradients()
{
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // get number of states
  const CFuint nbrStates = m_cellStates->size();

  // set gradients
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    m_cellGrads[iState] = &gradients[stateID];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::updateRHS()
{
  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // update rhs
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrRes = m_resUpdates.size();
  for (CFuint iRes = 0; iRes < nbrRes; ++iRes, ++resID)
  {
      rhs[resID] += resFactor*m_resUpdates[iRes];
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::setup()
{
  CFAUTOTRACE;

  // get the number of equations in the physical modes
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the volume term computer
  m_volTermComputer = getMethodData().getVolTermComputer();
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFD::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DiffVolTermRHSSpectralFD::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
