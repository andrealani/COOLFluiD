#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/DiffVolTermRHSSpectralFV.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DiffVolTermRHSSpectralFV, SpectralFVMethodData, SpectralFVModule> DiffVolTermRHSSpectralFVProvider("DiffVolTermRHS");

//////////////////////////////////////////////////////////////////////////////

DiffVolTermRHSSpectralFV::DiffVolTermRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
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

DiffVolTermRHSSpectralFV::~DiffVolTermRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::execute()
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

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // get the gradients in this cell (in variable m_cellGrads)
        setGradients();

        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);

        // reconstruct the gradients in the flux points
        m_volTermComputer->reconstructGradients(m_cellGrads,m_cellGrads.size());

        // compute the volume term
        m_volTermComputer->computeCellDiffVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

//   // WRITE RHS TO TECPLOT FILE
//   const CFuint nNodes = MeshDataStack::getActive()->getNbNodes();
//   const CFuint nCells = (*elemType)[0].getNbElems();
//
//   const CFuint nVar = 4;
//   CFreal** nodeC;
//   CFreal** nodeSol;
//   CFuint*  nodeNbCells;
//
//   // allocate memory
//   nodeSol = new CFreal*[nVar];
//   nodeC   = new CFreal*[nVar];
//   for (CFuint iVar = 0; iVar < nVar; ++iVar)
//   {
//     nodeSol[iVar] = new CFreal[nNodes];
//     nodeC[iVar] = new CFreal[nNodes];
//   }
//
//   nodeNbCells = new CFuint[nNodes];
//
//   // initialize solution in nodes and number of neighbouring cells
//   for (CFuint nodeID = 0; nodeID < nNodes; ++nodeID)
//   {
//     nodeNbCells[nodeID] = 0;
//
//     for (CFuint iVar = 0; iVar < nVar; ++iVar)
//     {
//       nodeSol[iVar][nodeID] = 0.0;
//     }
//   }
//
//   // get the datahandle of the rhs
//   DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
//   // get the local spectral FV data
//   vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
//
//   // loop over element types
//   for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
//   {
//     // get start and end indexes for this type of element
//     const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
//     const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
//
//     // get SV nodes local coordinates
//     SafePtr< vector< RealVector > > svNodeLocalCoords = svLocalData[m_iElemType]->getSVNodeCoords();
//
//     // get SV polynomial values at SV nodes
//     vector< vector< CFreal > > svPolyValsAtSVNodes
//                                   = svLocalData[m_iElemType]->getSVPolyValsAtNode(*svNodeLocalCoords);
//
//     // get the number of control volumes per spectral volume
//     const CFuint nbrCVs = svLocalData[m_iElemType]->getNbrOfCVs();
//
//     // loop over cells
//     for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
//     {
//       // build the GeometricEntity
//       geoData.idx = elemIdx;
//       m_cell = m_cellBuilder->buildGE();
//
//       // get the states in this cell
//       m_cellStates = m_cell->getStates();
//
//       vector< Node* >* nodes = m_cell->getNodes();
//       const CFuint nCellNodes = nodes->size();
//
//       // loop over nodes in cell
//       for (CFuint iNode = 0; iNode < nCellNodes; ++iNode)
//       {
//         const CFuint nodeID = (*nodes)[iNode]->getLocalID();
//         ++nodeNbCells[nodeID];
//
//         nodeC[0][nodeID] = (*(*nodes)[iNode])[0];
//         nodeC[1][nodeID] = (*(*nodes)[iNode])[1];
//
//         for (CFuint iVar = 0; iVar < nVar; ++iVar)
//         {
//           for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
//           {
//             const CFuint cvID = (*m_cellStates)[iCV]->getLocalID();
//             nodeSol[iVar][nodeID] += svPolyValsAtSVNodes[iNode][iCV]*rhs(cvID,iVar,nVar);
//           }
//         }
//       }
//
//       //release the GeometricEntity
//       m_cellBuilder->releaseGE();
//     }
//   }
//
//   // Open solution file
//   ofstream outFile;
//   outFile.open("rhs.plt");
//
//   if (!outFile.is_open())
//   {
//     cerr << "Error: failed to open rhs.plt!!!\nExiting..." << endl;
//     exit(1);
//   }
//
//   // Write header
//   outFile << "TITLE = \"rhs data\"" << endl;
//   outFile << "VARIABLES = \"x\", \"y\"";
//   for (CFuint iVar = 0; iVar < nVar; iVar++)
//     outFile << ", \"" << iVar+1 << "\"";
//   outFile << endl;
//   outFile << endl;
//
//   // Write zone for whole domain
//   outFile.precision(15);
//   outFile << "ZONE N=" << nNodes << ", E=" << nCells << ", F=FEPOINT, ET=TRIANGLE" << endl;
//   outFile << endl;
//
//   // Write solution in nodes
//   for (CFuint nodeID = 0; nodeID < nNodes; ++nodeID)
//   {
//     // Write coordinates to file
//     outFile << nodeC[0][nodeID] << " " << nodeC[1][nodeID] << " ";
//
//     // Write solution to file
//     for (CFuint iVar = 0; iVar < nVar; ++iVar)
//     {
//       outFile << nodeSol[iVar][nodeID]/nodeNbCells[nodeID] << " ";
//     }
//
//     outFile << endl;
//   }
//   outFile << endl;
//
//   // Write cell node connectivity to file
//   SafePtr< ConnectivityTable<CFuint> > cellNode = MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
//   for (CFuint cellID = 0; cellID < nCells; ++cellID)
//   {
//     for (CFuint iNode = 0; iNode < 3; ++iNode)
//     {
//       outFile << (*cellNode)(cellID,iNode)+1 << " ";
//     }
//     outFile << endl;
//   }
//
//   // Close the file
//   outFile.close();
//
//   // Free memory
//   for (CFuint iVar = 0; iVar < nVar; ++iVar)
//   {
//     delete[] nodeSol[iVar];
//     delete[] nodeC[iVar];
//   }
//   delete[] nodeSol;
//   delete[] nodeC;
//
//   delete[] nodeNbCells;
// CF_DEBUG_EXIT;

  CFTRACEEND;
 }

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::setVolumeTermData()
{
  // set the volume term data in the volume term computer
  m_volTermComputer->setVolumeTermData(m_iElemType);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::resizeResUpdatesAndCellGrads()
{
  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the number of CVs in current element type
  const CFuint nbrCVs = svLocalData[m_iElemType]->getNbrOfCVs();

  // resize m_resUpdates
  m_resUpdates.resize(nbrCVs*m_nbrEqs);

  // resize m_cellGrads
  m_cellGrads.resize(nbrCVs);
}

//////////////////////////////////////////////////////////////////////////////

void DiffVolTermRHSSpectralFV::setGradients()
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

void DiffVolTermRHSSpectralFV::updateRHS()
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

void DiffVolTermRHSSpectralFV::setup()
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

void DiffVolTermRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
DiffVolTermRHSSpectralFV::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
