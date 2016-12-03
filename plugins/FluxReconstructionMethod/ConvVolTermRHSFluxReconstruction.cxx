#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvVolTermRHSFluxReconstruction.hh"
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

MethodCommandProvider<ConvVolTermRHSFluxReconstruction, FluxReconstructionSolverData, FluxReconstructionModule> ConvVolTermRHSFluxReconstructionProvider("ConvVolTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvVolTermRHSFluxReconstruction::ConvVolTermRHSFluxReconstruction(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  socket_gradients("gradients"),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_solPntsLocalCoords(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_resUpdates(),
  m_gradUpdates(),
  m_nbrEqs(),
  m_dim()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

ConvVolTermRHSFluxReconstruction::~ConvVolTermRHSFluxReconstruction()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::execute()
{
  CFTRACEBEGIN;

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = false; //getMethodData().hasDiffTerm();

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

    // resize the variable m_resUpdates
    resizeResAndGradUpdates();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // if cell is parallel updatable or the gradients have to be computed,
      // set cell data and reconstruct states
      if ((*m_cellStates)[0]->isParUpdatable() || hasDiffTerm)
      {
        // set the current cell and compute the cell data in the volume term computer
        m_volTermComputer->setCurrentCell(m_cell);
        m_volTermComputer->computeCellData();

        // reconstruct the solution in the flux points
        m_volTermComputer->reconstructStates(*m_cellStates);
      }

      // if cell is parallel updatable, compute the volume term
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // compute the volume term
        m_volTermComputer->computeCellConvVolumeTerm(m_resUpdates);

        // update rhs
        updateRHS();
      }

      // if there is a diffusive term, compute the gradients
      if (hasDiffTerm)
      {
        computeGradients();
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

//   // WRITE GRADIENT TO TECPLOT FILE
//   const CFuint nNodes = MeshDataStack::getActive()->getNbNodes();
//   const CFuint nCells = (*elemType)[0].getNbElems();
//
//   const CFuint nVar = 2;
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
//   // get the gradients
//   DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
//   // get the local spectral FD data
//   vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();
//
//   // loop over element types
//   for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
//   {
//     // get start and end indexes for this type of element
//     const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
//     const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
//
//     // get cell nodes local coordinates
//     SafePtr< vector< RealVector > > cellNodeLocalCoords = sdLocalData[m_iElemType]->getCellNodeCoords();
//
//     // get solution polynomial values at cell nodes
//     vector< vector< CFreal > > solPolyValsAtCellNodes
//                                   = sdLocalData[m_iElemType]->getSolPolyValsAtNode(*cellNodeLocalCoords);
//
//     // get the number of solution points per cell
//     const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();
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
//         for (CFuint iVar = 0; iVar < nVar; ++iVar)
//         {
//           nodeC[iVar][nodeID] = (*(*nodes)[iNode])[iVar];
//
//           for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
//           {
//             const CFuint solID = (*m_cellStates)[iSol]->getLocalID();
//             nodeSol[iVar][nodeID] += solPolyValsAtCellNodes[iNode][iSol]*gradients[solID][1][iVar];
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
//   outFile.open("gradients.plt");
//
//   if (!outFile.is_open())
//   {
//     cerr << "Error: failed to open gradients.plt!!!\nExiting..." << endl;
//     exit(1);
//   }
//
//   // Write header
//   outFile << "TITLE = \"gradients solution data\"" << endl;
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
//     for (CFuint iNode = 0; iNode < 4; ++iNode)
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

void ConvVolTermRHSFluxReconstruction::setVolumeTermData()
{
  // set the volume term data in the volume term computer
  m_volTermComputer->setVolumeTermData(m_iElemType);

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();

  // get solution point local coordinates
  m_solPntsLocalCoords = sdLocalData[m_iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::resizeResAndGradUpdates()
{
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& sdLocalData = getMethodData().getFRLocalData();

  // get the number of solution points in current element type
  const CFuint nbrSolPnts = sdLocalData[m_iElemType]->getNbrOfSolPnts();

  // resize m_resUpdates
  m_resUpdates.resize(nbrSolPnts*m_nbrEqs);

  // resize m_gradUpdates
  m_gradUpdates.resize(0);
  m_gradUpdates.resize(nbrSolPnts);
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_gradUpdates[iSol].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradUpdates[iSol][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::updateRHS()
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
//     if (std::abs(rhs[resID]) > 1e-6)
//     {
//       CF_DEBUG_OBJ(rhs[resID]);
//     }
  }
//   CF_DEBUG_POINT;
//   CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::computeGradients()
{
  // compute the volume term contribution to the gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates);

  // add updates to gradients and divide by Jacobian determinant
  addGradVolTermsAndDivideByJacobDet();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::addGradVolTermsAndDivideByJacobDet()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // get jacobian determinants at solution points
  const std::valarray<CFreal> jacobDet =
      m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_solPntsLocalCoords);

  const CFuint nbrSolPnts = m_gradUpdates.size();
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    // get state ID
    const CFuint solID = (*m_cellStates)[iSol]->getLocalID();

    // inverse Jacobian determinant
    const CFreal invJacobDet = 1.0/jacobDet[iSol];

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[solID][iGrad] += m_gradUpdates[iSol][iGrad];
      gradients[solID][iGrad] *= invJacobDet;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of equations
  m_dim    = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

//   // get the volume term computer
//   m_volTermComputer = getMethodData().getVolTermComputer();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSFluxReconstruction::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvVolTermRHSFluxReconstruction::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
