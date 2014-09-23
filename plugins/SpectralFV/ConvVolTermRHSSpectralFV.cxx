#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/ConvVolTermRHSSpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"

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

MethodCommandProvider<ConvVolTermRHSSpectralFV, SpectralFVMethodData, SpectralFVModule> ConvVolTermRHSSpectralFVProvider("ConvVolTermRHS");

//////////////////////////////////////////////////////////////////////////////

ConvVolTermRHSSpectralFV::ConvVolTermRHSSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_rhs("rhs"),
  socket_gradients("gradients"),
  m_cellBuilder(CFNULL),
  m_volTermComputer(CFNULL),
  m_invVolFracCVs(CFNULL),
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

ConvVolTermRHSSpectralFV::~ConvVolTermRHSSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::configure ( Config::ConfigArgs& args )
{
  SpectralFVMethodCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::execute()
{
  CFTRACEBEGIN;

  // boolean telling whether there is a diffusive term
  const bool hasDiffTerm = getMethodData().hasDiffTerm();

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
/*CF_DEBUG_EXIT;*/

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
//         for (CFuint iVar = 0; iVar < nVar; ++iVar)
//         {
//           nodeC[iVar][nodeID] = (*(*nodes)[iNode])[iVar];
//
//           for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
//           {
//             const CFuint cvID = (*m_cellStates)[iCV]->getLocalID();
//             nodeSol[iVar][nodeID] += svPolyValsAtSVNodes[iNode][iCV]*gradients[cvID][1][iVar];
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

void ConvVolTermRHSSpectralFV::setVolumeTermData()
{
  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the inverse volume fractions of the CVs
  m_invVolFracCVs = svLocalData[m_iElemType]->getInvVolFracCV();

  // set the volume term data in the volume term computer
  m_volTermComputer->setVolumeTermData(m_iElemType);
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::resizeResAndGradUpdates()
{
  // get the local spectral FV data
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();

  // get the number of CVs in current element type
  const CFuint nbrCVs = svLocalData[m_iElemType]->getNbrOfCVs();

  // resize m_resUpdates
  m_resUpdates.resize(nbrCVs*m_nbrEqs);

  // resize m_gradUpdates
  m_gradUpdates.resize(0);
  m_gradUpdates.resize(nbrCVs);
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    m_gradUpdates[iCV].resize(m_nbrEqs);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_gradUpdates[iCV][iEq].resize(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::updateRHS()
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

void ConvVolTermRHSSpectralFV::computeGradients()
{
  // compute the volume term contribution to the gradients
  m_volTermComputer->computeGradientVolumeTerm(m_gradUpdates);

  // add updates to gradients and divide by CV volume
  addGradVolTermsAndDivideByCellVolume();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::addGradVolTermsAndDivideByCellVolume()
{
  // get the gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // get the cell volume
  const CFreal invCellVolume = 1.0/m_cell->computeVolume();

  const CFuint nbrCVs = m_gradUpdates.size();
  for (CFuint iCV = 0; iCV < nbrCVs; ++iCV)
  {
    // compute inverse CV volume
    const CFreal invVolume = (*m_invVolFracCVs)[iCV]*invCellVolume;

    // get CV ID
    const CFuint cvID = (*m_cellStates)[iCV]->getLocalID();

    // update gradients
    for (CFuint iGrad = 0; iGrad < m_nbrEqs; ++iGrad)
    {
      gradients[cvID][iGrad] += m_gradUpdates[iCV][iGrad];
      gradients[cvID][iGrad] *= invVolume;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::setup()
{
  CFAUTOTRACE;

  // dimensionality and number of equations
  m_dim    = PhysicalModelStack::getActive()->getDim();
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the volume term computer
  m_volTermComputer = getMethodData().getVolTermComputer();
}

//////////////////////////////////////////////////////////////////////////////

void ConvVolTermRHSSpectralFV::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> >
ConvVolTermRHSSpectralFV::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
