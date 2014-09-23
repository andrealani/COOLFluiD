#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "Common/PE.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Common/OSystem.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PathAppender.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "TecplotWriterNavierStokes/TecplotWriterNavierStokes.hh"
#include "TecplotWriterNavierStokes/WriteInstantAndAvgSolution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteInstantAndAvgSolution, TecWriterData, TecplotWriterNavierStokesModule>
writeInstantAndAvgSolutionProvider("WriteInstantAndAvgSolution");

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolution::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("UpdateVar","Update variables used in the space method.");
   options.addConfigOption< CFuint >("InitAvgStepsCounter","Initial value for counter of averaging steps.");
   options.addConfigOption< CFuint >("WriteToFileRate","Rate at which the solutions are written to a file.");
   options.addConfigOption< bool >("AddInitSolToAvgSol","Boolean telling whether to add the initial solution to the averaged solution.");
   options.addConfigOption< std::string >("AvgSolFileName","Name of the file where the averaged solution will be stored.");
}

//////////////////////////////////////////////////////////////////////////////

WriteInstantAndAvgSolution::WriteInstantAndAvgSolution(const std::string& name) : WriteSolution(name),
  socket_nodeAvgVars("nodeAvgVars"),
  m_updateToPrimVar(),
  m_updateVarStr(),
  m_primVarNames(),
  m_extraAvgVarNames(),
  m_nodalPrimState(CFNULL),
  m_nodeAvgVals(),
  m_avgStepsCounter(),
  m_writeToFileRate(),
  m_addCurrSolToAvgSol(),
  m_avgSolFileNameStr(),
  m_avgSolFileName()
{
  addConfigOptionsTo(this);

  m_updateVarStr = "Cons";
  setParameter("UpdateVar",&m_updateVarStr);

  m_avgStepsCounter = 0;
  setParameter("InitAvgStepsCounter",&m_avgStepsCounter);

  m_writeToFileRate = 1;
  setParameter("WriteToFileRate",&m_writeToFileRate);

  m_addCurrSolToAvgSol = false; // by default do not add the initial solution to the average
  setParameter("AddInitSolToAvgSol",&m_addCurrSolToAvgSol);

  m_avgSolFileNameStr = "avgSol-statistics.plt";
  setParameter("AvgSolFileName",&m_avgSolFileNameStr);
}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolution::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

  // get the nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // this is a sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
      socket_nstatesProxy.getDataHandle();
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // averaged variables in nodes
  DataHandle< RealVector> nodeAvgVars = socket_nodeAvgVars.getDataHandle();

  // get dimensionality, number of equations
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal nbEqsM1 = nbEqs-1;

  if (m_addCurrSolToAvgSol)
  {
    // factors for averaging
    const CFreal avgStepsCounterReal = static_cast<CFreal>(m_avgStepsCounter);
    const CFreal oldAvgSolWght = avgStepsCounterReal/(avgStepsCounterReal+1.0);
    const CFreal curInsSolWght = 1.0/(avgStepsCounterReal+1.0);

    // array to store the dimensional states
    RealVector dimState(nbEqs);
    State tempState;

    // number of nodes
    const CFuint nbrNodes = nodes.size();
    for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
    {
      // copy current state to tempState
      const RealVector& currState = *nodalStates.getState(iNode);
      for (CFuint ieq = 0; ieq < nbEqs; ++ieq)
      {
        tempState[ieq] = currState[ieq];
      }

      // convert to primitive variables
      *m_nodalPrimState = *m_updateToPrimVar->transform(&tempState);

      // make dimensional
      SafePtr< EulerVarSet > eulerVarSet = getMethodData().getUpdateVarSet().d_castTo<EulerVarSet>();
      const RealVector& refData = eulerVarSet->getModel()->getReferencePhysicalData();
      dimState[0      ] = refData[EulerTerm::RHO]*(*m_nodalPrimState)[0];
      dimState[nbEqsM1] = refData[EulerTerm::P]*(*m_nodalPrimState)[nbEqsM1];
      for (CFuint iDim = 0; iDim < dim; ++iDim)
      {
        dimState[iDim+1] = refData[EulerTerm::V]*(*m_nodalPrimState)[iDim+1];
      }

      // add contribution to averaged primitive variables
      RealVector& curNodeAvgVars = nodeAvgVars[iNode];
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        curNodeAvgVars[iEq] = oldAvgSolWght*curNodeAvgVars[iEq] + curInsSolWght*dimState[iEq];
      }

      // compute extra variables to be averaged
      CFuint varIdx = nbEqs;

      // rhorho
      curNodeAvgVars[varIdx] = oldAvgSolWght*curNodeAvgVars[varIdx] + curInsSolWght*dimState[0]*dimState[0]; ++varIdx;

      // 2D: uu, vv; 3D: uu, vv, ww
      for (CFuint iDim = 0; iDim < dim; ++iDim, ++varIdx)
      {
        curNodeAvgVars[varIdx] = oldAvgSolWght*curNodeAvgVars[varIdx] + curInsSolWght*dimState[iDim+1]*dimState[iDim+1];
      }
      // 2D: uv; 3D: uv, uw, vw
      for (CFuint iDim = 0; iDim < dim; ++iDim)
      {
        for (CFuint iDim2 = iDim+1; iDim2 < dim; ++iDim2, ++varIdx)
        {
          curNodeAvgVars[varIdx] = oldAvgSolWght*curNodeAvgVars[varIdx] + curInsSolWght*dimState[iDim+1]*dimState[iDim2+1];
        }
      }

      // pp
      curNodeAvgVars[varIdx] = oldAvgSolWght*curNodeAvgVars[varIdx] + curInsSolWght*dimState[nbEqsM1]*dimState[nbEqsM1];
    }

    // update m_avgStepsCounter
    ++m_avgStepsCounter;

  }

  // set m_addCurrSolToAvgSol to true (this command has been called at least once)
  m_addCurrSolToAvgSol = true;

  if (!(subSysStatus->getNbIter() % m_writeToFileRate))
  {
  CFLog(INFO,"WriteInstantAndAvgSolution::writeToFileStream --> actually writing files...");

  // open file for averaged solution
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / m_avgSolFileName;
  fpath = boost::filesystem::change_extension(fpath, std::string(".plt"));
  ostringstream add;
  if (PE::GetPE().IsParallel()) {
    add << "-P" << PE::GetPE().GetRank();
  }
  fpath = fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandleAvg =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& foutAvg = fhandleAvg->open(fpath);

  // get the nodes
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  // this is a sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
      socket_nstatesProxy.getDataHandle();
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // averaged variables in nodes
  DataHandle< RealVector> nodeAvgVars = socket_nodeAvgVars.getDataHandle();

  // we will assume that the number of nodes is the same as
  // the number of states but the connectivity might be different
  cf_assert(nodes.size() == nodalStates.getSize());

  // get dimensionality, number of equations and reference length
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // write Tecplot Header (for instantaneous and averaged solution)
  fout    << "TITLE      =  Unstructured grid data" << "\n";
  fout    << "VARIABLES  = ";
  foutAvg << "TITLE      =  Unstructured grid data" << "\n";
  foutAvg << "VARIABLES  = ";
  for (CFuint i = 0; i < dim; ++i) {
    fout    << " \"x" << i << '\"';
    foutAvg << " \"x" << i << '\"';
  }
  for (CFuint i = 0 ;  i < nbEqs; ++i) {
    std::string n = m_primVarNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout    << " " << n;
    foutAvg << " " << n;
  }
  const CFuint nbrExtraAvgVars = m_extraAvgVarNames.size();
  for (CFuint i = 0 ;  i < nbrExtraAvgVars; ++i) {
    std::string n = m_extraAvgVarNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    foutAvg << " " << n;
  }
  fout    << "\n";
  foutAvg << "\n";

  // array to store the dimensional states
  RealVector dimState(nbEqs);
  State tempState;

  // loop over trss
  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();
  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs)
  {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];

    if((trs->hasTag("inner")) && (trs->hasTag("cell")))
    {
      // we will assume that the number of nodes is the same as
      // the number of states but the connectivity might be different
      cf_assert(nodes.size() == nodalStates.getSize());

      // get element type data
      SafePtr<vector<ElementTypeData> > elementType =
        MeshDataStack::getActive()->getElementTypeData(trs->getName());

      // loop over the element types
      // and create a zone in the tecplot file
      // for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType)
      {
        ElementTypeData& eType = (*elementType)[iType];
        const CFuint nbCellsInType  = eType.getNbElems();

        // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
        if (nbCellsInType > 0)
        {
          const CFuint nbNodesInType  = eType.getNbNodes();
          std::valarray<CFuint> nodeIDs (nbNodesInType);

          // find which nodeIDs are used in the elements of this type
          vector<CFuint> nodesInType;
          for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {
            cf_assert(nbNodesInType == trs->getNbNodesInGeo(iCell));
            for (CFuint in = 0; in < nbNodesInType; ++in) {
              nodesInType.push_back(trs->getNodeID(iCell,in));
            }
          }

          // sort the vector so we can then remove duplicated nodes
          sort(nodesInType.begin(), nodesInType.end(), std::less<CFuint>());

          // remove duplicated nodes
          vector<CFuint>::iterator lastNode = unique(nodesInType.begin(),
                                                     nodesInType.end());
          nodesInType.erase(lastNode,nodesInType.end());

          // create a map from LocalIDs (in CPU) to IDs per ElementType
          typedef CFuint LocalID;
          typedef CFuint IDinType;
          CFMap<LocalID,IDinType> localToTypeID;
          localToTypeID.reserve(nodesInType.size());

          for (CFuint i = 0; i < nodesInType.size(); ++i) {
            // in the following, + 1 is due Tecplot numbering
            localToTypeID.insert(nodesInType[i],i + 1);
          }
          localToTypeID.sortKeys();

          // print zone header (for current instanteneous and averaged solution)
          // one zone per element type
          fout << "ZONE "
               << "  T=\"P" << PE::GetPE().GetRank()<< " ZONE" << iType << " " << eType.getShape() <<"\""
               << ", N=" << nodesInType.size()
               << ", E=" << nbCellsInType
               << ", F=FEPOINT"
               << ", ET=" << MapGeoEnt::identifyGeoEntTecplot(eType.getNbNodes(),eType.getGeoOrder(),dim)
               << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
               << ", AUXDATA TRS=\"" << trs->getName() << "\""
               << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
               << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
               << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
               << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
              << "\n";

          foutAvg << "ZONE "
                  << "  T=\"P" << PE::GetPE().GetRank()<< " ZONE" << iType << " " << eType.getShape() <<"\""
                  << ", N=" << nodesInType.size()
                  << ", E=" << nbCellsInType
                  << ", F=FEPOINT"
                  << ", ET=" << MapGeoEnt::identifyGeoEntTecplot(eType.getNbNodes(),eType.getGeoOrder(),dim)
                  << ", AUXDATA CPU=\"" << PE::GetPE().GetRank() << "\""
                  << ", AUXDATA TRS=\"" << trs->getName() << "\""
                  << ", AUXDATA Filename=\"" << getMethodData().getFilename().leaf() << "\""
                  << ", AUXDATA ElementType=\"" << eType.getShape() << "\""
                  << ", AUXDATA Iter=\"" << subSysStatus->getNbIter() << "\""
                  << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim() << "\""
                  << "\n";

          // print nodal coordinates and stored node variables (for current instanteneous and averaged solution)
          for (vector<CFuint>::iterator itr = nodesInType.begin(); itr != nodesInType.end(); ++itr) {

            // current node
            const CFuint nodeID = *itr;
            const Node& currNode = *nodes[nodeID];
            // node has to be printed with the right length
            for (CFuint iDim = 0; iDim < dim; ++iDim) {
              fout    << setw(20) << fixed << setprecision(12) << currNode[iDim]*refL << " ";
              foutAvg << setw(20) << fixed << setprecision(12) << currNode[iDim]*refL << " ";
            }

            // copy current state to tempState
            const RealVector& currState = *nodalStates.getState(nodeID);
            for (CFuint ieq = 0; ieq < nbEqs; ++ieq)
            {
              tempState[ieq] = currState[ieq];
            }

            // convert to primitive variables
            *m_nodalPrimState = *m_updateToPrimVar->transform(&tempState);

            // make dimensional
            Common::SafePtr< EulerVarSet > eulerVarSet = getMethodData().getUpdateVarSet().d_castTo<EulerVarSet>();
            const RealVector& refData = eulerVarSet->getModel()->getReferencePhysicalData();
            dimState[0      ] = refData[EulerTerm::RHO]*(*m_nodalPrimState)[0];
            dimState[nbEqsM1] = refData[EulerTerm::P]*(*m_nodalPrimState)[nbEqsM1];
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              dimState[iDim+1] = refData[EulerTerm::V]*(*m_nodalPrimState)[iDim+1];
            }

            // write instanteneous primitive variables to file
            fout << dimState << "\n";

            // write averaged primitive and extra variables to file
            RealVector& curNodeAvgVars = nodeAvgVars[nodeID];
            foutAvg << curNodeAvgVars << "\n";
          }

          // write  Connectivity
          for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {

            for(CFuint n = 0; n < nbNodesInType; ++n) {
              nodeIDs[n] = localToTypeID.find(trs->getNodeID(iCell, n));
            }

            // write their connectivity
            MapGeoEnt::writeTecplotGeoEntConn(fout,
                                              nodeIDs,
                                              eType.getGeoOrder(),
                                              dim);
            MapGeoEnt::writeTecplotGeoEntConn(foutAvg,
                                              nodeIDs,
                                              eType.getGeoOrder(),
                                              dim);
            fout    << "\n";
            foutAvg << "\n";
          }
        }
      }
    } //end if inner cells

    foutAvg.flush();
  } //end loop over trs
  fout.close();
  foutAvg.close();

  }
  else
  {
    CFLog(INFO,"WriteInstantAndAvgSolution::writeToFileStream --> NOT writing files this iteration...");
    fout.close();
    path cfgpath = getMethodData().getFilename();
    remove(cfgpath);
  }

  // write boundary surface data
  writeBoundarySurface();
}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolution::setup()
{
  CFAUTOTRACE;

  WriteSolution::setup();

  // set counter for number of averaging steps to zero
  m_avgStepsCounter = 0;

  // get number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();

  // get dimensionality and number of equations
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();
  const CFuint nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // resize m_nodalPrimState
  m_nodalPrimState = new State(RealVector(nbrEqs));

  // number of averaged variables in each node
  if (dim != 2 && dim != 3)
  {
    throw BadValueException (FromHere(),"WriteInstantAndAvgSolution::setup() --> dimensionality should be 2 or 3");
  }
  const CFuint nbrAvgVars = dim == 2 ? 9 : 13;

  // put variable names
  m_primVarNames.resize(nbrEqs);
  m_extraAvgVarNames.resize(nbrAvgVars-nbrEqs);
  if (dim == 2)
  {
    m_primVarNames[0] = "rho";
    m_primVarNames[1] = "u";
    m_primVarNames[2] = "v";
    m_primVarNames[3] = "p";

    m_extraAvgVarNames[0] = "rhorho";
    m_extraAvgVarNames[1] = "uu";
    m_extraAvgVarNames[2] = "vv";
    m_extraAvgVarNames[3] = "uv";
    m_extraAvgVarNames[4] = "pp";
  }
  else
  {
    m_primVarNames[0] = "rho";
    m_primVarNames[1] = "u";
    m_primVarNames[2] = "v";
    m_primVarNames[3] = "w";
    m_primVarNames[4] = "p";

    m_extraAvgVarNames[0] = "rhorho";
    m_extraAvgVarNames[1] = "uu";
    m_extraAvgVarNames[2] = "vv";
    m_extraAvgVarNames[3] = "ww";
    m_extraAvgVarNames[4] = "uv";
    m_extraAvgVarNames[5] = "uw";
    m_extraAvgVarNames[6] = "vw";
    m_extraAvgVarNames[7] = "pp";
  }

  // resize m_nodeAvgVals
  m_nodeAvgVals.resize(nbrAvgVars);

  // get averaged variables in nodes datahandle
  DataHandle< RealVector> nodeAvgVars = socket_nodeAvgVars.getDataHandle();

  // resize nodeAvgVars
  nodeAvgVars.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodeAvgVars[iNode].resize(nbrAvgVars);
  }

  // setup the update to primitive variables transformer
  m_updateToPrimVar->setup(1);

  // remove the extension and assign to filename
  m_avgSolFileName = boost::filesystem::basename(boost::filesystem::path(m_avgSolFileNameStr));
}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolution::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_nodalPrimState);

  DataHandle< RealVector> nodeAvgVars = socket_nodeAvgVars.getDataHandle();
  for (CFuint iNode = 0; iNode < nodeAvgVars.size(); ++iNode)
  {
    nodeAvgVars[iNode].resize(0);
  }
  nodeAvgVars.resize(0);

  WriteSolution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolution::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  TecWriterCom::configure(args);

  // get the physical model that we are dealing with
  // to pass it to the variable transformer
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  // create the transformer from update to primitive variables
  std::string provider =
      VarSetTransformer::getProviderName(physModel->getConvectiveName(), m_updateVarStr, "Prim");
  m_updateToPrimVar =
      Environment::Factory<VarSetTransformer>::getInstance()
      .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_updateToPrimVar.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteInstantAndAvgSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WriteSolution::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
    WriteInstantAndAvgSolution::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nodeAvgVars);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
