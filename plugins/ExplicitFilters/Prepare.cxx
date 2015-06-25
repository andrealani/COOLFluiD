#include "ExplicitFilters/ExplicitFilters.hh"
#include "Prepare.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include <iomanip>
#include "Framework/MapGeoEnt.hh"




//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Prepare, 
                      FilterData, 
                      ExplicitFiltersModule> 
PrepareProvider("Prepare");

//////////////////////////////////////////////////////////////////////////////

void Prepare::execute()
{
  CFAUTOTRACE;

  if(m_processRate == 0) return;
    
  CFuint nbIter = Framework::SubSystemStatusStack::getActive()->getNbIter();
  if((!(nbIter % m_processRate)) ||
    (nbIter == 1)) {


    CFLog(NOTICE,"-------------------------------------------------------------\n");
    Common::SafePtr<StencilComputer> stencilComputer = 
      getMethodData().getStencilComputer();
    Common::SafePtr<FilterStrategy> filter =
      getMethodData().getFilterStrategy();
    
    CFLog(NOTICE, "Computing Explicit Filtering stencils (first pass)");
    stencilComputer->compute();
    
    CFLog(NOTICE, "\nComputing Explicit Filtering weights (first pass)");
    filter->calculateAllWeights();
    
    if (m_outputBadFilters)    outputBadFilters();
    CFLog(NOTICE,"-------------------------------------------------------------\n");
    
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
Prepare::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result=
    FilterCom::needsSockets();
      
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Prepare::outputBadFilters() 
{
  CFAUTOTRACE;

  std::string output_file = "bad_filters.plt";
  std::string log_file = "bad_filters.log";


   // preparation of the output
  boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(output_file);
  file = Framework::PathAppender::getInstance().appendParallel( file );
  
  // preparation of the log
  boost::filesystem::path logfile = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(log_file);
  logfile = Framework::PathAppender::getInstance().appendParallel( logfile );

  CFLog(INFO, "Writing bad Filter stencils to: " << file << " \n");

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(file);
  
  Common::SelfRegistPtr<Environment::FileHandlerOutput> lhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& lout = lhandle->open(logfile);

  // write content
  writeBadFiltersToFileStream(fout,lout);
   
  //closing the file
  fhandle->close();
  lhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Prepare::writeBadFiltersToFileStream(std::ofstream& fout, std::ofstream& lout)
{
  CFAUTOTRACE;

  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  const CFreal refL = Framework::PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  
  CFuint zoneCounter=0;
  
  //  Tecplot Header
  fout << "TITLE      =  Unstructured grid bad cells" << "\n";
  fout << "VARIABLES  = ";
  for (CFuint i = 0; i < dim; ++i) {
    fout << " \"x" << i << '\"';
  }
  fout << " \"bad cells\"";
  fout << " \n";

  std::vector<Common::SafePtr<Framework::TopologicalRegionSet> > trsList =
    Framework::MeshDataStack::getActive()->getTrsList();

  for(CFuint iTrs= 0; iTrs < trsList.size(); ++iTrs)
  {
    Common::SafePtr<Framework::TopologicalRegionSet> trs = trsList[iTrs];

    if((trs->hasTag("inner")) && (trs->hasTag("cell")))
    {
      // we will assume that the number of nodes is the same as
      // the number of states but the connectivity might be different
      //   cf_assert(nodes.size() == nodalStates.getSize());

      Common::SafePtr<std::vector<Framework::ElementTypeData> > elementType =
        Framework::MeshDataStack::getActive()->getElementTypeData(trs->getName());
      
      // ========= Zone for good cells =========== //
      
      // loop over the element types
      // and create a zone in the tecplot file
      // for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType)
      {
        Framework::ElementTypeData& eType = (*elementType)[iType];
        
        const CFuint nbCellsInType  = eType.getNbElems();
        // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
        if (nbCellsInType > 0) 
        {          
          Framework::DataHandle<Framework::State*, Framework::GLOBAL > states = 
          socket_states.getDataHandle();
          Common::SafePtr<Common::ConnectivityTable<CFuint> > cellNodes =
          Framework::MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
          
          // find which nodeIDs are used in the elements of this type
          std::vector<CFuint> goodCellNodes(0);
          std::vector<CFuint> goodCells(0);
          
          
          for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {
            if ( getMethodData().getFilterFlag(iCell) == true) { // if good filter
              CFuint centreCellID=states[iCell]->getLocalID();   // this can be a state ID instead
              // get the number of faces of element "currStateID"
              goodCells.push_back(centreCellID);
              const CFuint nbNodes = cellNodes->nbCols(centreCellID);
              // Loop over all faces of currState
              for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
                const CFuint nodeID = (*cellNodes)(centreCellID, iNode);
                goodCellNodes.push_back(nodeID);
              }
            } 
          }
          CFuint nbGoodCells = goodCells.size();
          
          if (nbGoodCells > 0) {
            // sort the vector so we can then remove duplicated nodes
            sort(goodCellNodes.begin(), goodCellNodes.end(), std::less<CFuint>());
            // remove duplicated nodes
            std::vector<CFuint>::iterator lastNode = unique(goodCellNodes.begin(), goodCellNodes.end());
            
            goodCellNodes.erase(lastNode,goodCellNodes.end());
            
            
            // create a map from LocalIDs (in CPU) to IDs per ElementType
            typedef CFuint LocalID;
            typedef CFuint IDinType;
            Common::CFMap<LocalID,IDinType> localToTypeID;
            localToTypeID.reserve(goodCellNodes.size());
            
            for (CFuint i = 0; i < goodCellNodes.size(); ++i) {
              // in the following, + 1 is due Tecplot numbering
              localToTypeID.insert(goodCellNodes[i],i + 1);
            }
            localToTypeID.sortKeys();
            
	    const std::string nsp = getMethodData().getNamespace();
	    
            // print zone header
            // one sone per element type
            fout << "ZONE "
            << "  T=\"P" << Common::PE::GetPE().GetRank(nsp)<< " good filters " << iType << " " << eType.getShape() <<"\""
            << ", N=" << goodCellNodes.size()  // ---> not known
            << ", E=" << nbGoodCells  // ---> not known
            << ", F=FEPOINT"
            << ", ET=" << Framework::MapGeoEnt::identifyGeoEntTecplot
            (eType.getNbNodes(),
             eType.getGeoOrder(),
             Framework::PhysicalModelStack::getActive()->getDim())
            << ", AUXDATA CPU=\"" << Common::PE::GetPE().GetRank(nsp) << "\""
            << "\n";
            zoneCounter++;
            
            
            // print nodal coordinates and stored node variables
            
            for (std::vector<CFuint>::iterator itr = goodCellNodes.begin(); itr != goodCellNodes.end(); ++itr) 
            {
              // current node
              const CFuint nodeID = *itr;
              const Framework::Node& currNode = *nodes[nodeID];
              // node has to be printed with the right length
              for (CFuint iDim = 0; iDim < dim; ++iDim) {
                fout << std::setw(20) << std::fixed << std::setprecision(12)
                << currNode[iDim]*refL << " ";
              }
              fout << 0 << " \n";        
              
              //datahandle_output->printStateData(fout,stateID);
            }
            
            const CFuint nbNodesInType  = eType.getNbNodes();
            std::valarray<CFuint> nodeIDs (nbNodesInType);
            
            for(CFuint iCell=0; iCell<nbGoodCells; ++iCell) 
            {
              // collect nodes from this cell
              CFuint currStateID = goodCells[iCell];
              // get the number of faces of element "currStateID"
              const CFuint nbNodes = cellNodes->nbCols(currStateID);
              // Loop over all faces of currState
              for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
                const CFuint nodeID = (*cellNodes)(currStateID, iNode);
                nodeIDs[iNode] = localToTypeID.find(nodeID);
              }
              
              // write their connectivity
              Framework::MapGeoEnt::writeTecplotGeoEntConn(fout,
                                                           nodeIDs,
                                                           eType.getGeoOrder(),
                                                           Framework::PhysicalModelStack::getActive()->getDim());
              fout << " \n";
              lout << "zone " << zoneCounter << " tec_element " << iCell+1 << " CF_element " << currStateID << " \n";
            } // end write connectivity
            
          } // end nGoodCells > 0
          
        } // end nbCellsInType > 0
      } // end for each type
      
      
      // ========= Zone for bad cells =========== //
          
      // loop over the element types
      // and create a zone in the tecplot file
      // for each element type
      for (CFuint iType = 0; iType < elementType->size(); ++iType)
      {
        Framework::ElementTypeData& eType = (*elementType)[iType];
        
        const CFuint nbCellsInType  = eType.getNbElems();
        // Tecplot doesn't read zones with 0 elements !!! (this was a problem in parallel)
        if (nbCellsInType > 0) 
        {
           
          Framework::DataHandle<Framework::State*, Framework::GLOBAL > states = 
            socket_states.getDataHandle();
          Common::SafePtr<Common::ConnectivityTable<CFuint> > cellNodes =
            Framework::MeshDataStack::getActive()->getConnectivity("cellNodes_InnerCells");
          
          // find which nodeIDs are used in the elements of this type
          std::vector<CFuint> badCellNodes(0);
          std::vector<CFuint> badCells(0);
          
          
          for (CFuint iCell = eType.getStartIdx(); iCell < eType.getEndIdx(); ++iCell) {
            if ( getMethodData().getFilterFlag(iCell) == false) { // if bad filter
              CFuint centreCellID=states[iCell]->getLocalID();   // this can be a state ID instead
              // get the number of faces of element "currStateID"
              badCells.push_back(centreCellID);
              const CFuint nbNodes = cellNodes->nbCols(centreCellID);
              // Loop over all faces of currState
              for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
                const CFuint nodeID = (*cellNodes)(centreCellID, iNode);
                badCellNodes.push_back(nodeID);
              }
            } 
          }
          CFuint nbBadCells = badCells.size();
          
          if (nbBadCells > 0) {
            // sort the vector so we can then remove duplicated nodes
            sort(badCellNodes.begin(), badCellNodes.end(), std::less<CFuint>());
            // remove duplicated nodes
            std::vector<CFuint>::iterator lastNode = unique(badCellNodes.begin(), badCellNodes.end());
            
            badCellNodes.erase(lastNode,badCellNodes.end());
            
            
            // create a map from LocalIDs (in CPU) to IDs per ElementType
            typedef CFuint LocalID;
            typedef CFuint IDinType;
            Common::CFMap<LocalID,IDinType> localToTypeID;
            localToTypeID.reserve(badCellNodes.size());
            
            for (CFuint i = 0; i < badCellNodes.size(); ++i) {
              // in the following, + 1 is due Tecplot numbering
              localToTypeID.insert(badCellNodes[i],i + 1);
            }
            localToTypeID.sortKeys();
	    
	    const std::string nsp = getMethodData().getNamespace();
            
	    // print zone header
            // one sone per element type
            fout << "ZONE "
            << "  T=\"P" << Common::PE::GetPE().GetRank(nsp)<< " bad filters " << iType << " " << eType.getShape() <<"\""
            << ", N=" << badCellNodes.size()  // ---> not known
            << ", E=" << nbBadCells  // ---> not known
            << ", F=FEPOINT"
            << ", ET=" << Framework::MapGeoEnt::identifyGeoEntTecplot
            (eType.getNbNodes(),
             eType.getGeoOrder(),
             Framework::PhysicalModelStack::getActive()->getDim())
            << ", AUXDATA CPU=\"" << Common::PE::GetPE().GetRank(nsp) << "\""
            << "\n";
            zoneCounter++;
            
            
            // print nodal coordinates and stored node variables
            
            for (std::vector<CFuint>::iterator itr = badCellNodes.begin(); itr != badCellNodes.end(); ++itr) 
            {
              // current node
              const CFuint nodeID = *itr;
              const Framework::Node& currNode = *nodes[nodeID];
              // node has to be printed with the right length
              for (CFuint iDim = 0; iDim < dim; ++iDim) {
                fout << std::setw(20) << std::fixed << std::setprecision(12)
                << currNode[iDim]*refL << " ";
              }
              fout << 1 << " \n";        
              
              //datahandle_output->printStateData(fout,stateID);
            }
            
            const CFuint nbNodesInType  = eType.getNbNodes();
            std::valarray<CFuint> nodeIDs (nbNodesInType);
            
            for(CFuint iCell=0; iCell<nbBadCells; ++iCell) 
            {
              // collect nodes from this cell
              CFuint currStateID = badCells[iCell];
              // get the number of faces of element "currStateID"
              const CFuint nbNodes = cellNodes->nbCols(currStateID);
              // Loop over all faces of currState
              for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
                const CFuint nodeID = (*cellNodes)(currStateID, iNode);
                nodeIDs[iNode] = localToTypeID.find(nodeID);
              }
              
              // write their connectivity
              Framework::MapGeoEnt::writeTecplotGeoEntConn(fout,
                                                           nodeIDs,
                                                           eType.getGeoOrder(),
                                                           Framework::PhysicalModelStack::getActive()->getDim());
              fout << " \n";
              lout << "zone " << zoneCounter << " tec_element " << iCell+1 << " CF_element " << currStateID << " \n";
            } // end write connectivity
            
          } // end nbBadCells > 0
        
        } // end nbCellsInType > 0
      } // end for each type
    } //end if inner cells
  } //end loop over trs
  

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD
