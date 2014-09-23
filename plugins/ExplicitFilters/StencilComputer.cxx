#include <fstream>
#include <iomanip>
#include "StencilComputer.hh"
#include "ExplicitFilters/ExplicitFilters.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/MeshData.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Common/ConnectivityTable.hh"
#include "Common/CFMultiMap.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MapGeoEnt.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Framework/TrsGeoWithNodesBuilder.hh"
#include "Framework/GeometricEntityPool.hh"
#include <boost/progress.hpp>
#include <boost/scoped_ptr.hpp>
#include "FilterException.hh"
#include "MathTools/MathChecks.hh"

// #include "FiniteVolume/ComputeDummyStates.hh"
// #include "FiniteVolume/ComputeFaceNormalsFVMCC.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< StencilComputer,
                                   FilterData, 
                                   StencilComputer,
                                   // Framework::StencilComputerStrategy< FilterData > , 
                                   ExplicitFiltersModule > 
  // StencilComputer::PROVIDER
  stencilComputerProvider("StencilComputer");

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("Precompute","Defines if all filter stencils will be precomputed in setup phase");
  
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::configure(Config::ConfigArgs& args)
{
  Framework::StencilComputerStrategy<FilterData>::configure(args);
}
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
StencilComputer::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result=
    Framework::StencilComputerStrategy<FilterData>::providesSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
StencilComputer::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result=
    Framework::StencilComputerStrategy<FilterData>::needsSockets();
      
  result.push_back(&socket_gstates);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

StencilComputer::StencilComputer(const std::string& name) :
   Framework::StencilComputerStrategy<FilterData>(name),
   socket_gstates("gstates"),
   m_computed(false)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);  
	m_precompute = false;
  m_prepared = false;
	setParameter("Precompute",&m_precompute);
  
}

//////////////////////////////////////////////////////////////////////////////

StencilComputer::~StencilComputer()
{
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::setup()
{
  CFLog(INFO, "stencil computer setup \n");
  Framework::StencilComputerStrategy<FilterData>::setup();
    
  // set the size of the stencil data handle (1 stencil for each element)
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  
  CFLog(INFO, "FilterStencilsize = " << getMethodData().getStencils()->size() << "\n");
  
  // // Compute the stencils for explicit filtering
  // if(m_precompute){    
  //   CFLog(NOTICE, "Computing the stencils for explicit filtering... \n");
  //   compute();
  // }
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::prepareComputations()
{
  // Set DataHandles
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = 
  socket_states.getDataHandle();
  
  // Set TRS
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
  Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  Common::SafePtr<Framework::GeometricEntityPool<Framework::TrsGeoWithNodesBuilder> > geoBuilder =
  getMethodData().getGeoWithNodesBuilder();
  Framework::TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;
  
  // loop over ALL the cells to detected all the edges joining each
  // cell center and the centers of its FACE neighbors
  CFuint nbStencils = getMethodData().getStencils()->size();
  for (CFuint iStencil=0; iStencil<nbStencils; iStencil++) {
    
    
    // ============ get volume ============ //
    
    CFuint countNegativeVol = 0;
    
    // build the GeometricEntity
    geoData.idx = iStencil;
    Framework::GeometricEntity *const cell = geoBuilder->buildGE();
    CFreal volume = cell->computeVolume();
    
    if (volume < 0.0) {
      countNegativeVol++;
      std::cout.precision(14);
      CFout << "Cell [" << iStencil << "] with [" << cell->nbNodes() << "] nodes has negative volume [" << volume << "]\n";
    }
    
    if ( MathTools::MathChecks::isZero(volume)) {
      std::cout.precision(14);
      CFout << "Cell [" << iStencil << "] with [" << cell->nbNodes() << "] nodes has zero volume [" << volume << "]\n";
      // print coordinates
      for (CFuint i = 0; i < cell->nbNodes(); ++i) {
        CFout << *cell->getNode(i) << ", ";
      } CFout << "\n";
    }
    
    //release the GeometricEntity
    geoBuilder->releaseGE();
    
    // =========== get stencil width ========== //
    CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
    getMethodData().getStencil(iStencil)->setCellWidth(pow(volume,1./CFreal(dim)));
    
    // set radius before calculating the neighbors
    getMethodData().getStencil(iStencil)->setRadius(getMethodData().getStencilGridRatio() * getMethodData().getStencil(iStencil)->getCellWidth());
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::computeNeighbors(const CFuint& centreStateID,
                                       const CFuint& currStateID)
{
  CFAUTOTRACE;
  
  // Set DataHandles
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = 
    socket_states.getDataHandle();
  Framework::DataHandle<Framework::State*> gstates = 
    socket_gstates.getDataHandle();
   Framework::DataHandle<Framework::Node*, Framework::GLOBAL> nodes = 
    socket_nodes.getDataHandle();
  
  Common::SafePtr<CoordinateLinker > coordinateLinker = getMethodData().getCoordinateLinker();
  
  // Set TRS
  Common::SafePtr<Framework::MapGeoToTrsAndIdx> mapGeoToTrs = 
    Framework::MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  Common::SafePtr<Framework::TopologicalRegionSet> cells = 
    Framework::MeshDataStack::getActive()->getTrs("InnerCells");
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cellFaces =
    Framework::MeshDataStack::getActive()->getConnectivity("cellFaces");

  // Add the currState to the stencil of centreStateID
  // CFLog(INFO, "\t\tAdding Neighbour " << currStateID << "\t to stencil" << centreStateID << " \n");
  getMethodData().getStencil(centreStateID)->addElement(currStateID);
  
  // get the number of faces of element "currStateID"
  const CFuint nbNeighborFaces = cellFaces->nbCols(currStateID);
  
  // Loop over all faces of currState
  for (CFuint iFace = 0; iFace < nbNeighborFaces; ++iFace) {
    const CFuint faceID = (*cellFaces)(currStateID, iFace);
    const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(faceID);

    // Get the TRS containing the face
    Common::SafePtr<Framework::TopologicalRegionSet> trs = mapGeoToTrs->getTrs(faceID);
    const bool isBFace = mapGeoToTrs->isBGeo(faceID);

    Framework::State* neighborState = CFNULL;
    CFuint neighID = trs->getStateID(faceIdx,1);
	  
	  if(isBFace) {
      neighborState = gstates[neighID];
	  } else {
	    if (neighID==currStateID) neighID = trs->getStateID(faceIdx,0);
      neighborState = states[neighID];
	  }
    
    const CFuint neighborStateID = neighborState->getLocalID();
    
    
    if (getMethodData().getStencil(centreStateID)->isNotIncluded(neighborStateID,isBFace)) { // not yet in stencil
            
      CFreal distance = MathTools::MathFunctions::getDistance(coordinateLinker->getCoordinates(centreStateID),
                                                              coordinateLinker->getCoordinates(neighborStateID,isBFace));
      
      
      if (distance<getMethodData().getStencil(centreStateID)->getRadius()) {
        
        if (isBFace) {
          // At a boundary face, the ghost state is included. The ghostcell will not 
          // be used to search recursively for more neighbors  
          getMethodData().getStencil(centreStateID)->setDistanceToBoundary(distance);
          getMethodData().getStencil(centreStateID)->addElement(neighborStateID,isBFace);
        } 
        else {
          // The neighbor state is added to the stencil and used to search
          // recursively for more neighbors
          computeNeighbors(centreStateID,neighborStateID);
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::compute() {
  CFAUTOTRACE;
  
  
  if (!m_prepared) {
    prepareComputations();
    m_prepared = true;
  }
  
  // Recompute stencils and weights in case the stencil was not large enough
  CFuint maxNbStencilRecomputations = 3;
  CFuint iRecompute=0;
  CFuint nbStencilsToRecompute = nbStencilsToCompute();
  CFuint nbStencils = getMethodData().getStencils()->size();
  while (nbStencilsToRecompute != 0 && iRecompute < maxNbStencilRecomputations) {
    CFLog(INFO, "\nComputing " << nbStencilsToRecompute << " stencils...\n");
    
    // Progress bar
    boost::progress_display* progress(CFNULL);
    if (nbStencilsToRecompute > 500) {
      progress = new boost::progress_display(nbStencilsToRecompute);
    }
    
    // Recompute stencil with a larger radius
    CFreal enlargementFactor = 1.15;
    
    for (CFuint iStencil=0; iStencil<nbStencils; iStencil++) {
      if (getMethodData().getStencil(iStencil)->mustCompute()) {      
        if (progress != CFNULL) ++(*progress);
        
        // set radius before computations
        if (iRecompute) {
          getMethodData().getStencil(iStencil)->enlargeRadiusWithFactor(enlargementFactor);
        }
        
        // compute the neighbors of this element
        getMethodData().getStencil(iStencil)->clear();
        computeNeighbors(iStencil,iStencil);
        
        //If the distance to the boundary is too small, recalculate neighbors with larger radius
        CFreal distance = getMethodData().getStencil(iStencil)->getDistanceToBoundary();
        if (distance < 0.3*getMethodData().getStencil(iStencil)->getRadius()) {
          CFreal R = getMethodData().getStencil(iStencil)->getRadius();
          getMethodData().getStencil(iStencil)->setRadius( R + 1.5*distance);
          getMethodData().getStencil(iStencil)->clear();
          computeNeighbors(iStencil,iStencil);
        }

        try {
          postProcessStencil(iStencil);
          getMethodData().getStencil(iStencil)->setCompute(false);
        }
        catch (std::string str)
        {
          getMethodData().getStencil(iStencil)->setCompute(true);
        }
      }
    }
    nbStencilsToRecompute = nbStencilsToCompute();
    iRecompute++;
    delete progress;
  }    
  if (nbStencilsToRecompute)
    CFLog(INFO, "ERROR, no good stencils can be found for this type \n");
  
}

//////////////////////////////////////////////////////////////////////////////

CFuint StencilComputer::nbStencilsToCompute() 
{
  CFuint counter = 0;
  CFuint nbStencils = getMethodData().getStencils()->size();
  for (CFuint iStencil=0; iStencil<nbStencils; iStencil++) {
    if (getMethodData().getStencil(iStencil)->mustCompute()) {
      counter++;
    }
  }
  return counter;
}

//////////////////////////////////////////////////////////////////////////////

CFuint StencilComputer::nbInspectedStencilsToCompute() 
{
  Common::SafePtr<std::vector<CFuint> > inspectedCellIDs = getMethodData().getInspectedCellIDs();
  CFuint nbInspected = inspectedCellIDs->size();
  
  CFuint counter = 0;
  for (CFuint iInspected=0; iInspected<nbInspected; iInspected++) {
    CFuint iStencil = (*inspectedCellIDs)[iInspected];
    if (getMethodData().getStencil(iStencil)->mustCompute()) {
      counter++;
    }
  }
  return counter;
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::postProcessStencil(const CFuint& centreStateID) {
  /* Do nothing in general case, This is to be overloaded in child classes */
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::computeWithLargerRadius(Common::SafePtr<std::vector<CFuint> > stencilsToCompute, CFreal enlargementFactor) {

}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::computeInspected() {
  CFAUTOTRACE;
  
  
  if (!m_prepared) {
    prepareComputations();
    m_prepared = true;
  }
    
  // Recompute stencil with a larger radius
  CFreal enlargementFactor = 1.15;
  CFuint iRecompute=0;
  CFuint maxNbStencilRecomputations = 3;
  
    
  CFuint nbStencilsToRecompute = nbInspectedStencilsToCompute();
  while (nbStencilsToRecompute != 0 && iRecompute < maxNbStencilRecomputations) {
    CFLog(INFO, "\nComputing " << nbStencilsToRecompute << " stencils...\n");
    
    Common::SafePtr<std::vector<CFuint> > inspectedCellIDs = getMethodData().getInspectedCellIDs();
    CFuint nbInspected = inspectedCellIDs->size();
    for (CFuint iInspected=0; iInspected<nbInspected; iInspected++) {
      CFuint iStencil = (*inspectedCellIDs)[iInspected];
      if (getMethodData().getStencil(iStencil)->mustCompute()) {      
      
        // set radius before computations
        if (iRecompute) {
          getMethodData().getStencil(iStencil)->enlargeRadiusWithFactor(enlargementFactor);
        }
      
        // compute the neighbors of this element
        getMethodData().getStencil(iStencil)->clear();
        computeNeighbors(iStencil,iStencil);
      
        //If the distance to the boundary is too small, recalculate neighbors with larger radius
        CFreal distance = getMethodData().getStencil(iStencil)->getDistanceToBoundary();
        if (distance < 0.3*getMethodData().getStencil(iStencil)->getRadius()) {
          CFreal R = getMethodData().getStencil(iStencil)->getRadius();
          getMethodData().getStencil(iStencil)->setRadius( R + 1.5*distance);
          getMethodData().getStencil(iStencil)->clear();
          computeNeighbors(iStencil,iStencil);
        }

        try {
          postProcessStencil(iStencil);
          getMethodData().getStencil(iStencil)->setCompute(false);
        }
        catch (std::string str)
        {
          getMethodData().getStencil(iStencil)->setCompute(true);
        }
      }
    }  
    nbStencilsToRecompute = nbInspectedStencilsToCompute();
    iRecompute++;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::computeStencil(const CFuint& iStencil) {
  CFAUTOTRACE;
  
  if (!m_prepared) {
    prepareComputations();
    m_prepared = true;
  }
    
  // Recompute stencil with a larger radius
  CFreal enlargementFactor = 1.15;
  CFuint iRecompute=0;
  CFuint maxNbStencilRecomputations = 3;
      
  while (getMethodData().getStencil(iStencil)->mustCompute() && iRecompute < maxNbStencilRecomputations) {
    CFLog(INFO, "\nComputing stencil \n");
      
    // Only enlarge the radius after first failed attempt
    if (iRecompute) {
      getMethodData().getStencil(iStencil)->enlargeRadiusWithFactor(enlargementFactor);
    }

    // compute the neighbors of this element
    getMethodData().getStencil(iStencil)->clear();
    computeNeighbors(iStencil,iStencil);

    //If the distance to the boundary is too small, recalculate neighbors with larger radius
    CFreal distance = getMethodData().getStencil(iStencil)->getDistanceToBoundary();
    if (distance < 0.3*getMethodData().getStencil(iStencil)->getRadius()) {
      CFreal R = getMethodData().getStencil(iStencil)->getRadius();
      getMethodData().getStencil(iStencil)->setRadius( R + 1.5*distance);
      getMethodData().getStencil(iStencil)->clear();
      computeNeighbors(iStencil,iStencil);
    }

    try {
      postProcessStencil(iStencil);
      getMethodData().getStencil(iStencil)->setCompute(false);
    }
    catch (std::string str)
    {
      getMethodData().getStencil(iStencil)->setCompute(true);
    }
    iRecompute++;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::outputStencil() 
{
  CFAUTOTRACE;

   std::string output_file = "filter_stencil.plt";

   // preparation of the output
  boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(output_file);
  file = Framework::PathAppender::getInstance().appendParallel( file );

  CFLog(INFO, "Writing Explicit Filter stencils to: " << file << " \n");

  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  std::ofstream& fout = fhandle->open(file);

  // write content
  writeToFileStream(fout);
   
  //closing the file
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void StencilComputer::writeToFileStream(std::ofstream& fout)
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
  fout << " \"stencils\"";
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
      
      // ========= Zone for whole mesh =========== //
      
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
            if ( true) { // if good filter
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
            
            // print zone header
            // one sone per element type
            fout << "ZONE "
            << "  T=\"P" << Common::PE::GetPE().GetRank()<< " good filters " << iType << " " << eType.getShape() <<"\""
            << ", N=" << goodCellNodes.size()  // ---> not known
            << ", E=" << nbGoodCells  // ---> not known
            << ", F=FEPOINT"
            << ", ET=" << Framework::MapGeoEnt::identifyGeoEntTecplot
            (eType.getNbNodes(),
             eType.getGeoOrder(),
             Framework::PhysicalModelStack::getActive()->getDim())
            << ", AUXDATA CPU=\"" << Common::PE::GetPE().GetRank() << "\""
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
              // lout << "zone " << zoneCounter << " tec_element " << iCell+1 << " CF_element " << currStateID << " \n";
            } // end write connectivity
            
          } // end nGoodCells > 0
          
        } // end nbCellsInType > 0
      } // end for each type
      
      
      // ========= Zone for each stencil cells =========== //
        
      Common::SafePtr<std::vector<CFuint> > inspectedCellIDs = getMethodData().getInspectedCellIDs();
      
      for(CFuint iInspected=0; iInspected<inspectedCellIDs->size(); ++iInspected) {
        CFuint iStencil = (*inspectedCellIDs)[iInspected];
      
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
              if (iCell == iStencil) { // if we have the stencil
                Common::SafePtr<FilterStencil> stencil = getMethodData().getStencil(iStencil);
                for(CFuint i=0; i<stencil->getNbElements(); ++i) {
                  CFuint centreCellID=stencil->getElement(i);   // this can be a state ID instead
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
            
              // print zone header
              // one sone per element type
              fout << "ZONE "
              << "  T=\"P" << Common::PE::GetPE().GetRank()<< " filterstencils " << iType << " " << eType.getShape() <<"\""
              << ", N=" << badCellNodes.size()  // ---> not known
              << ", E=" << nbBadCells  // ---> not known
              << ", F=FEPOINT"
              << ", ET=" << Framework::MapGeoEnt::identifyGeoEntTecplot
              (eType.getNbNodes(),
               eType.getGeoOrder(),
               Framework::PhysicalModelStack::getActive()->getDim())
              << ", AUXDATA CPU=\"" << Common::PE::GetPE().GetRank() << "\""
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
                // lout << "zone " << zoneCounter << " tec_element " << iCell+1 << " CF_element " << currStateID << " \n";
              } // end write connectivity
            
            } // end nbBadCells > 0
        
          } // end nbCellsInType > 0
        } // end for each type
      }
      
      
      
      for(CFuint iInspected=0; iInspected<inspectedCellIDs->size(); ++iInspected) {
        CFuint iStencil = (*inspectedCellIDs)[iInspected];
      
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
              if (iCell == iStencil) { // if we have the stencil
                  CFuint centreCellID=iStencil;   // this can be a state ID instead
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
            
              // print zone header
              // one sone per element type
              fout << "ZONE "
              << "  T=\"P" << Common::PE::GetPE().GetRank()<< " filterstencilcentres " << iType << " " << eType.getShape() <<"\""
              << ", N=" << badCellNodes.size()  // ---> not known
              << ", E=" << nbBadCells  // ---> not known
              << ", F=FEPOINT"
              << ", ET=" << Framework::MapGeoEnt::identifyGeoEntTecplot
              (eType.getNbNodes(),
               eType.getGeoOrder(),
               Framework::PhysicalModelStack::getActive()->getDim())
              << ", AUXDATA CPU=\"" << Common::PE::GetPE().GetRank() << "\""
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
                fout << 2.0 << " \n";        
              
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
                // lout << "zone " << zoneCounter << " tec_element " << iCell+1 << " CF_element " << currStateID << " \n";
              } // end write connectivity
            
            } // end nbBadCells > 0
        
          } // end nbCellsInType > 0
        } // end for each type
      }
      
       // end for each inspected cell
    } //end if inner cells
  } //end loop over trs

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ExplicitFilters

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
