// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "Framework/Framework.hh"
#include "TecplotWriter/WriteTecplot.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::setDataSockets
(COOLFluiD::Framework::DataSocketSink<CFreal> volumesSocket,
 COOLFluiD::Framework::DataSocketSink<COOLFluiD::Framework::Node*, 
 COOLFluiD::Framework::GLOBAL> nodesSocket)
{

  socket_volumes = volumesSocket;
  socket_nodes = nodesSocket;

  m_isSocketsSet = true;
}

//////////////////////////////////////////////////////////////////////////////

WriteTecplot::WriteTecplot() : 
  socket_volumes("Null"),
  socket_nodes("Null"),
  l_variableNames(),
  l_variableVectorPointers(),
  l_variableScale(),
  m_geoBuilder(),
  m_isSocketsSet(false)
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::setup()
{
    if (l_strategy==0) //default
          l_strategy=1;
    m_geoBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::writeFile(std::string FileName)
{
  cf_assert(m_isSocketsSet);

  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  //std::cout << " writing file...\n" ;

  // get cells and geometry
  // get geometry
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbNodes = nodes.size();
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // output vector 
  std::vector<std::vector<CFreal> > outputVector;
  prepareNodeValues(l_strategy, outputVector);

  // open file
  SelfRegistPtr<Environment::FileHandlerOutput>* fhandle =
	 Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  ofstream& outputFile = (*fhandle)->open( boost::filesystem::path(FileName) );

  // write title and variables
  outputFile << "TITLE = \"" << l_title << "\"\n";
  outputFile << "VARIABLES = X Y"; 
  //X Y E ER EI EvR EvI RMSJouleHeatSourceCoupling elCondField T Jr Ji
  for (CFuint iVariables=0;iVariables< l_variableNames.size();iVariables++)
    outputFile << " " << l_variableNames[iVariables] ;
  outputFile << "\n"; 
  // write zones
  outputFile << "ZONE N=" << nbNodes
             << ", E=" << nbCells
             << ", F=FEPOINT"
             << ", ET=" << "QUADRILATERAL"
             << "\n";

//std::cout << "names.size: " << l_variableNames.size() << "  vectorPointers.size: " << l_variableVectorPointers.size() << "  nbNodes: " << nbNodes <<"\n";
  // write data
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const Node& currNode = *nodes[iNode];
    // Output to File
    outputFile.precision(12);
    outputFile << currNode[XX]
               << " "
               << currNode[YY] ;
    for (CFuint iVectors=0; iVectors < l_variableVectorPointers.size(); iVectors++ ) {
        outputFile << " " << outputVector[iVectors][iNode] ;
    }
    outputFile << "\n"; 
  }
  // writes connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();

    std::vector<Node*> nodes = *currCell->getNodes();

    // connectivity
    if (nodes.size() == 4) {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[3]->getLocalID()+1 <<"\n";
    }
    else {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[0]->getLocalID()+1 <<"\n";
    }

    m_geoBuilder.releaseGE();
  }

  // close file
  outputFile.close();
  delete fhandle;
  resetOutput();
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::writeFileStructure(std::string FileName)
{
  cf_assert(m_isSocketsSet);

  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  //std::cout << " writing file structure...\n" ;

  // get cells and geometry
  // get geometry
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbNodes = nodes.size();
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // open file
  SelfRegistPtr<Environment::FileHandlerOutput>* fhandle = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  ofstream& outputFile = (*fhandle)->open( boost::filesystem::path(FileName) );

  // write title and variables
  outputFile << "TITLE = \"" << l_title << "\"\n";
  outputFile << "VARIABLES = X Y\n"; 
  // write zones
  outputFile << "ZONE N=" << nbNodes
             << ", E=" << nbCells
             << ", F=FEPOINT"
             << ", ET=" << "QUADRILATERAL"
             << "\n";

  //
  // write data
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) {
    const Node& currNode = *nodes[iNode];
    // Output to File
    outputFile.precision(12);
    outputFile << currNode[XX]
               << " "
               << currNode[YY] ;
    outputFile << "\n"; 
  }

  //
  // writes connectivity
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();

    std::vector<Node*> nodes = *currCell->getNodes();

    // connectivity
    if (nodes.size() == 4) {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[3]->getLocalID()+1 <<"\n";
    }
    else {
      outputFile << nodes[0]->getLocalID()+1 << " " << nodes[1]->getLocalID()+1 <<
      " " << nodes[2]->getLocalID()+1 << " " << nodes[0]->getLocalID()+1 <<"\n";
    }

    m_geoBuilder.releaseGE();
  }

  // close file
  outputFile.close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::writeOutput(std::string FileName)
{
  cf_assert(m_isSocketsSet);

  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  // get cells and geometry
  // get geometry
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbNodes = nodes.size();
//  const CFuint nbCells = cells->getLocalNbGeoEnts();

  // output vector 
  std::vector<std::vector<CFreal> > outputVector;
  prepareNodeValues(l_strategy, outputVector);

  // open file
  SelfRegistPtr<Environment::FileHandlerOutput>* fhandle2 = 
	Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();

    { std::stringstream mm; 
    mm << boost::filesystem::path(FileName);
    mm >> FileName; } 

    CFuint iNode=0;

    string line;
    std::stringstream ss;
    ifstream myfile(FileName.c_str());
    if (myfile.is_open())
    {
      int counterRighe = 1;
      while (! myfile.eof() )
      {
        getline (myfile,line);
        if (counterRighe != 1) ss << "\n";  //newline
        if ((line.find("TITLE")<5)) ss << line << " " << l_title;
        else if ((line.find("VARIABLES")<5))  {
           ss << line;
           for (CFuint iVariables=0;iVariables< l_variableNames.size();iVariables++)
               ss << " " << l_variableNames[iVariables] ;
        }  
        else if ((line.find("ZONE")<5)) ss << line;
        else if (iNode < nbNodes)
        {
              ss.precision(12);
              ss << line ; 
              for (CFuint iVectors=0; iVectors < l_variableVectorPointers.size(); iVectors++ )
                  ss  << " " << outputVector[iVectors][iNode] ;
              iNode++;
        }
        else ss << line;
        counterRighe++ ;
      }
      myfile.close();
   
      ofstream& outputFile2 = (*fhandle2)->open(boost::filesystem::path(FileName));
      outputFile2 << ss.str();
      outputFile2.close(); 
    }
    else CFLog(VERBOSE, "!!! Unable to open file !!! \n  WriteTecplot::writeOutput(" << FileName << ")" << "\n"); 

    delete fhandle2;
    resetOutput();
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::prepareNodeValues(const CFuint strategy, std::vector<std::vector<CFreal> > &outputVector)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;

  //std::cout << " prepare node values...\n" ;

  // get geometry
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  StdTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  const CFuint nbNodes = nodes.size();
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  //
  //
  // prepare output vector
  outputVector.resize( l_variableNames.size() );
  for (CFuint i=0; i< l_variableNames.size(); i++) 
      outputVector[i].resize( nbNodes );

  //
  //
  // preparing output vector values:
  //    (please note that we need scaled values)
  //////// TWO POSSIBILITIES: ///////
  //
  // -1- we already have data in nodes
  // -2- we have only data in cells
  //
  // here we can use both. this first part is done for case -2-
  //
  // local vectors:
  RealVector volumeSum(0.,nbNodes);
  RealVector cellsAroundTheNode(0.,nbNodes);

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // set the builder data and build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
  
    // used to print electromagnetic field and Joule heat source in nodes
    std::vector<Node*> localNodes = *currCell->getNodes();
    for (CFuint iLocalNode = 0; iLocalNode < localNodes.size(); ++iLocalNode) {
        // getting node ID
        CFuint nodeID = localNodes[iLocalNode]->getLocalID(); //or getLocalID() ?
        // sum:
        volumeSum[nodeID] += volumes[iCell];
        cellsAroundTheNode[nodeID]++;      
        for (CFuint iVectors=0; iVectors < l_variableVectorPointers.size(); iVectors++) 
        {
          if ( l_inNodes[iVectors] != true ) {
            if (strategy == 1)
               outputVector[iVectors][nodeID] +=  l_variableVectorPointers[iVectors][iCell] * l_variableScale[iVectors] ;
            if (strategy == 2)
               outputVector[iVectors][nodeID] +=  l_variableVectorPointers[iVectors][iCell] * volumes[iCell] * l_variableScale[iVectors] ;
          }
        }
    }
    m_geoBuilder.releaseGE();
  }
  
  //
  // this second part is setting output 
  // using both case -1- and -2-
  //
  // preparing output variables for each node
  for (CFuint iNode = 0; iNode < nbNodes; ++iNode) { 
        for (CFuint iVectors=0; iVectors < l_variableVectorPointers.size(); iVectors++) {
          if ( l_inNodes[iVectors] ) {
            // IF
            // values stored in l_variableVectorPointers[iVectors][.]
            // are values in nodes and not in cells!!
            outputVector[iVectors][iNode] =  l_variableVectorPointers[iVectors][iNode] * l_variableScale[iVectors] ;
          } else {
            // IF
            // values stored in l_variableVectorPointers[iVectors][.]
            // are values in cell centers and not in nodes!!
            if (strategy == 1)
               outputVector[iVectors][iNode] *= 1/cellsAroundTheNode[iNode] ;
            if (strategy == 2)
               outputVector[iVectors][iNode] *= 1/volumeSum[iNode] ;
          }
        }
  }
}

//////////////////////////////////////////////////////////////////////////////

void WriteTecplot::ReIm_TO_ModPhase(const CFreal Re, const CFreal Im, CFreal &Modulo, CFreal &Phase)
{
  if ( Re == 0 && Im == 0 ) { Modulo = 0. ; Phase = 0. ; }
  else {
    Modulo = sqrt( Re*Re + Im*Im ) ; 

    if ( Re == 0 && Im > 0 ) Phase = 90;
    if ( Re == 0 && Im < 0 ) Phase = -90;
    if ( Re > 0 ) Phase = atan ( Im/Re ) / COOLFluiD::MathTools::MathConsts::CFrealPi() * 180;
    if ( Re < 0 ) Phase = 180 + (atan ( Im/Re ) / COOLFluiD::MathTools::MathConsts::CFrealPi() * 180) ;

    Phase = normalizeAngle( Phase );
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal WriteTecplot::normalizeAngle(CFreal DegreeAngle) 
{
    if ( DegreeAngle < -180 ) DegreeAngle += 360;
    if ( DegreeAngle > +180 ) DegreeAngle -= 360;
    if ( DegreeAngle < -180 || DegreeAngle > +180 ) normalizeAngle( DegreeAngle ) ;

    return DegreeAngle;
}

//////////////////////////////////////////////////////////////////////////////
/*
boost::filesystem::path LorentzForceSourceTerm::constructFilename(std::string fileName)
{
  const bool isParallel = PE::GetPE().IsParallel ();

  if (isParallel) {
    std::ostringstream fname;
    fname << boost::filesystem::basename(boost::filesystem::path(fileName))
          << "-" << PE::GetPE().GetRank()
          << boost::filesystem::extension(boost::filesystem::path(fileName));
  }

  return boost::filesystem::path(fileName);
}
*/
//////////////////////////////////////////////////////////////////////////////

  } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

