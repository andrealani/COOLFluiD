#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"

#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "RadiativeTransfer/RadiativeTransferOutputFVMCC.hh"
#include "FiniteVolume/CellCenterFVM.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////
      
MethodCommandProvider<RadiativeTransferOutputFVMCC, 
		      DataProcessingData, 
		      RadiativeTransferModule>
radiativeTransferFVMCCProvider("RadiativeTransferOutputFVMCC");

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >("OutputFileWall","Name of output file to write the wall values.");
  options.addConfigOption< bool >("AppendTime","Append time to file name.");
  options.addConfigOption< bool >("AppendIter","Append Iteration# to file name."); 
  options.addConfigOption< CFuint >("TID","Temperature ID within the state vector."); 
}
      
//////////////////////////////////////////////////////////////////////////////

RadiativeTransferOutputFVMCC::RadiativeTransferOutputFVMCC(const std::string& name) :
  DataProcessingCom(name),
  m_faceTrsGeoBuilder(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_nstates("nstates"),
  socket_faceAreas("faceAreas"), 
  socket_qradFluxWall("qradFluxWall"),
  m_fvmccData(CFNULL),
  m_mapTrsFaceToID(),
  m_iFace(0),
  m_currFace(CFNULL),
  m_midFaceNode(),
  m_valuesMat(),
  m_varNames()
{
  addConfigOptionsTo(this);
  
  m_nameOutputFileWall = "wall.plt";
  setParameter("OutputFileWall",&m_nameOutputFileWall);
  
  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);
  
  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
  
  m_TID = 0;
  setParameter("TID",&m_TID);
}

//////////////////////////////////////////////////////////////////////////////

RadiativeTransferOutputFVMCC::~RadiativeTransferOutputFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
RadiativeTransferOutputFVMCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_nstates);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_qradFluxWall);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::setup()
{
  DataProcessingCom::setup();
  
  m_faceTrsGeoBuilder.setup();
  
  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  // SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
  // cf_assert(fvmcc.isNotNull());
  // m_fvmccData = fvmcc->getData();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_midFaceNode.resize(dim);
  
  m_varNames.clear();
  for (CFuint i = 0; i < dim; ++i) {
    const std::string xdim = "x" + Common::StringOps::to_str(i);
    this->m_varNames.push_back(xdim);
  }
  m_varNames.push_back("T");        // temperature
  m_varNames.push_back("heatFRad"); // radiative heat flux
  
  const CFuint nbVariables = m_varNames.size() - dim;
  if (nbVariables < 1) {cf_always_assert(nbVariables > 0);}
  
  if (m_mapTrsFaceToID.size() == 0) {
    const vector< SafePtr<TopologicalRegionSet> >& trsList = getTrsList();
    
    // empty loop just to count faces
    CFuint totalNbFaces = 0;
    for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
      totalNbFaces += trsList[iTRS]->getLocalNbGeoEnts();
    }
    
    if (totalNbFaces > 0) {
      m_mapTrsFaceToID.reserve(totalNbFaces);
      m_valuesMat.resize(nbVariables, totalNbFaces);
      
      CFuint index = 0;  
      for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) {
	SafePtr<TopologicalRegionSet> trs = trsList[iTRS];
	const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  m_mapTrsFaceToID.insert(trs->getLocalGeoID(iFace), index++);
	}
      }
      cf_assert(index == totalNbFaces);
      
      m_mapTrsFaceToID.sortKeys();
    }
    
    DataHandle<CFreal> qradFluxWall = socket_qradFluxWall.getDataHandle();
    if (qradFluxWall.size() == 0) {
      qradFluxWall.resize(totalNbFaces);
    }
    cf_assert(qradFluxWall.size() == totalNbFaces);
  }
}

//////////////////////////////////////////////////////////////////////////////
     
void RadiativeTransferOutputFVMCC::unsetup()
{
  DataProcessingCom::unsetup();
  m_mapTrsFaceToID.clear();
}
      
//////////////////////////////////////////////////////////////////////////////
      
void RadiativeTransferOutputFVMCC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::executeOnTrs()
{
  CFAUTOTRACE;
  
  // AL: this needs to be reconsidered
  const std::string nsp = this->getMethodData().getNamespace();
  if (PE::GetPE().GetRank(nsp) == 0) {
    CFLog(VERBOSE, "RadiativeTransferOutputFVMCC::executeOnTrs() START\n");
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
  CFLog(VERBOSE, "RadiativeTransferOutputFVMCC::executeOnTrs() => TRS "
	<< currTrs->getName() << "\n");
  
  TrsGeoWithNodesBuilder::GeoData& geoData = m_faceTrsGeoBuilder.getDataGE();
  geoData.trs = currTrs;
  const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
  
  prepareOutputFileWall();
  
  for (m_iFace = 0; m_iFace < nbTrsFaces; ++m_iFace) {
    CFLogDebugMed("iFace = " << m_iFace << "\n");
    
    // build the GeometricEntity
    geoData.idx = m_iFace;
    
    m_currFace = m_faceTrsGeoBuilder.buildGE();
    computeWall();
    
    // release the geometric entity
    m_faceTrsGeoBuilder.releaseGE();
  }
  
  // all data are written on file at once to ease the parallel writing
  updateOutputFileWall();
  
  if (PE::GetPE().GetRank(nsp) == 0) {
    CFLog(VERBOSE, "RadiativeTransferOutputFVMCC::executeOnTrs() END\n");
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::computeWall()
{
  CFAUTOTRACE;
  
  cf_assert(m_mapTrsFaceToID.size() > 0);
  const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID());
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  
  vector<Node*>& faceNodes = *m_currFace->getNodes();
  const CFuint nbFaceNodes = faceNodes.size();
  cf_assert(m_TID < nstates[0].size());
  
  // average temperature on face
  CFreal temp = 0.;
  for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) {
    const CFuint nodeID = faceNodes[iNode]->getLocalID();
    cf_assert(nodeID < nstates.size());
    temp += nstates[nodeID][m_TID];
  }
  temp /= nbFaceNodes;
  
  updateValuesMat(0, index, temp);
  cf_assert(index < socket_qradFluxWall.getDataHandle().size());
  const CFreal heatFluxRad = socket_qradFluxWall.getDataHandle()[index];
  updateValuesMat(1, index, -heatFluxRad);
}
      
//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::updateOutputFileWall()
{
  const std::string nsp = this->getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // all processors will write their own data one after the other 
  for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(nsp); ++i) {
    if (i == PE::GetPE().GetRank(nsp)) {
      
      SafePtr<TopologicalRegionSet> currTrs = this->getCurrentTRS();
      if (currTrs->getLocalNbGeoEnts() > 0) {
	const CFuint dim = PhysicalModelStack::getActive()->getDim();
	boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
	  boost::filesystem::path(this->m_nameOutputFileWall + currTrs->getName());
	file = Framework::PathAppender::getInstance().appendAllInfo  
	  (file,this->m_appendIter,this->m_appendTime,false);   
	
	SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
	  Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
	// append to the existing file 
	// carefull here 
	ofstream& fout = fhandle->open(file, ios::app);
	
	DataHandle < Framework::Node*, Framework::GLOBAL > nodes = this->socket_nodes.getDataHandle();
	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
	SafePtr<vector<CFuint> > trsNodes = currTrs->getNodesInTrs();
	
	// this is to ensure consistency for hybrid meshes
	// take as shape the one corresponding to the maximum number of face nodes in the whole TRS 
	CFuint nbFaceNodes = 0;
	for (CFuint f = 0; f < nbTrsFaces; ++f) { 
	  nbFaceNodes = std::max(currTrs->getNbNodesInGeo(f), nbFaceNodes); 
	}
	const std::string shape = (nbFaceNodes == 3) ? "FETRIANGLE" : "FEQUADRILATERAL";
	
	// print zone header,
	// one zone per element type per cpu
	// therefore the title is dependent on those parameters
	fout << "ZONE "
	     << "  T=\"P" << PE::GetPE().GetRank("Default")<< " ZONE" << 0 << " " << shape <<"\""
	     << ", N=" << trsNodes->size()
	     << ", E=" << nbTrsFaces
	     << ", DATAPACKING=BLOCK"
	     << ", ZONETYPE=" << shape
	     << ", VARLOCATION=( [" << (dim + 1) << "-" << this->m_varNames.size() << "]=CELLCENTERED )";
	fout << "\n\n";
	
	const CFuint nbTrsNodes = trsNodes->size();
	const CFuint writeStride = 6; // this could be user defined
	
	for (CFuint iDim = 0; iDim < dim; ++iDim) {
	  for (CFuint n = 0; n < nbTrsNodes; ++n) {
	    fout.setf(ios::scientific,ios::floatfield);
	    fout.precision(12);
	    fout << (*nodes[(*trsNodes)[n]])[iDim];
	    ((n+1)%writeStride == 0) ? fout << "\n" : fout << " ";
	  }
	}  
	
	for (CFuint iVar = dim; iVar < this->m_varNames.size(); ++iVar) {
	  //   fout << "#### variable " << this->m_varNames[iVar] << "\n\n";
	  const CFuint varID = iVar-dim;
	  cf_assert(varID <  this->m_valuesMat.nbRows());
	  
	  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	    fout.setf(ios::scientific,ios::floatfield);
	    fout.precision(12);
	    
	    const CFuint index = this->m_mapTrsFaceToID.find(currTrs->getLocalGeoID(iFace));
	    cf_assert(index < this->m_valuesMat.nbCols());
	    
	    fout << this->m_valuesMat(varID, index);
	    ((iFace+1)%writeStride == 0) ? fout << "\n" : fout << " ";
	  }
	}
	
	CFMap<CFuint,CFuint> mapNodesID(nbTrsNodes);
	for (CFuint i = 0; i < nbTrsNodes; ++i) {
	  mapNodesID.insert((*trsNodes)[i],i+1);
	}
	mapNodesID.sortKeys();
	
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  const CFuint nbNodesInGeo = currTrs->getNbNodesInGeo(iFace);
	  for (CFuint in = 0; in < nbNodesInGeo; ++in) {
	    fout << mapNodesID.find(currTrs->getNodeID(iFace,in)) << " ";
	  }
	  
	  if (nbNodesInGeo < nbFaceNodes) {
	    // here you can only have the case 3 instead of 4 nodes
	    // output twice the last node ID
	    fout << mapNodesID.find(currTrs->getNodeID(iFace,nbNodesInGeo-1)) << " "; 
	  }
	  
	  fout << "\n";
	  fout.flush();
	}
	
	//closing the file
	fhandle->close();
      }
    }
    
    PE::GetPE().setBarrier(nsp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void RadiativeTransferOutputFVMCC::prepareOutputFileWall()
{
  CFAUTOTRACE;
  
 const std::string nsp = getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // only the first processor writes the header of the output file 
  if (PE::GetPE().GetRank (nsp) == 0) {
    SafePtr<TopologicalRegionSet> currTrs = this->getCurrentTRS();
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(this->m_nameOutputFileWall + currTrs->getName());
    file = Framework::PathAppender::getInstance().appendAllInfo  
      (file,this->m_appendIter,this->m_appendTime,false);   
    
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    // append to the existing file 
    ofstream& fout = fhandle->open(file);
    
    fout << "TITLE = Unstructured Surface Quantities" << "\n";
    fout << "VARIABLES = "; 
    for (CFuint i = 0; i < this->m_varNames.size(); ++i) {
      fout << this->m_varNames[i] << " ";
    }
    fout << "\n";
    fout.close();
  } 
  
  PE::GetPE().setBarrier(nsp);
}

//////////////////////////////////////////////////////////////////////////////
     
    } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
