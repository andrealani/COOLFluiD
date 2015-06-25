#include <fstream>
#include <iostream>

#include "Common/PE.hh"
#include "Common/BadValueException.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include "FiniteVolume/CellCenterFVM.hh"

#include "FiniteVolumeArcJet/FiniteVolumeArcJet.hh"
#include "FiniteVolumeArcJet/Radiation.hh"


/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeArcJet {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Radiation, DataProcessingData, FiniteVolumeArcJetModule>
RadiationProvider("Radiation");

//////////////////////////////////////////////////////////////////////////////

Radiation::Radiation(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_CellID("CellID"),
  socket_divq("divq"),
  socket_qx("qx"),
  socket_qy("qy"),
  socket_qz("qz"),
  socket_TempProfile("TempProfile"),
  //socket_radSource("radSource"),
  m_normal(),
  m_weight(),
  m_fieldSource(),
  m_fieldAbsor(),
  m_fieldAbSrcV(),
  m_fieldAbV(),
  m_In(),
  m_II(),
  m_opacities(),
  m_radSource(),
  m_Ttable(),
  m_Ptable(),
  m_sdone(),
  m_cdone(),
  m_dirs(),
  m_advanceOrder(),
  m_q(),
  m_divq(),
  m_qrAv(),
  m_divqAv()
{
  addConfigOptionsTo(this);
  
  m_nDirs = 8;
  this->setParameter("nDirs",&m_nDirs);
  
  m_dirName = Environment::DirPaths::getInstance().getWorkingDir();
  setParameter("DirName", &m_dirName);
  
  m_binTabName = "air-100.dat";
  setParameter("BinTabName", &m_binTabName); 
  
  m_writeToFile = true;
  setParameter("WriteToFile", &m_writeToFile);
  
  m_outTabName = "air-100.out";
  setParameter("OutTabName", &m_outTabName);
  
  m_useExponentialMethod = true;
  setParameter("UseExponentialMethod", &m_useExponentialMethod);
  
  m_radialData = false;
  setParameter("RadialData", &m_radialData);
  
  m_Nr = 100;
  setParameter("NbRadialPoints", &m_Nr);
}

//////////////////////////////////////////////////////////////////////////////

Radiation::~Radiation()
{
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >
    ("nDirs","Number of directions. Only allowed 8, 24, 48 and 80");
  options.addConfigOption< boost::filesystem::path >
    ("DirName","Name of the directories where the .dat file is located.");
  options.addConfigOption< string >
    ("BinTabName","Name of the .dat file");
  options.addConfigOption< bool >
    ("WriteToFile","Writing the table in ASCII");   
  options.addConfigOption< string >
    ("OutTabName","Name of the output file");    
  options.addConfigOption< bool >
    ("UseExponentialMethod","Exponential method for radiation. Explained in ICCFD7-1003");
  options.addConfigOption< bool >
    ("RadialData","radial q and divQ for the sphere case");
  options.addConfigOption< CFuint >
    ("NbRadialPoints","Number of Radial points in the sphere case");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Radiation::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_CellID);
  result.push_back(&socket_divq);
  result.push_back(&socket_qx);
  result.push_back(&socket_qy);
  result.push_back(&socket_qz);
  result.push_back(&socket_TempProfile);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Radiation::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::setup()
{
  CFAUTOTRACE;
  
  // Setting up the file containing the binary table with opacities
  
  std::cout <<"Radiation::setup() => Before setting the directory\n";
  std::cout <<"Radiation::setup() => m_dirName      = "<< m_dirName <<"\n"; 
  std::cout <<"Radiation::setup() => m_binTabName = "<< m_binTabName <<"\n";
  
  m_inFileHandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  m_binTableFile = m_dirName / boost::filesystem::path(m_binTabName);
  
  std::cout <<"Radiation::setup() => m_binTableFile = "<< m_binTableFile <<"\n";
  
  std::cout <<"Radiation::setup() => After setting the directory\n";
  
  // Check to force an allowed number of directions
  if ((m_nDirs != 8) && (m_nDirs != 24) && (m_nDirs != 48) &&
    (m_nDirs != 80)) {
    std::cout << "Radiation::setup() => This ndirs is not allowed. 8 directions is chosen \n";
    m_nDirs = 8;
  } 
  
  // Reading the table
  readOpacities();
  
  //Comments for debugging
//   CFreal p  = 1;
//   CFuint ib = 0;
  
//   CFuint stop = 1;
//   
//   for (CFuint i = 0; i < 152; ++i){
//     tableInterpolate(T[i], p, ib, val1, val2);
//     std::cout << i <<"\t\t"<< T[i] <<"\t\t"<< p <<"\t\t"<< val1 <<"\t\t"<< val2 <<"\n";
//     //if (stop == 1){break;}
//   }
//   
//   cf_assert(stop == 0);
  
  const CFuint DIM = 3;
  
  // Selecting the # of direction types depending on the option
  switch(m_nDirs) {
  case 8:     
    m_nDirTypes = 1;
    break;
  case 24:
    m_nDirTypes = 1;    
    break;
  case 48:
    m_nDirTypes = 2;
    break;
  case 80:
    m_nDirTypes = 3;
    break;
  default:		//Case nDirs == 8
    m_nDirTypes = 1;    
    break;
  }
  
  //Resizing the vectors
  m_weight.resize(m_nDirs);
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();
  
  m_sdone.resize(nbCells);
  m_cdone.resize(nbCells);
  if(m_useExponentialMethod){
    m_fieldSource.resize(nbCells);
    m_fieldAbsor.resize(nbCells);
  }
  else{
    m_fieldSource.resize(nbCells);
    m_fieldAbSrcV.resize(nbCells);
    m_fieldAbV.resize(nbCells);
  }
  m_In.resize(nbCells);
  m_II.resize(nbCells);
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++){
    m_sdone[iCell] = false;
    m_cdone[iCell] = false;
  }
  
  m_dirs.resize(m_nDirs);
  m_advanceOrder.resize(m_nDirs);
  m_q.resize(DIM);
  m_divq.resize(nbCells);
  for (CFuint i = 0; i< m_nDirs; i++) {
    m_dirs[i].resize(3);
    m_advanceOrder[i].resize(nbCells);
  }
  for (CFuint dir = 0; dir < DIM; ++dir){
    m_q[dir].resize(nbCells);
  }
  
  m_geoBuilder.setup();
  m_normal.resize(DIM, 0.); 
  
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  DataHandle<CFreal> divQ = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  DataHandle<CFreal> TempProfile = socket_TempProfile.getDataHandle();
  divQ.resize(nbCells);
  CellID.resize(nbCells);
  TempProfile.resize(nbCells);
  qx.resize(nbCells);
  qy.resize(nbCells);
  qz.resize(nbCells);
  divQ   = 0.0;
  qx     = 0.0;
  qy     = 0.0;
  qz     = 0.0;
  CellID = 0.0;
  TempProfile = 0.0;
  
  //Averages for the Sphere case
  m_qrAv.resize(m_Nr);
  m_divqAv.resize(m_Nr);
  m_qrAv   = 0;
  m_divqAv = 0;
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::getDirections()
{
  CFLog(DEBUG_MIN, "Radiation::getDirections() => Before computing the directions\n");
  CFreal pi = MathTools::MathConsts::CFrealPi();
  //std::cout << "Number of Directions = " << m_nDirs << "\n";
  switch(m_nDirs) {
  case 8:
    m_weight[0] = 4.*pi/m_nDirs;
    m_dirs[0][0] = 1./std::sqrt(3.);
    m_dirs[0][1] = 1./std::sqrt(3.);
    m_dirs[0][2] = 1./std::sqrt(3.);   
    break;
  case 24:
    m_weight[0] = 4.*pi/m_nDirs;
    m_dirs[0][0] = 0.2958759;
    m_dirs[0][1] = 0.2958759;
    m_dirs[0][2] = 0.9082483;
    break;
  case 48:
    m_weight[0] = 0.1609517;
    m_dirs[0][0] = 0.1838670;
    m_dirs[0][1] = 0.1838670;
    m_dirs[0][2] = 0.9656013;
    m_weight[1] = 0.3626469;
    m_dirs[1][0] = 0.1838670;
    m_dirs[1][1] = 0.6950514;
    m_dirs[1][2] = 0.6950514;
    break;
  case 80:
    m_weight[0] = 0.1712359;
    m_dirs[0][0] = 0.1422555;
    m_dirs[0][1] = 0.1422555;
    m_dirs[0][2] = 0.9795543;
    m_weight[1] = 0.0992284;
    m_dirs[1][0] = 0.1422555;
    m_dirs[1][1] = 1./std::sqrt(3.);
    m_dirs[1][2] = 0.8040087;
    m_weight[2] = 0.4617179;
    m_dirs[2][0] = 1./std::sqrt(3.);
    m_dirs[2][1] = 1./std::sqrt(3.);
    m_dirs[2][2] = 1./std::sqrt(3.);
    break;
  default:	// nDirs = 8
    m_weight[0] = 4.*pi/m_nDirs;
    m_dirs[0][0] = 1./std::sqrt(3.);
    m_dirs[0][1] = 1./std::sqrt(3.);
    m_dirs[0][2] = 1./std::sqrt(3.);
    break;
  }
  
  CFuint d = m_nDirTypes - 1; //Note that it has been changed, because the counters start at 0
  for (CFuint dirType = 0; dirType < m_nDirTypes; dirType++){
    for (CFuint p = 0; p <= 2; p++){
      CFuint l = p;	    //Note that it's different because the counter starts at 0
      CFuint m = (p+1) % 3; //Note a % b is the remainder of the division a/b
      CFuint n = (p+2) % 3;

      if (p == 0 || m_dirs[dirType][0] != m_dirs[dirType][1] ||
	    m_dirs[dirType][1] != m_dirs[dirType][2] || m_dirs[dirType][2] != m_dirs[dirType][0]) {
        CFLog(DEBUG_MIN, "Case1::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      if ( p+i+j+k != 0) {
		//Note that this is different because the counters are different
		d += 1;
		m_weight[d] = m_weight[dirType];
		m_dirs[d][0] = std::pow(-1.,i)*m_dirs[dirType][l];
		m_dirs[d][1] = std::pow(-1.,j)*m_dirs[dirType][m];
		m_dirs[d][2] = std::pow(-1.,k)*m_dirs[dirType][n];
		cout<<"Case1::dirTypes = " << dirType <<"\n";
		cout<< "l = " << l << " m = " << m << " n = " << n  <<"\n";
		cout<< "d = " << d <<"\n";
		cout<< "dirs[" << d <<"] = ("<<  m_dirs[d][0] <<", " << m_dirs[d][1] <<", "<<m_dirs[d][2]<<")\n";
	      }
	    }
	  }
	}
      }     
      if (m_dirs[dirType][0] != m_dirs[dirType][1] && m_dirs[dirType][1] != m_dirs[dirType][2] 
	  && m_dirs[dirType][2] != m_dirs[dirType][0]) {
	CFLog(DEBUG_MIN, "Case2::dirTypes = " << dirType <<"\n");
	CFLog(DEBUG_MIN, "l = " << l << "m = " << m << "n = " << n  <<"\n");
	for (int i = 0; i <= 1; i++) {
	  for (int j = 0; j <= 1; j++) {
	    for (int k = 0; k <= 1; k++) {
	      //Note that this is different because the counters are different
	      d += 1;
	      m_weight[d] = m_weight[dirType];
	      m_dirs[d][0] = std::pow(-1.,i)*m_dirs[dirType][l];
	      m_dirs[d][1] = std::pow(-1.,j)*m_dirs[dirType][m];
	      m_dirs[d][2] = std::pow(-1.,k)*m_dirs[dirType][n];
	      cout<<"Case2::dirTypes = " << dirType <<"\n";
	      cout<< "l = " << l << " m = " << m << " n = " << n  <<"\n";
	      cout<< "d = " << d <<"\n";
	      cout<< "dirs[" << d <<"] = ("<<  m_dirs[d][0] <<", " << m_dirs[d][1] <<", "<<m_dirs[d][2]<<")\n";
	    }
	  }
	}
      }          
    }
  }
  // Printing the Directions for debugging
  for (CFuint dir = 0; dir < m_nDirTypes; dir++) {
    CFLog(DEBUG_MIN, "Direction[" << dir <<"] = (" << m_dirs[dir][0] <<", " << m_dirs[dir][1] <<", " << m_dirs[dir][2] <<")\n");
    //cout << "Direction[" << dir <<"] = (" << m_dirs[dir][0] <<", " << m_dirs[dir][1] <<", " << m_dirs[dir][2] <<")\n";
  }
  
  CFLog(DEBUG_MIN, "Radiation::getDirections() => After computing the directions\n");
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::execute()
{
  CFAUTOTRACE;
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  DataHandle<CFreal> divQ   = socket_divq.getDataHandle();
  DataHandle<CFreal> qx = socket_qx.getDataHandle();
  DataHandle<CFreal> qy = socket_qy.getDataHandle();
  DataHandle<CFreal> qz = socket_qz.getDataHandle();
  DataHandle<CFreal> TempProfile   = socket_TempProfile.getDataHandle();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  
  // Compute the order of advance
  // Call the function to get the directions
  getDirections();
  getAdvanceOrder();
  // getAdvanceOrderMPI();
  
  m_q[0] = 0.0;
  m_q[1] = 0.0;
  m_q[2] = 0.0;  
  m_divq = 0.0;
  m_II   = 0.0;
  
  cout <<"weight = "<< m_weight <<"\n";
  
  
  for(CFuint ib = 0; ib < m_nbBins; ib++){
    //cout<<"bin = " << ib << endl;
    getFieldOppacities(ib);
    
    for(CFuint d = 0; d < m_nDirs; d++){
      //cout<<"direction = " << d << endl;
      //cout<<"m_advanceOrder["<< d<<"] = " << endl;
      //for (CFuint m = 0; m < nbCells; m++) {
        //cout<< m_advanceOrder[d][m] <<", "; 
      //}
      //cout<< endl;
      for (CFuint m = 0; m < nbCells; m++) {
	
	CFreal inDirDotnANeg = 0.;
	CFreal Ic            = 0.;
	
	CFuint iCell = std::abs(m_advanceOrder[d][m]);
	geoData.idx = iCell;    
	GeometricEntity* currCell = m_geoBuilder.buildGE();
	
	const CFuint elemID = currCell->getState(0)->getLocalID();	
	const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
	const CFuint nbFaces = faces.size();
	
	if(m_useExponentialMethod){
	  inDirDotnANeg = 0.;
	  CFreal dirDotnANeg   = 0;
	  CFreal Lc      = 0;
	  CFreal halfExp = 0;
		  
	  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	    
	    const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
	    const CFuint faceID = face->getID();
	    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
	    
	    const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
	    const CFuint DIM = 3;
	    for (CFuint dir = 0; dir < DIM; ++dir) {
	      m_normal[dir] = normals[startID+dir]*factor;
	    }
	    
	    
	    CFreal dirDotNA = m_normal[0]*m_dirs[d][0] + m_normal[1]*m_dirs[d][1] + m_normal[2]*m_dirs[d][2]; // The normal includes the area
	    
	    if(dirDotNA < 0.){
	      dirDotnANeg +=  dirDotNA;
	      State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
	      if(neighborState->isGhost() == false){
		inDirDotnANeg += m_In[neighborState->getLocalID()]*dirDotNA;
	      }
	      else {
		CFreal boundarySource = m_fieldSource[iCell];
		inDirDotnANeg += boundarySource*dirDotNA;
	      }
	    }
	  } 
	  Lc          = volumes[iCell]/(- dirDotnANeg); 
	  halfExp     = std::exp(-0.5*Lc*m_fieldAbsor[iCell]);
	  m_In[iCell] = (inDirDotnANeg/dirDotnANeg)*std::pow(halfExp,2) + (1. - std::pow(halfExp,2))*m_fieldSource[iCell];
	  Ic          = (inDirDotnANeg/dirDotnANeg)*halfExp + (1. - halfExp)*m_fieldSource[iCell];
	}
	else{
	  CFreal dirDotnAPos   = 0;
	  
	  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	    
	    
	    const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
	    const CFuint faceID = face->getID();
	    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
	    
	    const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
	    const CFuint DIM = 3;
	    for (CFuint dir = 0; dir < DIM; ++dir) {
	      m_normal[dir] = normals[startID+dir]*factor;
	    }
		
	    CFreal dirDotNA = m_normal[0]*m_dirs[d][0] + m_normal[1]*m_dirs[d][1] + m_normal[2]*m_dirs[d][2]; // The normal includes the area
	    
	    if(dirDotNA >= 0.){
	      dirDotnAPos += dirDotNA;
	    }
	    else {
	      State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
	      if(neighborState->isGhost() == false){
		inDirDotnANeg += m_In[neighborState->getLocalID()]*dirDotNA;
	      }
	      else {
		CFreal boundarySource = m_fieldSource[iCell];
		inDirDotnANeg += boundarySource*dirDotNA;
	      }
	    }
	  } 
	  m_In[iCell] = (m_fieldAbSrcV[iCell] - inDirDotnANeg)/(m_fieldAbV[iCell] + dirDotnAPos);
	  Ic = m_In[iCell];
	}
	
	m_q[0][iCell] += Ic*m_dirs[d][0]*m_weight[d];
	m_q[1][iCell] += Ic*m_dirs[d][1]*m_weight[d];
	m_q[2][iCell] += Ic*m_dirs[d][2]*m_weight[d];
	
	CFreal inDirDotnA = inDirDotnANeg;
	for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	  
	  const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
	  const CFuint faceID = face->getID();
	  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
	  
	  const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
	  const CFuint DIM = 3;
	  for (CFuint dir = 0; dir < DIM; ++dir) {
	    m_normal[dir] = normals[startID+dir]*factor;
	  }
	      
	  CFreal dirDotNA = m_normal[0]*m_dirs[d][0] + m_normal[1]*m_dirs[d][1] + m_normal[2]*m_dirs[d][2]; // The normal includes the area
	  
	  if(dirDotNA > 0.){
	    inDirDotnA += m_In[iCell]*dirDotNA;
	  }
	}

	m_divq[iCell] += inDirDotnA*m_weight[d];
	m_II[iCell]   += Ic*m_weight[d];
	
	m_geoBuilder.releaseGE();
      }  
      //cout<<"m_divq[1000] = "<< m_divq[1000] <<"\n";
    }
    cout <<"ib = \t"<< ib <<"\n";

  }
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++) {
    m_divq[iCell] = m_divq[iCell]/(volumes[iCell]); //converting area from m^3 into cm^3
    divQ[iCell] = m_divq[iCell];
    qx[iCell] = m_q[0][iCell];
    qy[iCell] = m_q[1][iCell];
    qz[iCell] = m_q[2][iCell];
  }
  //std::cout << "m_divq = " << m_divq <<"\n\n";
  //std::cout << "m_q[0] = " << m_q[0] <<"\n\n";
  //std::cout << "m_q[1] = " << m_q[1] <<"\n\n";
  //std::cout << "m_q[2] = " << m_q[2] <<"\n\n";

  if(m_radialData){
    writeRadialData();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::getFieldOppacities(CFuint ib)
{
  if(m_useExponentialMethod){
    m_fieldSource = 0.;
    m_fieldAbsor  = 0.;
  }
  else{
    m_fieldSource = 0.;
    m_fieldAbSrcV = 0.;
    m_fieldAbV    = 0.;
  }
  
  //cout <<"Radiation::getFieldOppacities==>entering"<<endl;
  
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  DataHandle<CFreal> TempProfile   = socket_TempProfile.getDataHandle();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint phiID = totalNbEqs-1;
  const CFuint TID = phiID-1;
  const CFuint pID = 0;
  
  for (CFuint iCell = 0; iCell < nbCells; iCell++) {
    geoData.idx = iCell;    
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const State *currState = currCell->getState(0); // please note the reference &
    //Get the field pressure and T commented because now we impose a temperature profile
    //CFreal p = (*currState)[pID];
    //CFreal T = (*currState)[TID];
    
    //Test temperature profile
    CFreal xPos = currCell->getState(0)->getCoordinates()[XX];
    CFreal yPos = currCell->getState(0)->getCoordinates()[YY];
    CFreal zPos = currCell->getState(0)->getCoordinates()[ZZ]; 
    CFreal r    = std::sqrt(xPos*xPos + yPos*yPos + zPos*zPos);
    
    CFreal p      = 1013250.; //Pressure in Pa
    CFreal Tmax   = 12000.;//12000;
    CFreal Tmin   = 1000.;//1000;
    CFreal deltaT = 0.0071;//0.0071;
    CFreal A      = std::pow(r*0.01/deltaT, 2);
    CFreal rmax   = 1.5;
    CFreal Amax   = std::pow(rmax*0.01/deltaT, 2);
    CFreal T      = Tmax - (Tmax - Tmin)*(1 - std::exp(-A))/(1 - std::exp(-Amax)); 
    TempProfile[iCell] = T;
    
    CFreal patm   = p/101325; //converting from Pa to atm
    
    CFreal val1 = 0;
    CFreal val2 = 0;
    
    tableInterpolate(T, patm, ib, val1, val2); 

    
    if(m_useExponentialMethod){
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbsor[iCell]  = 1e-30;
      }
      else {
	m_fieldSource[iCell] = val2/val1;
	m_fieldAbsor[iCell]  = val1;
      }
    }
    else{
      if (val1 <= 1e-30 || val2 <= 1e-30 ){
	m_fieldSource[iCell] = 1e-30;
	m_fieldAbV[iCell]    = 1e-30*volumes[iCell]; // Volumen converted from m^3 into cm^3
      }
      else {
	m_fieldSource[iCell] = val2/val1;
	m_fieldAbV[iCell]    = val1*volumes[iCell];
      }      
      m_fieldAbSrcV[iCell]   = m_fieldSource[iCell]*m_fieldAbV[iCell];
    }
    m_geoBuilder.releaseGE();
  }
  //cout <<"Radiation::getFieldOppacities==>exiting"<<endl;
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::getAdvanceOrder()
{
  // The order of advance calculation begins here
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint DIM = PhysicalModelStack::getActive()->getDim();
  cf_assert(DIM == DIM_3D);
  
  // Only for debugging purposes one can change the direction here
  //   m_dirs[0][0] = 1;
  //   m_dirs[0][1] = 0;
  //   m_dirs[0][2] = 0;
  
  CFLog(DEBUG_MIN, "Radiation::execute() => Computing order of advance. Before the loop over the directions\n");
  
  directions_loop:
  for (CFuint d = 0; d < m_nDirs; d++){
    CFLog(INFO, "Direction number =" << d <<"\n");
    CFuint mLast = 0;
    CFuint m = 0;
    CFuint stage = 1;
    // Initializing the sdone and cdone for each direction
    m_sdone.assign(nbCells, false);
    m_cdone.assign(nbCells, false);
    
    while (m < nbCells) { //The loop over the cells begins
      mLast = m;	  //Checking to see if it counts all the cells
      for (CFuint iCell = 0; iCell < nbCells; iCell++) {

	// std::cout<<"iCell = " << iCell <<"\n";
	if (m_sdone[iCell] == false) {

	  geoData.idx = iCell;
	  GeometricEntity* currCell = m_geoBuilder.buildGE();
	  const CFuint elemID = currCell->getState(0)->getLocalID();	
	  const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
	  const CFuint nbFaces = faces.size();
	  	  
	  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	    const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
	    const CFuint faceID = face->getID();
	    const CFuint startID = faceID*DIM;
	    
	    const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
	    for (CFuint dir = 0; dir < DIM; ++dir) {
	      m_normal[dir] = normals[startID+dir]*factor;
	    }
	        
	    const CFreal dotMult = m_normal[XX]*m_dirs[d][XX] + m_normal[YY]*m_dirs[d][YY] + m_normal[ZZ]*m_dirs[d][ZZ];
	    
	    State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
	    const CFuint neighborID = neighborState->getLocalID();
	    const bool neighborIsSdone =  m_sdone[neighborID] || neighborState->isGhost(); // AL: this could lead to a very subtle BUG
	    // first check if it is a ghost, otherwise it could fail if m_sdone for the neighborID of the ghost was already set to "true" and this is not a ghost
	    // ghosts have local IDs but they restart from 0 up to the maximum number of ghosts 
	    // const bool neighborIsSdone =  neighborState->isGhost() || m_sdone[neighborID]; 
	  
	    // cout << "is ghost = " << neighborState->isGhost() << endl;
	  
	    //If the dot product is negative and the neighbor is not done, it skips the cell and continue the loop
	    if (dotMult < 0. && neighborIsSdone == false) {goto cell_loop;}
	  }// end loop over the FACES
	  
	  // cout << "m_advanceOrder[" << d << "][" << m <<"] = " << iCell << endl;
	  
	  m_advanceOrder[d][m] = iCell;
	  CellID[iCell] = stage;
	  m += 1;
	  m_cdone[elemID] = true;
	}// end if(Cell is not done)
	
        cell_loop:
	m_geoBuilder.releaseGE();
      }// end of the loop over the CELLS
      
      CFLog(DEBUG_MIN, "m_advanceOrder["<< d <<"] = " << m_advanceOrder[d] << "\n");

      if (m == mLast) {		//Check that it wrote a cell in the current stage
	  std::cout << "No cell added to advance list in direction number = " << d <<". Problem with mesh.\n";
	  //cf_assert(m != mLast);
	  /// Test to rotate the directions
	  CFreal xAngleRotation;
	  CFreal yAngleRotation;
	  CFreal zAngleRotation;
	  std::cout << "Try rotating in x. Introduce an angle in degrees\n";
	  std::cout << "Theta = \n";
	  std::cin >> xAngleRotation;
	  std::cout << "Try rotating in y. Introduce an angle in degrees\n";
	  std::cout << "Phi = \n";
	  std::cin >> yAngleRotation;
	  std::cout << "Try rotating in z. Introduce an angle in degrees\n";
	  std::cout << "Psi = \n";
	  std::cin >> zAngleRotation;
	
	  xAngleRotation *= 3.141592654/180;
	  yAngleRotation *= 3.141592654/180;
	  zAngleRotation *= 3.141592654/180;
	
	  for(CFuint dirs = 0; dirs < m_nDirs; dirs++){
	    //Rotating over x
	    CFreal rot0 = m_dirs[dirs][0];
	    CFreal rot1 = m_dirs[dirs][1]*std::cos(xAngleRotation) - m_dirs[dirs][2]*std::sin(xAngleRotation);
	    CFreal rot2 = m_dirs[dirs][1]*std::sin(xAngleRotation) + m_dirs[dirs][2]*std::cos(xAngleRotation);
	    //Rotating over y
	    CFreal rot3 = rot0*std::cos(yAngleRotation) + rot2*std::sin(yAngleRotation);
	    CFreal rot4 = rot1;
	    CFreal rot5 = -rot0*std::sin(yAngleRotation) + rot2*std::cos(yAngleRotation);
	    //Rotating over z
	    CFreal rot6 = rot3*std::cos(zAngleRotation) - rot4*std::sin(zAngleRotation);
	    CFreal rot7 = rot3*std::sin(zAngleRotation) + rot4*std::cos(zAngleRotation);
	    CFreal rot8 = rot5;
	  
	    m_dirs[dirs][0] = rot6;
	    m_dirs[dirs][1] = rot7;
	    m_dirs[dirs][2] = rot8;
	    cout<< "dirs[" << dirs <<"] = ("<<  m_dirs[dirs][0] <<", " << m_dirs[dirs][1] <<", "<<m_dirs[dirs][2]<<")\n";
	  }
	  goto directions_loop;
	  CFLog(INFO, "No cell added to advance list. Problem with mesh\n");
	  cf_assert(m != mLast);
      }
      m_advanceOrder[d][m - 1] = - m_advanceOrder[d][m - 1];
      m_sdone = m_cdone;
      
      //      CFLog(INFO, "m  "<< m << " \n");
      CFLog(INFO, "End of the "<< stage << " stage \n");
      
      ++stage;
    }// end of the loop over the STAGES
    //Printing advanceOrder for debug porpuses
    CFLog(DEBUG_MIN, "m_advanceOrder["<< d <<"] = " << m_advanceOrder[d] <<"\n");  
  } //end for for directions
  cout<<"End of advance order calculation \n";
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::getAdvanceOrderMPI()
{
  // The order of advance calculation begins here
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint DIM = PhysicalModelStack::getActive()->getDim();
  cf_assert(DIM == DIM_3D);
    
  // Only for debugging purposes one can change the direction here
  //   m_dirs[0][0] = 1;
  //   m_dirs[0][1] = 0;
  //   m_dirs[0][2] = 0;
  
  CFLog(DEBUG_MIN, "Radiation::execute() => Computing order of advance. Before the loop over the directions\n");
  for (CFuint d = 0; d < m_nDirs; d++){
    CFLog(INFO, "Direction number =" << d <<"\n");
    CFuint mLast = 0;
    CFuint m = 0;
    CFuint stage = 1;
    // Initializing the sdone and cdone for each direction
    m_sdone.assign(nbCells, false);
    m_cdone.assign(nbCells, false);
    
    while (m < nbCells) { //The loop over the cells begins
      mLast = m;	  //Checking to see if it counts all the cells
      for (CFuint iCell = 0; iCell < nbCells; iCell++) {
	// 	std::cout<<"iCell = " << iCell <<"\n";
	// if the cell is parallel updatable means that all its faces will not be partition faces (something less to check), 
	// since there should be two surrounding layers of cells before reaching the actual partition faces
	if ((!m_sdone[iCell]) && states[iCell]->isParUpdatable()) { 
	  geoData.idx = iCell;
	  GeometricEntity* currCell = m_geoBuilder.buildGE();
	  const CFuint elemID = currCell->getState(0)->getLocalID();	
	  const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
	  const CFuint nbFaces = faces.size();
	  	  
	  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	    // 	    std::cout<<"iFace = " << iFace <<"\n";
	    const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
	    const CFuint faceID = face->getID();
	    const CFuint startID = faceID*DIM;
	    
	    const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
	    for (CFuint dir = 0; dir < DIM; ++dir) {
	      m_normal[dir] = normals[startID+dir]*factor;
	    }
	    
	    const CFreal dotMult = m_normal[XX]*m_dirs[d][XX] + m_normal[YY]*m_dirs[d][YY] + m_normal[ZZ]*m_dirs[d][ZZ];
	    
	    State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
	    const CFuint neighborID = neighborState->getLocalID();
	    const bool neighborIsSdone =  m_sdone[neighborID] || neighborState->isGhost(); // AL: this could lead to a very subtle BUG
	    // first check if it is a ghost, otherwise it could fail if m_sdone for the neighborID of the ghost was already set to "true" and this is not a ghost
	    // ghosts have local IDs but they restart from 0 up to the maximum number of ghosts 
	    // const bool neighborIsSdone =  neighborState->isGhost() || m_sdone[neighborID]; 
	    
	    //If the dot product is negative and the neighbor is not done, it skips the cell and continue the loop
	    if (dotMult < 0. && neighborIsSdone == false) {goto cell_loop;}
	  }// end loop over the FACES
	  
	  m_advanceOrder[d][m] = iCell;
	  CellID[iCell] = stage;
	  m += 1;
	  m_cdone[elemID] = true;
	}// end if(Cell is not done)
	
        cell_loop:
	m_geoBuilder.releaseGE();
      }// end of the loop over the CELLS
      
      CFLog(DEBUG_MIN, "m_advanceOrder["<< d <<"] = " << m_advanceOrder[d] << "\n");
      if (m == mLast) {		//Check that it wrote a cell in the current stage
	CFLog(INFO, "No cell added to advance list. Problem with mesh\n");
       	cf_assert(m != mLast); // overlap cells are ignored therefore this assertion could be failing  
      }

      cf_assert(m-1 < m_advanceOrder[d].size());
      m_advanceOrder[d][m - 1] = - m_advanceOrder[d][m - 1];
      m_sdone = m_cdone;
      CFLog(INFO, "End of the "<< stage << " stage \n");
      ++stage;
    }// end of the loop over the STAGES
    //Printing advanceOrder for debug porpuses
    CFLog(DEBUG_MIN, "m_advanceOrder["<< d <<"] = " << m_advanceOrder[d] <<"\n");  
  } //end for for directions
  cout<<"End of advance order calculation \n";
  
  cout << "Radiation::getAdvanceOrderMPI() => WAITING on Processor " << 
    PE::GetPE().GetRank(getMethodData().getNamespace()) << endl;
  for(;;){}
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::readOpacities()
{
  std::cout <<"Radiation::setup() => Before opening the binary File\n";
  fstream& is = m_inFileHandle->openBinary(m_binTableFile);
    
    // the three first numbers are #bins, #Temps, #pressures
    vector<double> data(3);
    is.read((char*)&data[0], 3*sizeof(double));
    
    m_nbBins  = ((int) data[0]);
    m_nbTemp  = ((int) data[1]);
    m_nbPress = ((int) data[2]);
    
    vector<double> Pressures(m_nbPress);
    is.read((char*)&Pressures[0], m_nbPress*sizeof(double));
    
    vector<double> Temperatures(m_nbTemp);
    is.read((char*)&Temperatures[0], m_nbTemp*sizeof(double));
    
    // in the table there are 143 zeros 
    vector<double> Zeros(143);
    is.read((char*)&Zeros[0], 143*sizeof(double));
    
    // Reading the table: for each pressure
    // each bin: value1(Temp1) value2(Temp1) ... value1(TempN) value2(TempN)
    vector< vector<double> > bins(m_nbPress*m_nbBins,vector<double>(2*m_nbTemp)); //initialize  
    for (CFuint ib =0; ib < m_nbPress*m_nbBins; ++ib){
      is.read((char*)&bins[ib][0], 2*m_nbTemp*sizeof(double));
    }  
    
    double end;
    is.read((char*)&end, sizeof(double));
    
    // Verifying that the last number is zero, so the table is finished
    cf_assert(int(end) == 0);
    
    is.close();
  //}
  std::cout <<"Radiation::setup() => After closing the binary File\n";
  
  // Storing the opacities and the radiative source into memory
  m_opacities.resize(m_nbPress*m_nbBins*m_nbTemp);
  m_radSource.resize(m_nbPress*m_nbBins*m_nbTemp);
    
  // Setting up the temperature and pressure tables
  m_Ttable.resize(m_nbTemp);
  m_Ptable.resize(m_nbPress);
  
  for(CFuint ib = 0; ib < m_nbBins; ++ib) {
    for(CFuint ip = 0; ip < m_nbPress; ++ip) {
      for(CFuint it = 0; it < m_nbTemp; ++it) {
	m_opacities[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] = bins[ip + ib*m_nbPress][2*it];
	m_radSource[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] = bins[ip + ib*m_nbPress][2*it + 1];
      }
    }
  }
  for(CFuint it = 0; it < m_nbTemp; ++it){
    m_Ttable[it] = Temperatures[it];
  }
  std::cout <<"\n";
  for(CFuint ip = 0; ip < m_nbPress; ++ip){
    m_Ptable[ip] = Pressures[ip];
  }
 
  // Writting the table into a file
  if(m_writeToFile){
    std::cout <<"Radiation::setup() => Writting file \n";
    boost::filesystem::path file = m_dirName / boost::filesystem::path(m_outTabName);
    file = Framework::PathAppender::getInstance().appendParallel( file );
    
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& fout = fhandle->open(file);
    
    fout <<"#Bins = "<< data[0] <<"\t#Temps = "<< data[1] <<"\t#Pressures = "<< data[2] << endl;
    
    fout << endl;
    fout <<"Pressures[atm]  = ";
    for(CFuint ip = 0; ip < m_nbPress; ++ip){
      fout << Pressures[ip] << " ";
    }
    fout << endl;
    fout << endl;
    fout <<"Temperatures[K] = ";
    for(CFuint it = 0; it < m_nbTemp; ++it){
      fout << Temperatures[it] << " ";
    }
    fout << endl;
    CFuint m = 0;
    fout.precision(18);
    for(CFuint ip = 0; ip < m_nbPress; ++ip){
      fout << endl;
      fout <<"Pressure = "<< Pressures[ip] << endl;
      fout <<"bin \t\t\t\t Temp \t\t\t\t val1 \t\t\t\t\t\t val2" << endl;
      CFuint beginLoop = m*ip;
      for(CFuint ib = 0; ib < m_nbBins; ++ib){
	for(CFuint it = 0; it < m_nbTemp; ++it){
	  fout << ib + 1 <<"\t\t\t\t" << Temperatures[it] <<"\t\t\t\t" << m_opacities[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] <<"\t\t\t\t" << m_radSource[it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp] << endl;
	}
      }
      m += m_nbBins; 
    }
    fhandle->close();     
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void Radiation::tableInterpolate(CFreal T, CFreal p, CFuint ib, CFreal& val1, CFreal& val2)
{
  
  //Find the lower bound fo the temperature and the pressure ranges
  //we assume that the temperature and pressure always fall in the bounds.
  //If they don't then the value are still interpolated from the nearest
  //two points in the temperature or pressure list
  CFuint it = m_nbTemp - 2;
  for (CFuint i = 1; i < (m_nbTemp - 2); i++){
    if(m_Ttable[i] > T) { it = i - 1; break;}
  }
  
  CFuint ip = m_nbPress - 2;
  for (CFuint i = 1; i < (m_nbPress - 2); i++){
    if(m_Ptable[i] > p) { ip = i - 1; break;}
  }
  
  //Linear interpolation for the pressure
  
  const CFuint iPiBiT           = it + ib*m_nbTemp + ip*m_nbBins*m_nbTemp;
  const CFuint iPplus1iBiT      = it + ib*m_nbTemp + (ip + 1)*m_nbBins*m_nbTemp;
  const CFuint iPiBiTplus1      = (it + 1) + ib*m_nbTemp + ip*m_nbBins*m_nbTemp;
  const CFuint iPplus1iBiTplus1 = (it + 1) + ib*m_nbTemp + (ip + 1)*m_nbBins*m_nbTemp;
  
  // Linear interpolation for the pressure
  // Interpolation of the opacities
  const CFreal bt1op = (m_opacities[iPplus1iBiT] - m_opacities[iPiBiT])*
		    (p - m_Ptable[ip])/(m_Ptable[ip + 1] - m_Ptable[ip]) + m_opacities[iPiBiT];
  
  const CFreal bt2op = (m_opacities[iPplus1iBiTplus1] - m_opacities[iPiBiTplus1])*
		    (p - m_Ptable[ip])/(m_Ptable[ip + 1] - m_Ptable[ip]) + m_opacities[iPiBiTplus1];
  
  // Interpolation of the source
  const CFreal bt1so = (m_radSource[iPplus1iBiT] - m_radSource[iPiBiT])*
		    (p - m_Ptable[ip])/(m_Ptable[ip + 1] - m_Ptable[ip]) + m_radSource[iPiBiT];
  
  const CFreal bt2so = (m_radSource[iPplus1iBiTplus1] - m_radSource[iPiBiTplus1])*
		    (p - m_Ptable[ip])/(m_Ptable[ip + 1] - m_Ptable[ip]) + m_radSource[iPiBiTplus1];    
  
  // Logarithmic interpolation for the temperature
  // Protect against log(0) and x/0 by switching to linear interpolation if either
  // bt1 or bt2 == 0.  (Note we can't allow log of negative numbers either)
  // Interpolation of the opacities   
  if(bt1op <= 0 || bt2op <= 0){
    val1 = (bt2op - bt1op)*(T - m_Ttable[it])/(m_Ttable[it + 1] - m_Ttable[it]) + bt1op;
//    cout <<"\nOption1 \n";
//    cout <<"T = "<< T <<"\tTi+1 = "<<m_Ttable[it + 1]<<"\tTi = "<<m_Ttable[it] <<"\n";
//    cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  else {
    val1 = std::exp((T - m_Ttable[it])/(m_Ttable[it + 1] - m_Ttable[it])*std::log(bt2op/bt1op))*bt1op;
//     cout <<"\nOption2 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<m_Ttable[it + 1]<<"\tTi = "<<m_Ttable[it] <<"\n";
//     cout <<"val1 = " << val1 <<"\tbt2op ="<< bt2op <<"\tbt1op ="<< bt1op <<"\n";
  }
  // Interpolation of the source
  if(bt1so <= 0 || bt2so <= 0){
    val2 = (bt2so - bt1so)*(T - m_Ttable[it])/(m_Ttable[it + 1] - m_Ttable[it]) + bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<m_Ttable[it + 1]<<"\tTi = "<<m_Ttable[it] <<"\n";
//     cout <<"val1 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  else {
    val2 = std::exp((T - m_Ttable[it])/(m_Ttable[it + 1] - m_Ttable[it])*std::log(bt2so/bt1so))*bt1so;
//     cout <<"\nOption3 \n";
//     cout <<"T = "<< T <<"\tTi+1 = "<<m_Ttable[it + 1]<<"\tTi = "<<m_Ttable[it] <<"\n";
//     cout <<"val2 = " << val2 <<"\tbt2so ="<< bt2so <<"\tbt1so ="<< bt1so <<"\n";
  }
  
  //cf_assert(ib == 0);
  
}
//////////////////////////////////////////////////////////////////////////////
void Radiation::writeRadialData()
{
  std::cout<<"Writting radial data for the spherical test case \n";

  boost::filesystem::path file = m_dirName / boost::filesystem::path("radialData.plt");
  file = Framework::PathAppender::getInstance().appendParallel( file );
  
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& outputFile = fhandle->open(file);
  
  outputFile << "TITLE  = Radiation radial data for a sphere\n";
  outputFile << "VARIABLES = r  qr divq nbPoints\n";

  DataHandle<CFreal> CellID = socket_CellID.getDataHandle();
  
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;
 
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  CFuint nbPoints;
  CFreal rCoord;
  const CFreal Radius = 1.5;
  
  for(CFuint ir = 0; ir < m_Nr; ir++){
    nbPoints = 0;
    rCoord = (ir + 0.5)*Radius/m_Nr; //middle point between ir and (ir + 1)
    
    for(CFuint iCell = 0; iCell < nbCells; iCell++){
      geoData.idx = iCell;
      GeometricEntity* currCell = m_geoBuilder.buildGE();
      
      Node& coordinate = currCell->getState(0)->getCoordinates();
      CFreal x = coordinate[XX];
      CFreal y = coordinate[YY];
      CFreal z = coordinate[ZZ];
      
      CFreal rCell = std::sqrt(x*x + y*y + z*z);
      
      if(rCell >= ir*Radius/m_Nr && rCell < (ir + 1)*Radius/m_Nr){
	nbPoints++;
	m_divqAv[ir] += m_divq[iCell];
	m_qrAv[ir]   += (m_q[0][iCell]*x + m_q[1][iCell]*y + m_q[2][iCell]*z)/rCell; //*rCell*rCell; Multiply by r**2 for area-weighted average
      }
      m_geoBuilder.releaseGE();
    }
    if(nbPoints > 0){
      m_divqAv[ir] /= nbPoints;
      m_qrAv[ir]   /= nbPoints; //m_qrAv[ir]   /= nbPoints*rCoord*rCoord; //area-weighted average radial flux
      outputFile << rCoord << " " << m_qrAv[ir] << " " << m_divqAv[ir] << " " <<  nbPoints << "\n";
    }
  }
  fhandle->close();
  std::cout<<"Finished writting the radial data\n";
  
}

//////////////////////////////////////////////////////////////////////////////

void Radiation::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolumeArcJet

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

