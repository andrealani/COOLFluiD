#include <fstream>

#include "Common/CFLog.hh"
#include "Common/DebugFunctions.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/CFEnv.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/PhysicalConsts.hh"
#include "RadiativeTransfer/RadiationLibrary/Models/ArcJet/ArcJetRadiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ArcJetRadiator,
			    Radiator,
			    RadiativeTransferModule,
			    1>
ArcJetRadiatorProvider("ArcJetRadiator");

//////////////////////////////////////////////////////////////////////////////

void ArcJetRadiator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string > ("BinTabName","Name of the dat file");
}

//////////////////////////////////////////////////////////////////////////////

ArcJetRadiator::ArcJetRadiator(const std::string& name) :
  Radiator(name),
  m_binTabName()
{

  addConfigOptionsTo(this);
  
  m_binTabName = "air-100.dat";

  setParameter("BinTabName", &m_binTabName ); 
   
}
      
//////////////////////////////////////////////////////////////////////////////

ArcJetRadiator::~ArcJetRadiator()
{
}

//////////////////////////////////////////////////////////////////////////////

void ArcJetRadiator::configure ( Config::ConfigArgs& args )
{
  Radiator::configure(args);
  ConfigObject::configure(args);
}
      
//////////////////////////////////////////////////////////////////////////////
      
void ArcJetRadiator::setup()
{
  Radiator::setup();
  
  m_inFileHandle  = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  //read opacities
  readOpacities();
  
  //get the statesID used for this Radiator
  m_radPhysicsPtr->getCellStateIDs( m_statesID );
  
  // set up the TRS states
  m_pstates.reset
    (new Framework::DofDataHandleIterator<CFreal, State, GLOBAL>(m_states, &m_statesID));
}
  
//////////////////////////////////////////////////////////////////////////////

void ArcJetRadiator::unsetup()
{

}

/////////////////////////////////////////////////////////////////////////////

void ArcJetRadiator::setupSpectra(CFreal wavMin, CFreal wavMax)
{ 

  m_wavMin = wavMin;
  m_wavMax = wavMax;
  m_dWav = 1.;

  genData();
}
      
//////////////////////////////////////////////////////////////////////////////
inline void ArcJetRadiator::getSpectralIdxs(CFreal lambda, CFuint *idx)
{
  *idx = static_cast< CFuint >( lambda );
}

/// m_data: matrix storing absorption and emission coefficients
/// wavelength       = m_radCoeff(local state ID, spectral point idx*3)
/// emission coeff   = m_radCoeff(local state ID, spectral point idx*3+1)
/// absorption coeff = m_radCoeff(local state ID, spectral point idx*3+2)
///
CFreal ArcJetRadiator::getEmission(CFreal lambda, RealVector &s_o)
{
  CFuint spectralIdx;
  CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();

  getSpectralIdxs(lambda, &spectralIdx);

  return m_data( stateIdx, spectralIdx*3+1 );

}

//////////////////////////////////////////////////////////////////////////////
CFreal ArcJetRadiator::getAbsorption(CFreal lambda, RealVector &s_o)
{
  CFuint spectralIdx;
  CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();

  getSpectralIdxs(lambda, &spectralIdx);
  
  return m_data( stateIdx, spectralIdx*3   );
}

void ArcJetRadiator::computeEmissionCPD()
{
  CFuint nbCells = m_statesID.size();
  CFuint nbCpdPoints = m_nbBins;

  m_spectralLoopPowers.resize( nbCells );
  m_cpdEms.resize( nbCells * nbCpdPoints );

  CFreal emInt;
  
  for (CFuint s = 0; s < nbCells; ++s) {
    //first pass to get the total emission coefficient
    //use this information to get the cell's total Radiative Power

    //cout<<"firrst pass, cell "<<s<<' '<<m_dWav<<endl;
    emInt=0.;
    m_cpdEms[ s*nbCpdPoints + 0 ] = 0;
    for (CFuint j=1; j < nbCpdPoints; j++ ){
      //cout<<"spectral point "<<j<<" of "<<nbCpdPoints<<endl;
      emInt +=  m_data(s, (j-1)*3+1 );
      m_cpdEms[ s*nbCpdPoints + j ] = emInt;
    }
//    cout<<"done, cell "<<s<<' '<<emInt <<endl;
    m_spectralLoopPowers[s]=emInt*getCellVolume( m_pstates->getStateLocalID(s) ) * 4.0 * 3.1415926535897;

    //cout<<"second pass, cell "<<s<<endl;
    //second pass to get the [0->1] comulative probably distribution
    for (CFuint j=1; j < nbCpdPoints; j++ ){
      m_cpdEms[ s*nbCpdPoints + j ] /= emInt;
    }
  }

//  CFuint m_nbPoints = 100;
//  if (  PE::GetPE().GetRank()  == 0) {
//  //test for the first cell
//  cout<<endl<<" wav = [";
//  for (CFuint j=0; j < m_nbPoints; j++ ){
//    cout<<m_data(0,j*3+0)<<' ';
//  }
//  cout<<" ];" <<endl<<" em = [";
//  for (CFuint j=0; j < m_nbPoints; j++ ){
//    cout<<m_data(0,j*3+1)<<' ';
//  }
//  cout<<" ];"<<endl<<" am = [";
//  for (CFuint j=0; j < m_nbPoints; j++ ){
//    cout<<m_data(0,j*3+2)<<' ';
//  }
//
//
////  cout<<" ];"<<endl<<" cmd = [ ";
////  for (CFuint j=0; j < nbCpdPoints; j++ ){
////    cout<<m_cpdEms[ j ] <<' ';
////  }
////
////  cout<<" ];"<<endl<<" integral= "<<m_spectralLoopPowers[0]/getCellVolume( m_pstates->getStateLocalID(0))<<' '<<endl;
////
////  cout<<"state idx= "<<m_radPhysicsHandlerPtr->getCurrentCellTrsIdx()<<endl;
////  cout<<"PDF = [";
////  RealVector so(3);
////  CFreal wav;
////  for(CFuint i=0;i<1000;++i){
////    for(CFuint j=0;j<10;++j){
////      getRandomEmission(wav,so);
////      cout<< wav <<' ';
////    }
////    cout<<" ..."<<endl;
////  }
////  cout<<" ];"<<endl;
//}


}

//////////////////////////////////////////////////////////////////////////////

CFreal ArcJetRadiator::getSpectraLoopPower()
{
  return m_spectralLoopPowers[ m_radPhysicsHandlerPtr->getCurrentCellTrsIdx() ];
}


/////////////////////////////////////////////////////////////////////////////

void ArcJetRadiator::getRandomEmission(CFreal &lambda, RealVector &s_o)
{
  //cout<<"get emission"<<endl;
  static CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  static CFuint dim2 = m_radPhysicsHandlerPtr->isAxi() ? 3 : dim;

  CFuint nbCpdPoints = m_nbBins;
  CFuint stateIdx = m_radPhysicsHandlerPtr->getCurrentCellTrsIdx();
  
  //cout<<"cpd, stateIdx "<<stateIdx<<" :";
  //vector<CFreal>::iterator it;
  //for(it = it_start; it<=it_end; ++it ){
    //cout<<*it<<' ';
  //}
  //cout<<endl;
  
  CFreal rand = m_rand.uniformRand();
  CFreal* it_start = &m_cpdEms[ stateIdx*nbCpdPoints ];
  CFreal* it_end = &m_cpdEms[ (stateIdx+1)*nbCpdPoints-1 ];
  CFreal* it_upp = std::upper_bound(it_start, it_end, rand);
  CFuint spectralIdx = it_upp - it_start;
  
  //cout<<stateIdx<<dx<<' '<<spectralIdx1<<' '<<spectralIdx2<<endl;

  lambda = m_data( stateIdx, spectralIdx*3 );

  //cout<<x0<<' '<<x1<<' '<<y0<<' '<<y1<<' '<<rand<<' '<<lambda<<endl;

  m_rand.sphereDirections(dim2, s_o);
}

void ArcJetRadiator::readOpacities()
{
// assuming the following format:
  // nbBins
  // nbTems
  // nbPress
  // Pressure[0] (Pa) 
  // Pressure[1] (Pa)
  // ...
  // Pressure[nbPress-1] (Pa)
  // Temperature[0] (K) 
  // Temperature[1] (K) 
  // ...
  // Temperature[nbTemps-1] (K) 
  // noting[0] 
  // noting[1]
  // ...
  // nothing[142]
  // m_opacities[1][1][1] ( [Press][Temp][Bin] ) 
  // m_radSource[1][1][1]
  // m_opacities[1][2][1]
  // m_radSource[1][2][1]
  // ...
  // m_opacities[1][nbtemps-1][1]
  // m_radsource[1][nbtemps-1][1]
  // ...
  // m_opacities[2][1][2]
  // m_radSource[2][1][2]
  // m_opacities[2][nbtemps-1][1]
  // m_radsource[2][nbtemps-1][1]
  //...
  // m_opacities[nbPress][nbTemps-1][1]
  // m_radSource[nbPress][nbTemps-1][1]
  //...
  // m_opacities[1][1][2]
  // m_radSource[1][1][2]
  // m_opacities[1][2][2]
  // m_radSource[1][2][2]
  // ...
  // ...
  // m_opacities[nbPress-1][nbTemps-1][nbBins-1]
  // m_radSource[nbPress-1][nbTemps-1][nbBins-1]
  // EOF

  boost::filesystem::path path(m_binTabName);
  m_binTableFile = Environment::DirPaths::getInstance().getWorkingDir() / path;

  fstream& is = m_inFileHandle->openBinary(m_binTableFile);
	
  // Reading nbBins, nbTemps, nbPressures
  vector<double> data(3);
  is.read((char*)&data[0], 3*sizeof(double));
  
  m_nbBins  = ((int) data[0]);
  m_nbTemp  = ((int) data[1]);
  m_nbPress = ((int) data[2]);
  

  // Alocating memory for the 3d matrixes 
  m_opacities.resize(m_nbPress);
  m_radSource.resize(m_nbPress);
  for (CFuint i = 0; i < m_nbPress; ++i) {
    m_opacities[i].resize(m_nbTemp);
    m_radSource[i].resize(m_nbTemp);
    for (CFuint j = 0; j < m_nbTemp; ++j){
      m_opacities[i][j].resize(m_nbBins);
      m_radSource[i][j].resize(m_nbBins);
    }
  } 
  
  // Reading the temperature and pressure tables
  m_Ptable.resize(m_nbPress);
  is.read((char*)&m_Ptable[0], m_nbPress*sizeof(double));
  
  m_Ttable.resize(m_nbTemp);
  is.read((char*)&m_Ttable[0], m_nbTemp*sizeof(double));
  
  // wavelengths: just the IDs 
  m_wavTable.resize(m_nbBins);
  for(CFuint i=0; i<m_nbBins ; ++i){
     m_wavTable[i] = i;
  }

  // in the table there are 143 zeros 
  vector<double> Zeros(143);
  is.read((char*)&Zeros[0], 143*sizeof(double));

  //reading the contents of the table

  //file is small, so it can be read in one go
  std::vector<CFreal> temp(m_nbTemp*2*m_nbPress*m_nbBins);
  is.read((char*)&temp[0], m_nbTemp*2*m_nbPress*m_nbBins*sizeof(double)); 
  
  CFuint step = m_nbPress * m_nbTemp*2;
  for(CFuint b = 0; b < m_nbBins; ++b){
    for(CFuint p = 0; p < m_nbPress; ++p ){
      for(CFuint t=0; t<m_nbTemp; ++t){
        m_radSource[p][t][b] = temp[ t*2   + p*m_nbPress + b*step];
        m_opacities[p][t][b] = temp[ t*2+1 + p*m_nbPress + b*step];
      }
    }
  }


//  cout<< "em00 = [ ";
//  for(CFuint b = 0; b < m_nbBins; ++b){
//        cout << m_radSource[0][19][b] <<' ';
//  }
//  cout<<" ];"<<endl;

}

//////////////////////////////////////////////////////////////////////////////
void ArcJetRadiator::tableInterpolate(CFreal T, CFreal P, CFuint ib, CFreal& val1, CFreal& val2)
{
  
  //Binary searach for the upper and lower bound fo the temperature and the 
  //pressure ranges
  //we assume that the temperature and pressure always fall in the bounds.
  
  std::vector<CFreal>::iterator it_uppT, it_uppP;
  it_uppT = std::upper_bound(m_Ttable.begin(), m_Ttable.end(), T );
  it_uppP = std::upper_bound(m_Ptable.begin(), m_Ptable.end(), P );

  it_uppP = ( it_uppP == m_Ptable.begin() )? (++it_uppP) : it_uppP;
  it_uppT = ( it_uppT == m_Ttable.begin() )? (++it_uppT) : it_uppT;

  CFuint iT2 = it_uppT - m_Ttable.begin();
  CFuint iT1 = iT2 -1 ;

  CFuint iP2 = it_uppP - m_Ptable.begin();
  CFuint iP1 = iP2 -1 ;


  // Bilinear interpolation
  CFreal T1 = m_Ttable[iT1];
  CFreal T2 = m_Ttable[iT2];
  
  CFreal P1 = m_Ptable[iP1];
  CFreal P2 = m_Ptable[iP2];

  // Opacity
  CFreal K11 = m_opacities[iP1][iT1][ib];
  CFreal K12 = m_opacities[iP1][iT2][ib];
  CFreal K21 = m_opacities[iP2][iT1][ib];
  CFreal K22 = m_opacities[iP2][iT2][ib];

  // Source Term
  CFreal E11 = m_radSource[iP1][iT1][ib];
  CFreal E12 = m_radSource[iP1][iT2][ib];
  CFreal E21 = m_radSource[iP2][iT1][ib];
  CFreal E22 = m_radSource[iP2][iT2][ib];

  // distances 
  CFreal pp1 = P - P1; //x direction: P
  CFreal p2p = P2 - P; //y direction: T
  CFreal tt1 = T - T1;
  CFreal t2t = T2 - T;
  
  // scaling factor
  CFreal frac = 1/( (T2-T1)*(P2-P1) );

  // Interpolated emission coeff / source term
  val1 = frac * (   E11 * p2p * t2t 
		  + E21 * pp1 * t2t
		  + E12 * p2p * tt1
		  + E22 * pp1 * tt1
		); 

  // Interpolated absorption coeff / opacity
  val2 = frac * (   K11 * p2p * t2t 
		  + K21 * pp1 * t2t
		  + K12 * p2p * tt1
		  + K22 * pp1 * tt1
		); 


//  cout<< T1 << ' ' <<T2 << ' '<<P1<<' '<<P2 <<endl
//  cout<< iT1 << ' ' <<iT2 << ' '<<iP1<<' '<<iP2 <<endl;
//  cout<<"val1: "<<val1<<" val2: "<<val2<<endl;
//  cout<<E11<<' '<<E12<<' '<<E21<<' '<<E22<<endl;
//  cout<<K11<<' '<<K12<<' '<<K21<<' '<<K22<<endl;
}

void ArcJetRadiator::genData(){
  
  const CFuint tempID = m_radPhysicsHandlerPtr->getTempID();
  const CFuint pressID = 0;
  const CFuint nbStates = m_statesID.size();
 
//  CFLog(INFO, "nbStates: "<<nbStates<<'\n' );

  // loop dataBin Wavelenghts to find the min and max id's

  CFuint minIdx = 0;
  CFuint maxIdx = 0;

  for (CFuint i=0; i< m_wavTable.size() ; ++i){
    minIdx = m_wavTable[i] < m_wavMin ? i   : minIdx;
    maxIdx = m_wavTable[i] < m_wavMax ? i+1 : maxIdx;
  }

  const CFuint nbPoints = maxIdx - minIdx;

  cf_assert(nbPoints > 0 && 
		  " Requested spectral loop out of range of the radiation Tables" );

  m_data.resize(nbStates, (nbPoints+2)*3*2 );
  
  CFreal kappa, ems, tempState, pressState;
  for (CFuint itState=0; itState< nbStates; ++itState){
    for (CFuint itBin=minIdx; itBin< maxIdx; ++itBin ){
      //CFLog(INFO, "itState: "<<itState<<'\n');
      CFreal *const currState = m_pstates->getState( itState );
      tempState = currState[tempID];
      pressState = currState[pressID]*1e-5; //Pa to bar

      cf_assert(tempState  >= 0 );
      cf_assert(pressState  >= 0 );

      tableInterpolate(tempState, pressState, itBin, kappa, ems);
 
      cf_assert(kappa >= 0 && "negative interpoldated absorption coeff");
      cf_assert(ems >= 0 && "negative interpoldated emmission coeff");

      CFreal wav = (itBin == minIdx) ? m_wavMin : m_wavTable[itBin];
             wav = (itBin == maxIdx) ? m_wavMax : m_wavTable[itBin];
      
      m_data(itState, itBin*3 + 0 ) = wav;
      m_data(itState, itBin*3 + 1 ) = ems; //*1e4; // converted from cm^3 to m^3 
      m_data(itState, itBin*3 + 2 ) = kappa;

      //CFLog(INFO, "out! \n";)

    }
  }

}


} // namespace RadiativeTransfer

} // namespace 
