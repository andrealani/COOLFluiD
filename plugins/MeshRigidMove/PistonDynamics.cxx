#include "Common/PE.hh"
#include "PistonDynamics.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "MeshRigidMove/MeshRigidMove.hh"
#include "RigidMoveData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Node.hh"
#include "MathTools/MathConsts.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathChecks.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PistonDynamics, RigidMoveData, MeshRigidMoveModule> pistonDynamicsProvider("PistonDynamics");

//////////////////////////////////////////////////////////////////////////////
 
std::vector<Common::SafePtr<BaseDataSocketSink> >
 PistonDynamics::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_volumes);
  return result;
}
//////////////////////////////////////////////////////////////////////////////

void PistonDynamics::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< CFreal >("Dpiston","Diameter of the piston");
 options.addConfigOption< CFreal >("Lpiston","Length of the piston");
 options.addConfigOption< CFreal >("Mpiston","Mass of the piston");
 options.addConfigOption< CFreal >("Xmin","Mini abs");
 options.addConfigOption< CFreal >("Xmax","Max_abs");
 options.addConfigOption< CFreal >("Xleft","Xleft");
 //options.addConfigOption< std::string  >("BackFace","Back Face of Piston");
}

//////////////////////////////////////////////////////////////////////////////
PistonDynamics::PistonDynamics(const std::string& name) :
   RigidMoveCom(name),
   socket_nodes("nodes"),
   socket_states("states"),
   socket_gstates("gstates"),
   socket_volumes("volumes"),
   _faceTrsGeoBuilder(),
   _geoWithNodesBuilder(),
   _Acceleration(0.0),
   _Pfront(0.0),
   _Pback(100.0),
   _Vpiston(0.0),
   _LogP(0.0)
{
  addConfigOptionsTo(this);
  D_piston = 0.01;
  setParameter("Dpiston",&D_piston);
  L_piston = 0.03;
  setParameter("Lpiston",&L_piston);
  M_piston = 1.;
  setParameter("Mpiston",&M_piston);
  Min_piston = 0.1;
  setParameter("Xmin",&Min_piston);
  Max_piston = 1.;
  setParameter("Xmax",&Max_piston);
  _totalDisplacement   = 1.;
  setParameter("Xleft",&_totalDisplacement);
}

//////////////////////////////////////////////////////////////////////////////

PistonDynamics::~PistonDynamics()
{
}
//////////////////////////////////////////////////////////////////////////////////

void PistonDynamics::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes  = socket_nodes.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();
  RealVector coord(nbDim);

  const CFreal dt = SubSystemStatusStack::getActive()->getDTDim();

  // Compute acceleration, speed, displacement of piston
  _Acceleration = getAcceleration();
  const CFreal speed = _Vpiston + _Acceleration*dt;
   CFreal displacement = speed*dt;
   CFreal PistL = _totalDisplacement - Min_piston;
   CFreal PistR =  (_totalDisplacement + L_piston);
   CFout << "Piston Speed: " << speed <<" Right face position  "<< PistR << "Limit "<< Max_piston  <<"\n";
   const CFreal Limit = Max_piston-1e-3;
  for (CFuint i = 0; i < nbNodes; ++i) {
    for (CFuint j = 0; j < nbDim; ++j) {
      // Get the old coordinates
      coord[j] = (*nodes[i])[j];
    }

    
  if (PistR>Limit){ 
  CFout << "Piston has reached its limit postiion " << PistR <<"\n";
   assert(PistR<Limit);
   }

 
    /// Test modification (piston testcase)
//    if (acceleration > 0){
        if((coord[0] > Min_piston) && (coord[0] < PistR)){
           coord[0] += displacement*(coord[0]-Min_piston)/PistL; 
         }
     
       if((coord[0] >= PistR) && (coord[0] < Max_piston)){
          coord[0] += displacement*(Max_piston-coord[0])/(Max_piston-PistR); 
        }
      //  if((coord[0] > Min_piston) && (coord[0] < PistR)){
      //     coord[0] += displacement*(Max_piston-coord[0])/(Max_piston-PistR);
      //   }
//     }
  
//    if (acceleration < 0){
//        if((coord[0] > Min_piston) && (coord[0] < PistR)){
//          coord[0] += displacement*(coord[0]-Min_piston)/(PistL-Min_piston); 
 //        }
     
  //     if((coord[0] >= PistR) && (coord[0] < Max_piston)){
 //         coord[0] += displacement*(Max_piston-coord[0])/(Max_piston-PistR); 
   //     }
   // }



    for (CFuint j = 0; j < nbDim; ++j) {
      // Set the new coordinates
     // if (coord[j] <= -0.5){
         (*nodes[i])[j] = coord[j];
    //} //end oord[j] <= 3
   }
  }

  // Compute the new position of piston
  _totalDisplacement += displacement;
  _Vpiston = speed; 
  WritePistonData();
  
  // recompute volumes
  computeVolumes();
}

///////////////////////////////////////////////////////////////////////////////////////
void PistonDynamics::WritePistonData()
{
// only the first processor writes Piston Data
  PE::GetPE().setBarrier();
    if (PE::GetPE().GetRank() == 0) {   
   CFreal PistL = _totalDisplacement - Min_piston;
   CFreal PistR = PistL+L_piston;
   const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
   const CFreal Time = SubSystemStatusStack::getActive()->getCurrentTimeDim();
    ofstream fout("PistonData.plt", ios::app);
   if (iter <1){
      fout << "TITLE = Unstructured grid data" << endl;
      fout << "VARIABLES  =  \"Time\" \"Position\" \"PL\" \"PR\"   \"Acceleration\"  \"Speed\" \"Pressure_Back\" \"Pressure_Front\"   " << endl;
      fout.precision(14);fout << Time   << "  "   << _totalDisplacement << "  " << PistL << "  "<< PistR << "  " << _Acceleration << " "  << _Vpiston    << "  " << _Pback << "   " <<  _Pfront   << endl;
      //fout.precision(14);fout << Time   << "  "   << _totalDisplacement << "  " << PistL << "  "<< PistR << "  " << PistL-PistR << " "  << _Vpiston    << "  " << _Pback << "   " <<  _Pfront   << endl;
     }
     else {
      fout.precision(14); fout << Time   << "  "   << _totalDisplacement << "  " << PistL << "  "<< PistR << "  " << _Acceleration << " "   << _Vpiston    << "  " << _Pback << "    " << _Pfront << endl;
      //fout.precision(14); fout << Time   << "  "   << _totalDisplacement << "  " << PistL << "  "<< PistR << "  " << PistL-PistR << " "   << _Vpiston    << "  " << _Pback << "    " << _Pfront << endl;
    }
      fout << " " << endl;
}
 PE::GetPE().setBarrier();
}
////////////////////////////////////////////////////////////////////////////////////////////////
void PistonDynamics::WriteSpaceTimeData(const CFreal& Xpostion,
  				        const CFreal& LogP,
                                              CFuint& Imax, 
                                              CFuint& nbface, 
                                       std::string filename)  

{
       const  CFreal m_saveRate = 3;
       const CFreal Time = SubSystemStatusStack::getActive()->getCurrentTimeDim();
       const CFuint nbiter = SubSystemStatusStack::getActive()->getNbIter();

       long Jmax = static_cast<CFuint>(nbiter)+1;
              ofstream fout(filename.c_str(),ios::app); 
        //    ofstream fout(filename.c_str(),ios_base::app); 
       if ((nbiter == 0) && (nbface <1)){
             // ofstream fout(filename.c_str(),ios::trunc); 
     	fout << "TITLE = Space-Time-Diagram.plt" << endl;
        fout << "VARIABLES  =  \"Position\" \"Time\" \"LogP\"  " << endl;
        fout << "ZONE  T= \"SLUG-0\" "<< endl;
        fout << "I=" << Imax  << ","<< "\n";
        fout << "J=" << Jmax << "      ,"<< "\n";
        fout << "K=" << 1 <<"," <<"\n";
        fout << "F=POINT   "<< "\n";
        fout.precision(14);
        fout << Xpostion   << "  "   << Time << "  " << LogP << " " << endl;
        }
       else {
        fout.precision(14); 
        fout << Xpostion   << "  "   << Time << "  " << LogP << " " << endl;
   	 }
        //fout.close();
       //if (nbface == Imax-1){
       //if (nbF == Imax-1){
       if ((nbiter >=0) && (nbface ==0)){
        fstream fout(filename.c_str());
        fout.seekp(99,std::ios::beg);
        fout << Jmax;
       // fout.close();
      }
     
}
////////////////////////////////////////////////////////////////////////////////////////////////
CFreal PistonDynamics::getAcceleration()
{
//Get the pressure at the bac of the piston
// _Pback

// Get teh pressure at teh front of the piston
// _Pfront

  // The frontal piston area
    const CFreal A_piston = 3.14159 * D_piston* D_piston * 0.25;
  // The pressure forces
   const CFreal  Fp = A_piston *  getdeltaP();
  // The effective frontal aero of the seal
    const CFreal A_seal = 3.14159 * D_piston * L_piston;
  // Definition of the velocity tolerence
    const CFreal Vtol = 1e-6;
  // Definition of the coefficient of friction 
  // of the seal material on the tube wall
    const CFreal muf = 0.2; 
  // The maximum magnitude of frictional force is 
    const CFreal Ffmax = muf*A_seal*_Pfront;
  // The friciton force:
    CFreal Ff = (MathFunctions::signum(_Vpiston))*Ffmax;
  //      cout << "Ff" << Ff << endl;
   if (((std::abs(_Vpiston)) < Vtol) && ((std::abs(Fp))< Ffmax)){
        Ff = -1.0*A_piston*(_Pback -_Pfront); 
   //     cout << "Ffin" << Ff << endl;
   }
     //CFreal Acceleration = 1/M_piston*(Fp+Ff)
     //;
    return ((1/M_piston)*(Fp));
 
   //return Acceleration;
} 

//////////////////////////////////////////////////////////////////////////////
CFreal PistonDynamics::getdeltaP()
{
 
      CFuint count = 0; 
      CFuint count1 = 0;
     // we assume we are dealing with FiniteVolume
     //Loop over the TRS's and add `the "TRSName" + "-boundaryNormals" datasocketsink to the _dynamicSockets
    vector<SafePtr<TopologicalRegionSet> >&  trsList = getTrsList();
     for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS) { 

         SafePtr<TopologicalRegionSet> currTrs = trsList[iTRS];
         // create locally 
         const std::string trsName = currTrs->getName(); 

         SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = _faceTrsGeoBuilder.getGeoBuilder();
         geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
         FaceTrsGeoBuilder::GeoData& geoData = _faceTrsGeoBuilder.getDataGE();
         geoData.trs = trsList[iTRS];
  
    ///   Number of face in TRS on the partition or Total ?????????????????? 
         CFuint nbFaces = trsList[iTRS]->getLocalNbGeoEnts();
         //cout << "nbFaces Local" << nbFaces  << endl;
       //const CFuint nbFacesG = trsList[iTRS]->getGlobalNbGeoEnts();
        CFuint nbFacesG = 0;
        
        CFreal av_perssure = 0;
        CFreal  Total_av_perssure = 0;
         // cout << "nbace" << nbFaces << endl;
        for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
 
        // build ths current face
        geoData.idx = iFace;
        GeometricEntity *const currFace = _faceTrsGeoBuilder.buildGE();
        const CFuint faceID = currFace->getID();
    
        const vector<Node*>& nodesInFace = *currFace->getNodes();
        const State& lState = *currFace->getState(0);
        const State& rState = *currFace->getState(1);
        const Node& lNode = lState.getCoordinates();
        const Node& rNode = rState.getCoordinates();
                     
        //av_perssure = lState[0];
        if (lState[0]>0){
        av_perssure += lState[0];
        }
        else{
        av_perssure = 0.;
         }
         
 // if divide by save rate is integer then it is tim,e to save
//if (!(SubSystemStatusStack::getActive()->getNbIter() % m_saveRate))
      //  if ((iTRS > trsList.size()-3) && (currFace->getState(0)->isParUpdatable())&& (lState[0]>0)) {
        if ((iTRS > trsList.size()-3) && (lState[0]>0)) {
         const  CFreal LogP = std::log(lState[0]);
         const  CFreal Xposition = lNode[0];
         CFreal PistL = _totalDisplacement; // - Min_piston;
         //if (Xposition<=PistL){
         if (iTRS == 2){
          std::string NameL = "Space-Time-Diagram_Left_" + Common::StringOps::to_str(PE::GetPE().GetRank()) + ".plt";   
          WriteSpaceTimeData(Xposition,LogP,nbFaces,count,NameL);
         count++;        
 	   }
         else{   
          std::string NameR = "Space-Time-Diagram_Right_" + Common::StringOps::to_str(PE::GetPE().GetRank()) + ".plt";   
          WriteSpaceTimeData(Xposition,LogP,nbFaces,count1,NameR);
           count1++;        
           }
           
        }  
        // release the face
        _faceTrsGeoBuilder.releaseGE();
        }
       
      if (iTRS < trsList.size()-1) {
       if (PE::GetPE().GetProcessorCount() > 1) {
        MPI_Allreduce(&av_perssure, &Total_av_perssure, 1, MPI_DOUBLE, MPI_SUM,PE::GetPE().GetCommunicator());  
       MPI_Allreduce(&nbFaces, &nbFacesG, 1, MPI_INTEGER, MPI_SUM,PE::GetPE().GetCommunicator());  
           }
       else {
          Total_av_perssure = av_perssure;
          nbFacesG = nbFaces;  
         }
       if (iTRS==1){
       //   Total_av_perssure = av_perssure;
        // cout << "nbFaces TRS>0" << nbFaces << endl;
       //_Pfront = Total_av_perssure/1.;
       _Pfront = Total_av_perssure/nbFacesG;
        cout <<  trsName  << _Pfront << endl; 
        }
       if (iTRS == 0) {
        // Total_av_perssure = av_perssure;
      //   cout << "nbFaces" << nbFaces << endl;
       //_Pback =  Total_av_perssure/1.;
       _Pback =  Total_av_perssure/nbFacesG;
        cout <<  trsName << _Pback <<endl; 
       }
    }
 }     
       CFreal deltaP = _Pback - _Pfront;
       return deltaP;  
}

//////////////////////////////////////////////////////////////////////////////

void PistonDynamics::setup()
{
  _faceTrsGeoBuilder.setup();
  _geoWithNodesBuilder.setup();
}

//////////////////////////////////////////////////////////////////////////////

void PistonDynamics::computeVolumes()
{
  // this has to be consistent with StdSetup
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
   getTrs("InnerCells");

  TrsGeoWithNodesBuilder::GeoData& geoData = _geoWithNodesBuilder.getDataGE();
  geoData.trs = cells;
  
  CFuint countNegativeVol = 0;
  const CFuint nbElems = cells->getLocalNbGeoEnts();
  for (CFuint iElem = 0; iElem < nbElems; ++iElem) {
    // build the GeometricEntity
    geoData.idx = iElem;
    GeometricEntity *const cell = _geoWithNodesBuilder.buildGE();
    volumes[iElem] = cell->computeVolume();
    
    if (volumes[iElem] < 0.0)
      {
      countNegativeVol++;
      cout.precision(14);
      CFout << "Cell [" << iElem << "] with [" << cell->nbNodes() << "] nodes has negative volume [" << volumes[iElem] << "]\n";
      volumes[iElem] = 0.0;
    }

    if ( MathChecks::isZero(volumes[iElem]) )
    {
      cout.precision(14);
      CFout << "Cell [" << iElem << "] with [" << cell->nbNodes() << "] nodes has zero volume [" << volumes[iElem] << "]\n";
      // print coordinates
      for (CFuint i = 0; i < cell->nbNodes(); ++i) {
       CFout << *cell->getNode(i) << ", ";
      } CFout << "\n";
      volumes[iElem] = 0.0;
    }
    
    cf_assert(volumes[iElem] > 0.);
    
    //release the GeometricEntity
    _geoWithNodesBuilder.releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD
