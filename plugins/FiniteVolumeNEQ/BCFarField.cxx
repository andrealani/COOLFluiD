#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolumeNEQ/BCFarField.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include <iostream>
#include <vector>


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BCFarField,
                      CellCenterFVMData,
                      FiniteVolumeNEQModule>
aBCFarFieldProvider("BCFarField");

//////////////////////////////////////////////////////////////////////////////

void BCFarField::defineConfigOptions(Config::OptionList& options)
{
 
  options.addConfigOption< std::string >("FarFieldFile","Name of File containing the Farfield.");

  // options.addConfigOption< CFreal >("MeanDensity","Meanflow density (constant for incompressible flows.)");
  // options.addConfigOption< bool >("Interpolate","Interpolate the meanflow to the LEE grid?)");
  // options.addConfigOption< bool >("Write","Write the meanflow to file?)");
  // options.addConfigOption< bool >("WriteCT","Write the connectivity table between grids?)");
  
}

//////////////////////////////////////////////////////////////////////////////
BCFarField::BCFarField(const std::string& name) :
  FVMCC_BC(name)
  //,socket_states("states"),
  // socket_meanflow("meanflow"),
  // socket_fluent_states("fluent_states"),
  //socket_fluent_coords("fluent_coords"),
  // nbNodesPerElemFluent()
{
  addConfigOptionsTo(this);
  
  m_nameInputFile = "FarField.plt";
  setParameter("FarFieldFile",&m_nameInputFile);

  //  setParameter("MeanFlowFile",&m_file_name);
  
  //  MeanDensity = 0.0;
  //setParameter("MeanDensity",&MeanDensity);
  
  //_Interpolate = false;
  //setParameter("Interpolate",&_Interpolate);
  
  //  _Write = false;
  //setParameter("Write",&_Write);
  
  //_WriteCT = false;
  // setParameter("WriteCT",&_WriteCT);
}

//////////////////////////////////////////////////////////////////////////////

BCFarField::~BCFarField()
{
}

//////////////////////////////////////////////////////////////////////////////

//std::vector<Common::SafePtr<BaseDataSocketSink> >
//BCFarField::needsSockets()
//{
// std::vector<Common::SafePtr<BaseDataSocketSink> > result;
// result.push_back(&socket_states);
//  return result;
//}

//////////////////////////////////////////////////////////////////////////////

void BCFarField::setup()
{
FVMCC_BC::setup();

  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(m_nameInputFile);
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

  ifstream& file = fhandle->open(fname);

//  Better to declare it at the CFcase
  CFuint nbStatesFile = 30;
  CFuint ID_file;
  CFreal Ar_Coords[2];
  CFreal Ar_States[3];

  for(CFuint i=0; i < nbStatesFile; i++){

  file >> ID_file;
  file >> Ar_Coords[0] >> Ar_Coords[1];
  file >> Ar_States[0] >> Ar_States[1] >> Ar_States[2];

  V_Coords.push_back(vector<CFreal>() );
  V_Coords[i].push_back(Ar_Coords[0]);
  V_Coords[i].push_back(Ar_Coords[1]);

  V_States.push_back(vector<CFreal>() );
  V_States[i].push_back(Ar_States[0]);
  V_States[i].push_back(Ar_States[1]);
  V_States[i].push_back(Ar_States[2]);
  }

  fhandle->close();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "BCFarfield::execute() called for TRS: "
	<< trs->getName() << "\n");

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = getMethodData().getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.isBFace = true;
  geoData.trs = trs;

  CFreal error= 1.5e-3;

  const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    CFLogDebugMed( "iFace = " << iFace << "\n");

//    CF_DEBUG_POINT;

    // build the GeometricEntity
    geoData.idx = iFace;

    GeometricEntity *const face = geoBuilder->buildGE();
//    setGhostState(face);;

//  CF_DEBUG_POINT;

  State& innerState = *face->getState(0);
//  State& ghostState = *face->getState(1);

  CFreal CoordX = innerState.getCoordinates()[XX];
  CFreal CoordY = innerState.getCoordinates()[YY];
//  CFuint LocalID = innerState.getLocalID();

//  CF_DEBUG_OBJ(LocalID);

  CFreal Diffx=0., Diffy=0.;

    for (CFuint k=0 ; k<nbStatesFile; k++ ){
      Diffx= abs(CoordX-V_Coords[k][0]);
      Diffy= abs(CoordY-V_Coords[k][1]);

      if ((Diffx < error) && (Diffy < error)){
//      cout << "Found the correct matching. The file ID is: " << k << " the FaceID: " << iFace << " and the LocalID is : " << LocalID << endl;
//      cout << V_States[k][0] << " " << V_States[k][1] << " " << V_States[k][2] << endl;

      V_Pairing.push_back(k);
      
//      cout << V_Pairing[iFace] << endl;

      break;
        
      }
    }

  geoBuilder->releaseGE();

  }

//////////////////////////////////////////////////////////////////////////////
// Here I get the partial densities from mutation for every element in the field

//  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
//  MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

//  const CFuint nbCells = cells->nbRows();

//  for(CFuint i_Cell = 0; i_Cell < nbCells; i_Cell++){
  
//    geoData.idx = i_Cell;

//    SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
//    SafePtr<CellCenterFVM> fvmcc = spaceMethod.d_castTo<CellCenterFVM>();
//    cf_assert(fvmcc.isNotNull());
  
//    m_cellBuilder = fvmcc->getData()->getCellTrsGeoBuilder();

//    GeometricEntity *const element = m_cellBuilder->buildGE();

//    State& everyState = *element->getState(0);
//    CFuint LocalID = everyState.getLocalID();

//    cout << LocalID << " " << everyState[0] << endl;

  
//  }


//////////////////////////////////////////////////////////////////////////////

// DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  //DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  //meanflow.resize(states.size());

  //const CFuint nb_vars = PhysicalModelStack::getActive()->getNbEq();

  //for (CFuint i = 0; i < states.size(); ++i)
  // meanflow[i].resize(nb_vars);
  
//  for (CFuint i_check=0; i_check < nbTrsFaces; i_check++){
//   cout << V_Pairing[i_check] << endl;
//   }

}

//////////////////////////////////////////////////////////////////////////////

//void BCFarField::executeOnTrs()
//{
      //CFAUTOTRACE;

  
      // SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
      //CFLogDebugMin( "BCFarfield::execute() called for TRS: "
      //	<< trs->getName() << "\n");
      //CF_DEBUG_POINT;
      // Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
      // geoBuilder = getMethodData().getFaceTrsGeoBuilder();

      // SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
      // geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

      // FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
      // geoData.isBFace = true;
      //  geoData.trs = trs;
      // CF_DEBUG_POINT;
      //const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
      //for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
      //CFLogDebugMed( "iFace = " << iFace << "\n");

    // build the GeometricEntity
    //geoData.idx = iFace;

      // GeometricEntity *const face = geoBuilder->buildGE();
      //CF_DEBUG_POINT;
      //setGhostState(face);
      //CF_DEBUG_POINT;
      //geoBuilder->releaseGE();

      // }

  // SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  //  CFLogDebugMin( "Meanflow::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  //  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  //  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  // DataHandle<RealVector> fluent_states = socket_fluent_states.getDataHandle();
  // DataHandle<RealVector> fluent_coords = socket_fluent_coords.getDataHandle();
  
  // const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  // const CFuint nbstates=PhysicalModelStack::getActive()->getNbEq();

/* read data from Fluent-Tecplot file */

  // boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(m_file_name);
  // Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

  // ifstream& file = fhandle->open(fname);

  // std::string dummy;
  // char ch;
  //CFuint NNode = 0;

/* header */
  // getline (file,dummy);


  // for(int i=1; i< 2*nbdim + 2; i++) {
  //	getline (file,dummy);
  // }

  // getline (file,dummy);
  // getline (file,dummy);

  // file >> ch >> ch; 
  //file >> NNode;
  
  // fluent_coords.resize(NNode);
  // fluent_states.resize(NNode);

  //for (CFuint i = 0; i < NNode; ++i) {
  //fluent_states[i].resize(nbstates);
  //fluent_coords[i].resize(nbdim);
  // }

  // for(CFint i=1; i<4; i++)
  //getline (file,dummy);

/* variables */
  // for(CFuint dim=0; dim<nbdim; dim++) {
  //for(CFuint i=0; i<NNode; i++)
  //  file >> fluent_coords[i][dim];
  // }
  
  //density is not loaded
  // for(CFuint st=1; st<nbstates; st++) {
  ///for(CFuint i=0; i<NNode; i++)
  // file >> fluent_states[i][st];
  // }

/* connectivity */
//   nbNodesPerElemFluent.resize(NElem);
//   for(CFuint i=0; i<NElem; i++)
//     nbNodesPerElemFluent[i]=3;
// 
//   Common::ConnectivityTable<CFuint> FluentMesh;
//   FluentMesh.resize(nbNodesPerElemFluent);
// 
// 
// 
//   for(CFuint i=0; i<NElem; i++) {
//     for(CFuint j=0; j<3; j++)
//       file >> (FluentMesh) (i, j);
//   }

  // fhandle->close();

  //if (_Interpolate) {
  //  CTable.resize(states.size());
  
/* interpolate data to the RDS grid */
  // for (CFuint iState = 0; iState < states.size(); ++iState) {
  // Node& coord = states[iState]->getCoordinates();
  //  CFreal mindist = 10000.0;
  // for (CFuint FluentState = 0; FluentState < NNode; ++FluentState) {
  //   RealVector d=fluent_coords[FluentState]-coord;
  //   CFreal dist=d.norm2();
  //  if(dist<mindist) {
  //	CTable[iState] = FluentState;
  //    mindist = dist;
  //    RealVector& meanflow_state = meanflow[iState];
  //    meanflow_state[nbdim+1] = fluent_states[FluentState][nbdim+1]+101325.;
  //	meanflow_state[0] = MeanDensity;
  //	for (CFuint iDim = 1; iDim < nbdim+1; iDim++)
  //	  meanflow_state[iDim] = fluent_states[FluentState][iDim];
  //  }
  //}
  // }
  
// write connectivity table for the source interpolation
 
// if(_WriteCT) {
//  writeConnectivityFile();
// }
  //
  //}
  //else {
  //for (CFuint iState = 0; iState < states.size(); ++iState) {
  //  RealVector& meanflow_state = meanflow[iState];
  //  meanflow_state[nbdim+1] = fluent_states[iState][nbdim+1]+101325.;
  //  meanflow_state[0] = MeanDensity;
  //	for (CFuint iDim = 1; iDim < nbdim+1; iDim++)
  //	  meanflow_state[iDim] = fluent_states[iState][iDim];
  //  }
  //}

  //if(_Write) {
  // writeMeanflowFile();
  //}


  // SafePtr<LinEulerTerm> lterm = PhysicalModelStack::getActive()->getImplementor()->	getConvectiveTerm().d_castTo<LinEulerTerm>();

  // lterm->setMeanFlowArray(meanflow);
  
  
      //}

//////////////////////////////////////////////////////////////////////////////
void BCFarField::setGhostState(Framework::GeometricEntity *const face){

  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  CFuint FaceID = ghostState.getLocalID();
//  CF_DEBUG_OBJ(FaceID);


  for (CFuint i=0; i < 5; i++){
    ghostState[i]= innerState[i];
  }


  ghostState[5] = max(40.0, ( 2*V_States[V_Pairing[FaceID]][0]-innerState[5] ) );
  ghostState[6] = 2*V_States[V_Pairing[FaceID]][1] - innerState[6];
  ghostState[7] = max(1500.0, ( 2*V_States[V_Pairing[FaceID]][2]-innerState[7] ) );

  for (CFuint k=0; k < 8; k++){
//  CF_DEBUG_OBJ(innerState[k]);
//  CF_DEBUG_OBJ(ghostState[k]);
  }

}

//////////////////////////////////////////////////////////////////////////////
void BCFarField::unsetup()
{


  // DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  // DataHandle<RealVector> fluent_states = socket_fluent_states.getDataHandle();
  // DataHandle<RealVector> fluent_coords = socket_fluent_coords.getDataHandle();

  // meanflow.resize(0);
  // fluent_states.resize(0);
  // fluent_coords.resize(0);
}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

