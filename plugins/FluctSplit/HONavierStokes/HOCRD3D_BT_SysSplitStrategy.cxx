#include "Common/BadValueException.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/HONavierStokes/FluctSplitHONavierStokes.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/HONavierStokes/HOCRD3D_BT_SysSplitStrategy.hh"

#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/EulerTerm.hh"

#define BONANNI_DEBUG 0
#define BONANNI_DEBUG_NEW 0

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

template < int ID0, int ID1, int ID2 >
void compute_normal_to_face ( const std::vector<State*>& states, RealVector& normal )
{
  CFreal x0X;
  CFreal x0Y;
  CFreal x0Z;
  CFreal x1X;
  CFreal x1Y;
  CFreal x1Z;
  CFreal x2X;
  CFreal x2Y;
  CFreal x2Z;
  

    x0X = states[ID0]->getCoordinates()[XX];
    x0Y = states[ID0]->getCoordinates()[YY];
    x0Z = states[ID0]->getCoordinates()[ZZ];
      

   x1X = states[ID1]->getCoordinates()[XX];
      x1Y = states[ID1]->getCoordinates()[YY];
      x1Z = states[ID1]->getCoordinates()[ZZ];

     x2X = states[ID2]->getCoordinates()[XX];
     x2Y = states[ID2]->getCoordinates()[YY];
     x2Z = states[ID2]->getCoordinates()[ZZ];


  const CFreal v1x = x1X - x0X;
  const CFreal v1y = x1Y - x0Y;
  const CFreal v1z = x1Z - x0Z;

  const CFreal v2x = x2X - x1X;
  const CFreal v2y = x2Y - x1Y;
  const CFreal v2z = x2Z - x1Z;

  normal[XX] = 0.5 * (  v1y*v2z - v1z*v2y );
  normal[YY] = 0.5 * ( -v1x*v2z + v1z*v2x );
  normal[ZZ] = 0.5 * (  v1x*v2y - v1y*v2x );

}

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRD3D_BT_SysSplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHONavierStokesModule>
aHOCRD3DBT_Provider("HOCRD3D_BT");



//////////////////////////////////////////////////////////////////////////////

void HOCRD3D_BT_SysSplitStrategy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal>("Delta","Delta of variable.");
  options.addConfigOption< CFreal>("Length","Reference Length.");
  options.addConfigOption< CFreal>("Speed","Reference Speed.");
  options.addConfigOption< std::string>("VarName","Variable name.");
  options.addConfigOption< bool >  ("StoreThetas","Store the thetas for visualization");
  options.addConfigOption< bool >  ("UmaxTheta","Use the maximum of the thetas");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("MinTheta","Minimum theta, used to keep a minimum diffusion");
  options.addConfigOption< std::string>("Shockdetector","Which shock detetecto to use");

}


//////////////////////////////////////////////////////////////////////////////


//.................... MODIFIED ...............................
HOCRD3D_BT_SysSplitStrategy::HOCRD3D_BT_SysSplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_thetas("thetas",false),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_scaledFaceNormals(),
  matrix_face_norms(4,3),
  matrix_node_norms(4,3),
  vector_face_areas(4),
  vector_node_areas(4)
{
  addConfigOptionsTo(this);
   const CFuint dim = DIM_3D;
   const CFuint nbfaces = 4;
   const CFuint nbnodes = 4;

   /// @todo this is a memory leak, they are not destroyed in the destructor of InwardNormalsData
   CFreal * faceNormals = new CFreal[nbfaces*dim];
   CFreal * faceAreas   = new CFreal[nbfaces];
   CFreal * nodeNormals = new CFreal[nbnodes*dim];
   CFreal * nodeAreas   = new CFreal[nbnodes];
   // place them in the inwardnormals
   m_subcell_normals =
     new InwardNormalsData(faceNormals,faceAreas,nodeNormals,nodeAreas,nbfaces,0);

   _deltaP = 0.0;
  setParameter("Delta",&_deltaP);

  _length = 1.0;
  setParameter("Length",&_length);

  _speed = 0.0;
  setParameter("Speed",&_speed);

  _varName = "p";
  setParameter("VarName",&_varName);


  m_store_thetas = false;
  setParameter("StoreThetas",&m_store_thetas);

  m_use_max_theta = true;
  setParameter("UmaxTheta",&m_use_max_theta);


  m_min_theta = 0.;
  setParameter("MinTheta",&m_min_theta);


  // Choice between Jirka shock capturing and 
  // The improved one of Antonino, the one of Antonino
  // is implemented only in 2D
  _sh_detector = "Jirka";
  this->setParameter("Shockdetector",&_sh_detector);


}
//.............................................................




//////////////////////////////////////////////////////////////////////////////

HOCRD3D_BT_SysSplitStrategy::~HOCRD3D_BT_SysSplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD3D_BT_SysSplitStrategy::unsetup()
{
  deletePtr(m_subcell_normals);


  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }
for (CFuint i = 0; i < m_k1Plus.size(); ++i)
  {
    deletePtr(m_k1Plus[i]);
    deletePtr(m_k2Plus[i]);
    deletePtr(m_k3Plus[i]);
    deletePtr(m_k4Plus[i]);
    deletePtr(m_k5Plus[i]);
    deletePtr(m_k6Plus[i]);
    deletePtr(m_k7Plus[i]);

    deletePtr(m_k[i]);
    deletePtr(m_kMin[i]);
  }
  deletePtr(m_inverter);

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
      deletePtr(m_qdExtraVars[i]);
  }
 

}

//////////////////////////////////////////////////////////////////////////////

void HOCRD3D_BT_SysSplitStrategy::setup()
{



  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  if (getMethodData().isMultipleSplitter())
    throw Common::BadValueException(FromHere(), "Cannot use HOCRD on curved elements with multiple splitters");


  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i)
  {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }


  //   number of quadrature points used to compute the fluctuation
  //   const CFuint nbQdPts = 3;


  m_scaledFaceNormals.resize(4);
  for (CFuint  i = 0; i < m_scaledFaceNormals.size(); ++i)
  {
    m_scaledFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }



  Facenormal0.resize(PhysicalModelStack::getActive()->getDim());
  Facenormal1.resize(PhysicalModelStack::getActive()->getDim()); 
  Facenormal2.resize(PhysicalModelStack::getActive()->getDim());
  Facenormal3.resize(PhysicalModelStack::getActive()->getDim());
 // physical data evaluated in the quadrature points
  m_pdata.resize(10);
  for (CFuint  i = 0; i < 10; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  

  // table storing the states and residual of the sub-element
  substates.resize(4);   // 4 states in each sub element
  subresidual.resize(4); // 4 residuals in each sub element 
//.............................................................

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);
  subresidual[3].resize(nbEqs);
//.................... MODIFIED ...............................
  m_phisubT.resize(7); // P2 tetra has 7 sub tetrahedra 
  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);
  m_phisubT[4] = new RealVector(nbEqs);
  m_phisubT[5] = new RealVector(nbEqs);
  m_phisubT[6] = new RealVector(nbEqs);
  m_phiN1.resize(4);
  m_phiN2.resize(4);
  m_phiN3.resize(4);
  m_phiN4.resize(4);
  m_phiN5.resize(4);
  m_phiN6.resize(4);
  m_phiN7.resize(4);
  m_phi.resize(4);

  m_phiN1[0].resize(nbEqs);
  m_phiN1[1].resize(nbEqs);
  m_phiN1[2].resize(nbEqs);
  m_phiN1[3].resize(nbEqs);

  m_phiN2[0].resize(nbEqs);
  m_phiN2[1].resize(nbEqs);
  m_phiN2[2].resize(nbEqs);
  m_phiN2[3].resize(nbEqs);


  m_phiN3[0].resize(nbEqs);
  m_phiN3[1].resize(nbEqs);
  m_phiN3[2].resize(nbEqs);
  m_phiN3[3].resize(nbEqs);

  m_phiN4[0].resize(nbEqs);
  m_phiN4[1].resize(nbEqs);
  m_phiN4[2].resize(nbEqs);
  m_phiN4[3].resize(nbEqs);

  m_phiN5[0].resize(nbEqs);
  m_phiN5[1].resize(nbEqs);
  m_phiN5[2].resize(nbEqs);
  m_phiN5[3].resize(nbEqs);

  m_phiN6[0].resize(nbEqs);
  m_phiN6[1].resize(nbEqs);
  m_phiN6[2].resize(nbEqs);
  m_phiN6[3].resize(nbEqs);

  m_phiN7[0].resize(nbEqs);
  m_phiN7[1].resize(nbEqs);
  m_phiN7[2].resize(nbEqs);
  m_phiN7[3].resize(nbEqs);

  m_phi[0].resize(nbEqs);
  m_phi[1].resize(nbEqs);
  m_phi[2].resize(nbEqs);
  m_phi[3].resize(nbEqs);


//.............................................................

//.................... MODIFIED ...............................
  subelemtable.resize(7,4);     // sub element table : contain the faces of each sub-element
  subelemfacedir.resize(7,4); // contains the information about the sign of face normal
  subelem_nodetable.resize(7,4); // contains the nodes in each element

  subelemtable(0,0) = 16;             subelemfacedir(0,0) = -1.;
  subelemtable(0,1) = 1;              subelemfacedir(0,1) = 1.;
  subelemtable(0,2) = 14;             subelemfacedir(0,2) = 1.;
  subelemtable(0,3) = 10;             subelemfacedir(0,3) = 1.;

  subelemtable(1,0) = 6;              subelemfacedir(1,0) = 1.;
  subelemtable(1,1) = 18;             subelemfacedir(1,1) = -1.;
  subelemtable(1,2) = 13;             subelemfacedir(1,2) = 1.;
  subelemtable(1,3) = 9;              subelemfacedir(1,3) = 1.;

  subelemtable(2,0) = 19;             subelemfacedir(2,0) = 1.;
  subelemtable(2,1) = 0;              subelemfacedir(2,1) = 1.;
  subelemtable(2,2) = 16;             subelemfacedir(2,2) = 1.;
  subelemtable(2,3) = 11;             subelemfacedir(2,3) = 1.;

  subelemtable(3,0) = 18;             subelemfacedir(3,0) = 1.;
  subelemtable(3,1) = 7;              subelemfacedir(3,1) = 1.;
  subelemtable(3,2) = 21;             subelemfacedir(3,2) = 1.;
  subelemtable(3,3) = 8;              subelemfacedir(3,3) = 1.;

  subelemtable(4,0) = 20;             subelemfacedir(4,0) = -1.;
  subelemtable(4,1) = 21;             subelemfacedir(4,1) = -1.;
  subelemtable(4,2) = 19;             subelemfacedir(4,2) = -1.;
  subelemtable(4,3) = 15;             subelemfacedir(4,3) = 1.;

  subelemtable(5,0) = 5;              subelemfacedir(5,0) = 1.;
  subelemtable(5,1) = 2;              subelemfacedir(5,1) = 1.;
  subelemtable(5,2) = 12;             subelemfacedir(5,2) = 1.;
  subelemtable(5,3) = 17;             subelemfacedir(5,3) = -1.;

  subelemtable(6,0) = 17;             subelemfacedir(6,0) = 1.;
  subelemtable(6,1) = 3;              subelemfacedir(6,1) = 1.;
  subelemtable(6,2) = 20;             subelemfacedir(6,2) = 1.;
  subelemtable(6,3) = 4;              subelemfacedir(6,3) = 1.;


  // sub face table : contain the node that contain each face
  subfacetable.resize(22,3); // 22 sub faces with 3 states each

//..... External 1 .................... External 2 ..................... External 3 ..................... External 4 ...
//...... FACE 0 OK..........       //...... FACE 4 OK..........       //...... FACE 8 OK..........      //...... FACE 12 OK  ..........
  subfacetable(0,0) = 0;           subfacetable(4,0) = 0;           subfacetable(8,0) = 0;          subfacetable(12,0) = 3;
  subfacetable(0,1) = 4;           subfacetable(4,1) = 9;           subfacetable(8,1) = 6;          subfacetable(12,1) = 7;
  subfacetable(0,2) = 7;           subfacetable(4,2) = 8;           subfacetable(8,2) = 5;          subfacetable(12,2) = 8;

//...... FACE 1 OK ..........       //...... FACE 5 OK..........       //...... FACE 9 OK..........      //...... FACE 13 OK ..........
  subfacetable(1,0) = 1;           subfacetable(5,0) = 3;           subfacetable(9,0) = 2;          subfacetable(13,0) = 2;
  subfacetable(1,1) = 7;           subfacetable(5,1) = 8;           subfacetable(9,1) = 5;          subfacetable(13,1) = 8;
  subfacetable(1,2) = 4;           subfacetable(5,2) = 9;           subfacetable(9,2) = 6;          subfacetable(13,2) = 5;

//...... FACE 2 OK..........       //...... FACE 6 OK..........       //...... FACE 10 OK..........     //...... FACE 14 OK ..........
  subfacetable(2,0) = 3;           subfacetable(6,0) = 2;           subfacetable(10,0) = 1;         subfacetable(14,0) = 1;
  subfacetable(2,1) = 9;           subfacetable(6,1) = 6;           subfacetable(10,1) = 4;         subfacetable(14,1) = 5;
  subfacetable(2,2) = 7;           subfacetable(6,2) = 8;           subfacetable(10,2) = 5;         subfacetable(14,2) = 7;
//...... FACE 3 OK..........       //...... FACE 7 OK..........       //...... FACE 11 OK ..........     //...... FACE 15 OK ..........
  subfacetable(3,0) = 0;           subfacetable(7,0) = 0;           subfacetable(11,0) = 0;         subfacetable(15,0) = 5;
  subfacetable(3,1) = 7;           subfacetable(7,1) = 8;           subfacetable(11,1) = 5;         subfacetable(15,1) = 8;
  subfacetable(3,2) = 9;           subfacetable(7,2) = 6;           subfacetable(11,2) = 4;         subfacetable(15,2) = 7;

//..................... Internal ...................
//...... FACE 16 OK..........   //...... FACE 19 OK..........
  subfacetable(16,0) = 4;       subfacetable(19,0) = 0;
  subfacetable(16,1) = 5;       subfacetable(19,1) = 7;
  subfacetable(16,2) = 7;       subfacetable(19,2) = 5;

//...... FACE 17 OK ..........   //...... FACE 20 OK..........
  subfacetable(17,0) = 7;       subfacetable(20,0) = 0;
  subfacetable(17,1) = 8;       subfacetable(20,1) = 8;
  subfacetable(17,2) = 9;       subfacetable(20,2) = 7;

//...... FACE 18 OK..........   //...... FACE 21 OK..........
  subfacetable(18,0) = 5;       subfacetable(21,0) = 0;
  subfacetable(18,1) = 6;       subfacetable(21,1) = 5;
  subfacetable(18,2) = 8;       subfacetable(21,2) = 8;

//.............................................................

  // Tetra 0 OK
  subelem_nodetable(0,0) = 1;
  subelem_nodetable(0,1) = 5;
  subelem_nodetable(0,2) = 4;
  subelem_nodetable(0,3) = 7;
  
   // Tetra 1 OK
  subelem_nodetable(1,0) = 5;
  subelem_nodetable(1,1) = 2;
  subelem_nodetable(1,2) = 6;
  subelem_nodetable(1,3) = 8;
  
   // Tetra 2 OK
  subelem_nodetable(2,0) = 4;
  subelem_nodetable(2,1) = 5;
  subelem_nodetable(2,2) = 0;
  subelem_nodetable(2,3) = 7;
  
   // Tetra 3 OK
  subelem_nodetable(3,0) = 0;
  subelem_nodetable(3,1) = 5;
  subelem_nodetable(3,2) = 6;
  subelem_nodetable(3,3) = 8;
  
   // Tetra 4 OK
  subelem_nodetable(4,0) = 5;
  subelem_nodetable(4,1) = 7;
  subelem_nodetable(4,2) = 8;
  subelem_nodetable(4,3) = 0;
  
   // Tetra 5 OK
  subelem_nodetable(5,0) = 7;
  subelem_nodetable(5,1) = 8;
  subelem_nodetable(5,2) = 9;
  subelem_nodetable(5,3) = 3;
  
   // Tetra 6 OK
  subelem_nodetable(6,0) = 0;
  subelem_nodetable(6,1) = 8;
  subelem_nodetable(6,2) = 9;
  subelem_nodetable(6,3) = 7;
  
  faceflux.resize(22); // one flux per sub face


  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(nbEqs);

  // Surface quadrature rule (so on triangles)
  qd0.resize(4); // quadrature points per face
  qd1.resize(4); // quadrature points per face
  qd2.resize(4); // quadrature points per face


  
  qd0[0] = 0.333;      qd1[0] = 0.333;      qd2[0] = 0.333;
  qd0[1] = 0.6;      qd1[1] = 0.2;      qd2[1] = 0.2;
  qd0[2] = 0.2;      qd1[2] = 0.6;      qd2[2] = 0.2;
  qd0[3] = 0.2;      qd1[3] = 0.2;      qd2[3] = 0.6;
  

  wqd.resize(4); 

//....O(h^3)............................
  wqd[0] = -27.0/48.0;
  wqd[1] = 25.0/48.0;
  wqd[2] = 25.0/48.0;
  wqd[3] = 25.0/48.0;
//......................................
//--------------------------------------------------------------



//..... HERE I JUST ALLOCATE MEMORY ................................

  qdstates.resize(4); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State(); 
  qdstates[2] = new State(); 
   qdstates[3] = new State();

  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(10);
  for (CFuint i = 0; i < 10; ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }


  // NOTE: resize the beta's storages in the MethodData

  const CFuint nbSubTriangles = 7;  //..now I have 7 sub tetrahedra..
  const CFuint nbNodesInTriangle = 4;

  DistributionData::BetaMatrices& betaMats = getMethodData().getDistributionData().betaMats;
  betaMats.resize(nbSubTriangles);

  for (CFuint t = 0; t < nbSubTriangles; ++t) 
   {
    betaMats[t].resize(nbNodesInTriangle);

      for (CFuint i = 0; i < nbNodesInTriangle; ++i) 
       {
        betaMats[t][i].resize(nbEqs, nbEqs);
       }
   }

  RealMatrix currFaceMat = getMethodData().getDistributionData().FacesMat_tet;
  currFaceMat.resize(22, 3);

  getMethodData().getDistributionData().computeBetas = true;
//...........................................................



//....................Resizing ..........
   nf0.resize(3);
   nf1.resize(3);
   nf2.resize(3);
   nf3.resize(3);
   nf4.resize(3);
   nf5.resize(3);
   nf6.resize(3);
   nf7.resize(3);
   nf8.resize(3);
   nf9.resize(3);
   nf10.resize(3);
   nf11.resize(3);
   nf12.resize(3);
   nf13.resize(3);
   nf14.resize(3);
   nf15.resize(3);
   nf16.resize(3);
   nf17.resize(3);
   nf18.resize(3);
   nf19.resize(3);
   nf20.resize(3);
   nf21.resize(3);

   FacesAreas.resize(22);
   FacesMat.resize(22,3);


  m_kPlus.resize(4);
  m_k.resize(4);
  m_kMin.resize(4);
  m_eValues.resize(4);
  m_k1Plus.resize(4);
  m_k2Plus.resize(4);
  m_k3Plus.resize(4);
  m_k4Plus.resize(4);
  m_k5Plus.resize(4);
  m_k6Plus.resize(4);
  m_k7Plus.resize(4);

  m_eValuesP.resize(nbEqs);
  m_eValuesM.resize(nbEqs);
  m_rightEv.resize(nbEqs,nbEqs);
  m_leftEv.resize(nbEqs,nbEqs);
  for (CFuint i = 0; i < 4; ++i) {
    m_kPlus[i] = new RealMatrix(nbEqs, nbEqs);
    m_kMin[i] = new RealMatrix(nbEqs, nbEqs);
    m_k[i] = new RealMatrix(nbEqs, nbEqs);
    m_eValues[i] = new RealVector(nbEqs);
    m_k1Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k2Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k3Plus[i] = new RealMatrix(nbEqs,nbEqs );
    m_k4Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k5Plus[i] = new RealMatrix(nbEqs, nbEqs);
    m_k6Plus[i] = new RealMatrix(nbEqs,nbEqs );
    m_k7Plus[i] = new RealMatrix(nbEqs, nbEqs);


  }
  Fn0u0.resize(nbEqs);
   Fn1u1.resize(nbEqs);
   Fn2u2.resize(nbEqs);
   Fn3u3.resize(nbEqs);
   Fn0u4.resize(nbEqs);
   Fn1u4.resize(nbEqs);
   Fn1u5.resize(nbEqs);
   Fn2u5.resize(nbEqs);
   Fn0u6.resize(nbEqs);
   Fn2u6.resize(nbEqs);
   Fn1u7.resize(nbEqs);
   Fn3u7.resize(nbEqs);
   Fn3u8.resize(nbEqs);
   Fn2u8.resize(nbEqs);
   Fn3u9.resize(nbEqs);
   Fn0u9.resize(nbEqs);


  m_uTemp.resize(nbEqs);
  m_sumKplus.resize(nbEqs,nbEqs);
  m_invK.resize(nbEqs, nbEqs);
  m_inverter = MatrixInverter::create(nbEqs, false);
  m_adimNormal.resize(3);
  m_sumKplusU.resize(nbEqs);

  _cterm = Framework::PhysicalModelStack::getActive()->getImplementor()->
          getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>();

  _cterm->resizePhysicalData(_pData);

  _grad.resize(Framework::PhysicalModelStack::getActive()->getDim());
  // set _choiceVar to pressure if default value is required
  if (_varName == "rho") {
                _varID = Physics::NavierStokes::EulerTerm::RHO;
                         }
  else {
                _varID = Physics::NavierStokes::EulerTerm::P;
       }
  if ( m_store_thetas && !socket_thetas.isConnected())
        throw Common::BadValueException (FromHere(),"User required storing of thetas but socket is not connected");






}





///////////////////////////////////////////////////////////////////////////////////////////////////

void HOCRD3D_BT_SysSplitStrategy::computeFluctuation(vector<RealVector>& residual)
{

  DistributionData& distdata = getMethodData().getDistributionData();
   CellID = distdata.cellID;

   //.................... CELL STATES...............................
   vector<State*>& states = *distdata.states;
   //...............................................................
   InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[distdata.cellID]);

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i) { residual[i] = 0.0; }


  cf_assert (states.size() == 10); // P2 tetrahedra for solution space
  cf_assert (residual.size() >= states.size()); // state residual

//*******************************************************************************************
  
  
 // First we need to compute the normals of each face and put them in vectors
  compute_normal_to_face <0,4,7> ( states, nf0 ); 
  compute_normal_to_face <1,7,4> ( states, nf1 );   
  compute_normal_to_face <3,9,7> ( states, nf2 );  
  compute_normal_to_face <0,7,9> ( states, nf3 );   
  compute_normal_to_face <0,9,8> ( states, nf4 );  
  compute_normal_to_face <3,8,9> ( states, nf5 );    
  compute_normal_to_face <2,6,8> ( states, nf6 );    
  compute_normal_to_face <0,8,6> ( states, nf7 );    
  compute_normal_to_face <0,6,5> ( states, nf8 );   
  compute_normal_to_face <2,5,6> ( states, nf9 );    
  compute_normal_to_face <1,4,5> ( states, nf10 ); 
  compute_normal_to_face <0,5,4> ( states, nf11 );  
  compute_normal_to_face <3,7,8> ( states, nf12 );   
  compute_normal_to_face <2,8,5> ( states, nf13 );   
  compute_normal_to_face <1,5,7> ( states, nf14 );   
  compute_normal_to_face <5,8,7> ( states, nf15 );   
  compute_normal_to_face <4,5,7> ( states, nf16 );    
  compute_normal_to_face <7,8,9> ( states, nf17 );  
  compute_normal_to_face <5,6,8> ( states, nf18 );   
  compute_normal_to_face <0,7,5> ( states, nf19 );    
  compute_normal_to_face <0,8,7> ( states, nf20 );    
  compute_normal_to_face <0,5,8> ( states, nf21 );  


  FacesMat(0,XX) = nf0[XX];
  FacesMat(0,YY) = nf0[YY];
  FacesMat(0,ZZ) = nf0[ZZ];

  FacesMat(1,XX) = nf1[XX];
  FacesMat(1,YY) = nf1[YY];
  FacesMat(1,ZZ) = nf1[ZZ];

  FacesMat(2,XX) = nf2[XX];
  FacesMat(2,YY) = nf2[YY];
  FacesMat(2,ZZ) = nf2[ZZ];

  FacesMat(3,XX) = nf3[XX];
  FacesMat(3,YY) = nf3[YY];
  FacesMat(3,ZZ) = nf3[ZZ];

  FacesMat(4,XX) = nf4[XX];
  FacesMat(4,YY) = nf4[YY];
  FacesMat(4,ZZ) = nf4[ZZ];

  FacesMat(5,XX) = nf5[XX];
  FacesMat(5,YY) = nf5[YY];
  FacesMat(5,ZZ) = nf5[ZZ];

  FacesMat(6,XX) = nf6[XX];
  FacesMat(6,YY) = nf6[YY];
  FacesMat(6,ZZ) = nf6[ZZ];

  FacesMat(7,XX) = nf7[XX];
  FacesMat(7,YY) = nf7[YY];
  FacesMat(7,ZZ) = nf7[ZZ];

  FacesMat(8,XX) = nf8[XX];
  FacesMat(8,YY) = nf8[YY];
  FacesMat(8,ZZ) = nf8[ZZ];

  FacesMat(9,XX) = nf9[XX];
  FacesMat(9,YY) = nf9[YY];
  FacesMat(9,ZZ) = nf9[ZZ];

  FacesMat(10,XX) = nf10[XX];
  FacesMat(10,YY) = nf10[YY];
  FacesMat(10,ZZ) = nf10[ZZ];

  FacesMat(11,XX) = nf11[XX];
  FacesMat(11,YY) = nf11[YY];
  FacesMat(11,ZZ) = nf11[ZZ];

  FacesMat(12,XX) = nf12[XX];
  FacesMat(12,YY) = nf12[YY];
  FacesMat(12,ZZ) = nf12[ZZ];

  FacesMat(13,XX) = nf13[XX];
  FacesMat(13,YY) = nf13[YY];
  FacesMat(13,ZZ) = nf13[ZZ];

  FacesMat(14,XX) = nf14[XX];
  FacesMat(14,YY) = nf14[YY];
  FacesMat(14,ZZ) = nf14[ZZ];

  FacesMat(15,XX) = nf15[XX];
  FacesMat(15,YY) = nf15[YY];
  FacesMat(15,ZZ) = nf15[ZZ];

  FacesMat(16,XX) = nf16[XX];
  FacesMat(16,YY) = nf16[YY];
  FacesMat(16,ZZ) = nf16[ZZ];

  FacesMat(17,XX) = nf17[XX];
  FacesMat(17,YY) = nf17[YY];
  FacesMat(17,ZZ) = nf17[ZZ];

  FacesMat(18,XX) = nf18[XX];
  FacesMat(18,YY) = nf18[YY];
  FacesMat(18,ZZ) = nf18[ZZ];

  FacesMat(19,XX) = nf19[XX];
  FacesMat(19,YY) = nf19[YY];
  FacesMat(19,ZZ) = nf19[ZZ];

  FacesMat(20,XX) = nf20[XX];
  FacesMat(20,YY) = nf20[YY];
  FacesMat(20,ZZ) = nf20[ZZ];

  FacesMat(21,XX) = nf21[XX];
  FacesMat(21,YY) = nf21[YY];
  FacesMat(21,ZZ) = nf21[ZZ];
  

  distdata.FacesMat_tet= FacesMat;
    
  

// compute the areas of each face
  FacesAreas[0] =  sqrt(nf0[XX]*nf0[XX]+nf0[YY]*nf0[YY]+nf0[ZZ]*nf0[ZZ]) ;
   FacesAreas[1] = sqrt(nf1[XX]*nf1[XX]+nf1[YY]*nf1[YY]+nf1[ZZ]*nf1[ZZ]) ;
   FacesAreas[2] = sqrt(nf2[XX]*nf2[XX]+nf2[YY]*nf2[YY]+nf2[ZZ]*nf2[ZZ]) ;
   FacesAreas[3] = sqrt(nf3[XX]*nf3[XX]+nf3[YY]*nf3[YY]+nf3[ZZ]*nf3[ZZ]) ;
   FacesAreas[4] = sqrt(nf4[XX]*nf4[XX]+nf4[YY]*nf4[YY]+nf4[ZZ]*nf4[ZZ]) ;
   FacesAreas[5] = sqrt(nf5[XX]*nf5[XX]+nf5[YY]*nf5[YY]+nf5[ZZ]*nf5[ZZ]) ;
   FacesAreas[6] = sqrt(nf6[XX]*nf6[XX]+nf6[YY]*nf6[YY]+nf6[ZZ]*nf6[ZZ]) ;
   FacesAreas[7] = sqrt(nf7[XX]*nf7[XX]+nf7[YY]*nf7[YY]+nf7[ZZ]*nf7[ZZ]) ;
   FacesAreas[8] = sqrt(nf8[XX]*nf8[XX]+nf8[YY]*nf8[YY]+nf8[ZZ]*nf8[ZZ]) ;
   FacesAreas[9] = sqrt(nf9[XX]*nf9[XX]+nf9[YY]*nf9[YY]+nf9[ZZ]*nf9[ZZ]) ;
   FacesAreas[10] = sqrt(nf10[XX]*nf10[XX]+nf10[YY]*nf10[YY]+nf10[ZZ]*nf10[ZZ]) ;
   FacesAreas[11] = sqrt(nf11[XX]*nf11[XX]+nf11[YY]*nf11[YY]+nf11[ZZ]*nf11[ZZ]) ;
   FacesAreas[12] = sqrt(nf12[XX]*nf12[XX]+nf12[YY]*nf12[YY]+nf12[ZZ]*nf12[ZZ]) ;
   FacesAreas[13] = sqrt(nf13[XX]*nf13[XX]+nf13[YY]*nf13[YY]+nf13[ZZ]*nf13[ZZ]) ;
   FacesAreas[14] = sqrt(nf14[XX]*nf14[XX]+nf14[YY]*nf14[YY]+nf14[ZZ]*nf14[ZZ]) ;
   FacesAreas[15] = sqrt(nf15[XX]*nf15[XX]+nf15[YY]*nf15[YY]+nf15[ZZ]*nf15[ZZ]) ;
   FacesAreas[16] = sqrt(nf16[XX]*nf16[XX]+nf16[YY]*nf16[YY]+nf16[ZZ]*nf16[ZZ]) ;
   FacesAreas[17] = sqrt(nf17[XX]*nf17[XX]+nf17[YY]*nf17[YY]+nf17[ZZ]*nf17[ZZ]) ;
   FacesAreas[18] = sqrt(nf18[XX]*nf18[XX]+nf18[YY]*nf18[YY]+nf18[ZZ]*nf18[ZZ]) ;
   FacesAreas[19] = sqrt(nf19[XX]*nf19[XX]+nf19[YY]*nf19[YY]+nf19[ZZ]*nf19[ZZ]) ;
   FacesAreas[20] = sqrt(nf20[XX]*nf20[XX]+nf20[YY]*nf20[YY]+nf20[ZZ]*nf20[ZZ]) ;
   FacesAreas[21] =sqrt(nf21[XX]*nf21[XX]+nf21[YY]*nf21[YY]+nf21[ZZ]*nf21[ZZ]) ;



//.......for Fluxes file...........................................................
Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
//.................................................................................


  (*m_phisubT[0]) = 0.0;
  (*m_phisubT[1]) = 0.0;
  (*m_phisubT[2]) = 0.0;
  (*m_phisubT[3]) = 0.0;
  (*m_phisubT[4]) = 0.0;
  (*m_phisubT[5]) = 0.0;
  (*m_phisubT[6]) = 0.0;



  /*****         tetra 0:          *****/

  
//in substates[] I store the states of the sub-elements
  substates[0] = states[1];
  substates[1] = states[5];
  substates[2] = states[4];
  substates[3] = states[7];

  matrix_node_norms (0,XX) = nf16[XX];
  matrix_node_norms (0,YY) = nf16[YY];
  matrix_node_norms (0,ZZ) = nf16[ZZ];
 
  matrix_node_norms (1,XX) = -nf1[XX];
  matrix_node_norms (1,YY) = -nf1[YY];
  matrix_node_norms (1,ZZ) = -nf1[ZZ];

  matrix_node_norms (2,XX) = -nf14[XX];
  matrix_node_norms (2,YY) = -nf14[YY];
  matrix_node_norms (2,ZZ) = -nf14[ZZ];

  matrix_node_norms (3,XX) = -nf10[XX];
  matrix_node_norms (3,YY) = -nf10[YY];
  matrix_node_norms (3,ZZ) = -nf10[ZZ];
 
  vector_node_areas[0] = FacesAreas[16];
  vector_node_areas[1] = FacesAreas[1];
  vector_node_areas[2] = FacesAreas[14];
  vector_node_areas[3] = FacesAreas[10];

     // Finally the number of the tetrahedron is not used
     CFuint nb_tetra = 0;
    computeHOCurvedFluctuation(nb_tetra);

  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;


//...............................................................
  // The kis and the distribution are computed in this class
   computeK(substates,m_k1Plus);

//...............................................................
 


  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &distdata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
  distributeN(m_k1Plus,m_phiN1);


  computeBlendingCoeff(false, m_theta1);
  
//..............................................................




   /*****         tetra 1:          *****/


  substates[0] = states[5];
  substates[1] = states[2];
  substates[2] = states[6];
  substates[3] = states[8];
 
  
  matrix_node_norms (0,XX) = -nf6[XX];
  matrix_node_norms (0,YY) = -nf6[YY];
  matrix_node_norms (0,ZZ) = -nf6[ZZ];

  matrix_node_norms (1,XX) = nf18[XX];
  matrix_node_norms (1,YY) = nf18[YY];
  matrix_node_norms (1,ZZ) = nf18[ZZ];

  matrix_node_norms (2,XX) = -nf13[XX];
  matrix_node_norms (2,YY) = -nf13[YY];
  matrix_node_norms (2,ZZ) = -nf13[ZZ];

  matrix_node_norms (3,XX) = -nf9[XX];
  matrix_node_norms (3,YY) = -nf9[YY];
  matrix_node_norms (3,ZZ) = -nf9[ZZ];

  vector_node_areas[0] = FacesAreas[6];
  vector_node_areas[1] = FacesAreas[18];
  vector_node_areas[2] = FacesAreas[13];
  vector_node_areas[3] = FacesAreas[9];


  nb_tetra = 1;
    computeHOCurvedFluctuation(nb_tetra);
  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

  //...............................................................
  computeK(substates,m_k2Plus);
  //...............................................................
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);
  // transform fluxes of subelement to distribution variables
  distributeN(m_k2Plus,m_phiN2);
  computeBlendingCoeff(false, m_theta2);


   /*****         tetra 2:           *****/


  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[0];
  substates[3] = states[7];
  

  matrix_node_norms (0,XX) = -nf19[XX];
  matrix_node_norms (0,YY) = -nf19[YY];
  matrix_node_norms (0,ZZ) = -nf19[ZZ];

  matrix_node_norms (1,XX) = -nf0[XX];
  matrix_node_norms (1,YY) = -nf0[YY];
  matrix_node_norms (1,ZZ) = -nf0[ZZ];

  matrix_node_norms (2,XX) = -nf16[XX];
  matrix_node_norms (2,YY) = -nf16[YY];
  matrix_node_norms (2,ZZ) = -nf16[ZZ];

  matrix_node_norms (3,XX) = -nf11[XX];
  matrix_node_norms (3,YY) = -nf11[YY];
  matrix_node_norms (3,ZZ) = -nf11[ZZ];


  vector_node_areas[0] = FacesAreas[19];
  vector_node_areas[1] = FacesAreas[0];
  vector_node_areas[2] = FacesAreas[16];
  vector_node_areas[3] = FacesAreas[11];

  nb_tetra = 2;
  computeHOCurvedFluctuation(nb_tetra);
  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

 computeK(substates,m_k3Plus);
//...............................................................



  // transform fluxes of subelement to distribution variables
 *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

   distributeN(m_k3Plus,m_phiN3);
   computeBlendingCoeff(false, m_theta3);


//..............................................................


   /*****         tetra 3:            *****/


  substates[0] = states[0];
  substates[1] = states[5];
  substates[2] = states[6];
  substates[3] = states[8];

  matrix_node_norms (0,XX) = -nf18[XX];
  matrix_node_norms (0,YY) = -nf18[YY];
  matrix_node_norms (0,ZZ) = -nf18[ZZ];

  matrix_node_norms (1,XX) = -nf7[XX];
  matrix_node_norms (1,YY) = -nf7[YY];
  matrix_node_norms (1,ZZ) = -nf7[ZZ];

  matrix_node_norms (2,XX) = -nf21[XX];
  matrix_node_norms (2,YY) = -nf21[YY];
  matrix_node_norms (2,ZZ) = -nf21[ZZ];

  matrix_node_norms (3,XX) = -nf8[XX];
  matrix_node_norms (3,YY) = -nf8[YY];
  matrix_node_norms (3,ZZ) = -nf8[ZZ];


  vector_node_areas[0] = FacesAreas[18];
  vector_node_areas[1] = FacesAreas[7];
  vector_node_areas[2] = FacesAreas[21];
  vector_node_areas[3] = FacesAreas[8];

  nb_tetra = 3;
    computeHOCurvedFluctuation(nb_tetra);
  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

//...............................................................
  computeK(substates,m_k4Plus);
//...............................................................


  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
   distributeN(m_k4Plus,m_phiN4);
   computeBlendingCoeff(false, m_theta4);


//..............................................................


   /*****         tetra 4:          *****/



  substates[0] = states[5];
  substates[1] = states[7];
  substates[2] = states[8];
  substates[3] = states[0];
  
  matrix_node_norms (0,XX) =  nf20[XX];
  matrix_node_norms (0,YY) =  nf20[YY];
  matrix_node_norms (0,ZZ) =  nf20[ZZ];

  matrix_node_norms (1,XX) =  nf21[XX];
  matrix_node_norms (1,YY) =  nf21[YY];
  matrix_node_norms (1,ZZ) =  nf21[ZZ];

  matrix_node_norms (2,XX) =  nf19[XX];
  matrix_node_norms (2,YY) =  nf19[YY];
  matrix_node_norms (2,ZZ) =  nf19[ZZ];

  matrix_node_norms (3,XX) = -nf15[XX];
  matrix_node_norms (3,YY) = -nf15[YY];
  matrix_node_norms (3,ZZ) = -nf15[ZZ];


  vector_node_areas[0] =  FacesAreas[20];
  vector_node_areas[1] =  FacesAreas[21];
  vector_node_areas[2] =  FacesAreas[19];
  vector_node_areas[3] =  FacesAreas[15];

nb_tetra = 4;
    computeHOCurvedFluctuation(nb_tetra);
  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

//...............................................................
  computeK(substates,m_k5Plus);
//...............................................................

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[4]);
  distributeN(m_k5Plus,m_phiN5);
  computeBlendingCoeff(true, m_theta5);

 
   /*****         tetra 5:            *****/


  substates[0] = states[7];
  substates[1] = states[8];
  substates[2] = states[9];
  substates[3] = states[3];

  matrix_node_norms (0,XX) = -nf5[XX];
  matrix_node_norms (0,YY) = -nf5[YY];
  matrix_node_norms (0,ZZ) = -nf5[ZZ];

  matrix_node_norms (1,XX) = -nf2[XX];
  matrix_node_norms (1,YY) = -nf2[YY];
  matrix_node_norms (1,ZZ) = -nf2[ZZ];

  matrix_node_norms (2,XX) = -nf12[XX];
  matrix_node_norms (2,YY) = -nf12[YY];
  matrix_node_norms (2,ZZ) = -nf12[ZZ];

  matrix_node_norms (3,XX) = nf17[XX];
  matrix_node_norms (3,YY) = nf17[YY];
  matrix_node_norms (3,ZZ) = nf17[ZZ];


  vector_node_areas[0] = FacesAreas[5];
  vector_node_areas[1] = FacesAreas[2];
  vector_node_areas[2] = FacesAreas[12];
  vector_node_areas[3] = FacesAreas[17];
 
nb_tetra = 5;
    computeHOCurvedFluctuation(nb_tetra);
  // NOTE: compute the residual and the upwind parameters k in this cell
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

  computeK(substates,m_k6Plus);


  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[5]);
  
  distributeN(m_k6Plus,m_phiN6);
  
  computeBlendingCoeff(false, m_theta6);

  
   /*****         tetra 6:          *****/


  substates[0] = states[0];
  substates[1] = states[8];
  substates[2] = states[9];
  substates[3] = states[7];



  matrix_node_norms (0,XX) = -nf17[XX];
  matrix_node_norms (0,YY) = -nf17[YY];
  matrix_node_norms (0,ZZ) = -nf17[ZZ];

  matrix_node_norms (1,XX) = -nf3[XX];
  matrix_node_norms (1,YY) = -nf3[YY];
  matrix_node_norms (1,ZZ) = -nf3[ZZ];

  matrix_node_norms (2,XX) = -nf20[XX];
  matrix_node_norms (2,YY) = -nf20[YY];
  matrix_node_norms (2,ZZ) = -nf20[ZZ];

  matrix_node_norms (3,XX) = -nf4[XX];
  matrix_node_norms (3,YY) = -nf4[YY];
  matrix_node_norms (3,ZZ) = -nf4[ZZ];


  vector_node_areas[0] = FacesAreas[17];
  vector_node_areas[1] = FacesAreas[3];
  vector_node_areas[2] = FacesAreas[20];
  vector_node_areas[3] = FacesAreas[4];

  nb_tetra = 6;
  computeHOCurvedFluctuation(nb_tetra);
  distdata.tStates = computeConsistentStates(&substates);
  distdata.subStates = &substates;

  computeK(substates,m_k7Plus);

  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[6]);
	  
  distributeN(m_k7Plus,m_phiN7);

  computeBlendingCoeff(false, m_theta7);

 // Then we take the max of the thetas
 // this is outside because is used for the artificial viscosity
 m_theta = max(max(max(max( max( max(m_theta1, m_theta2), m_theta3), m_theta4),m_theta5),m_theta6),m_theta7);
 if (m_use_max_theta)
    { 
       if (m_theta > 0.6) {
           m_theta1 = max(m_min_theta, m_theta);
           m_theta2 = m_theta1;
           m_theta3 = m_theta1;
           m_theta4 = m_theta1;
           m_theta5 = m_theta1;
           m_theta6 = m_theta1;
           m_theta7 = m_theta1;
      }
      else {
       m_theta1 = max(m_min_theta, m_theta1);
       m_theta2 = max(m_min_theta, m_theta2);
       m_theta3 = max(m_min_theta, m_theta3);
       m_theta4 = max(m_min_theta, m_theta4);
       m_theta5 = max(m_min_theta, m_theta5);
       m_theta6 = max(m_min_theta, m_theta6);
       m_theta7 = max(m_min_theta, m_theta7);

      }


     }
 else {
       m_theta1 = max(m_min_theta, m_theta1);
       m_theta2 = max(m_min_theta, m_theta2);
       m_theta3 = max(m_min_theta, m_theta3);
       m_theta4 = max(m_min_theta, m_theta4);
       m_theta5 = max(m_min_theta, m_theta5);
       m_theta6 = max(m_min_theta, m_theta6);
       m_theta7 = max(m_min_theta, m_theta7);

     }

 if (m_store_thetas && !distdata.isPerturb)
    {
             Framework::DataHandle< CFreal > thetas = socket_thetas.getDataHandle();
             for (CFuint iEq = 0; iEq < 5; ++iEq)
             {
                thetas( distdata.cellID*7 + 0, iEq, 5) = m_theta1;
                thetas( distdata.cellID*7 + 1, iEq, 5) = m_theta2;
                thetas( distdata.cellID*7 + 2, iEq, 5) = m_theta3;
                thetas( distdata.cellID*7 + 3, iEq, 5) = m_theta4;
                thetas( distdata.cellID*7 + 4, iEq, 5) = m_theta5;
                thetas( distdata.cellID*7 + 5, iEq, 5) = m_theta6;
                thetas( distdata.cellID*7 + 6, iEq, 5) = m_theta7;              


          }
      }

     /*****         tetra 0:          *****/
     //in substates[] I store the states of the sub-elements
     substates[0] = states[1];
     substates[1] = states[5];
     substates[2] = states[4];
     substates[3] = states[7];

     // NOTE: compute the residual and the upwind parameters k in this cell
     distdata.tStates = computeConsistentStates(&substates);
     distdata.subStates = &substates;

     SafePtr<vector<RealMatrix> >& currBetaMatrix = distdata.currBetaMat;
     currBetaMatrix = &distdata.betaMats[0];
     // transform fluxes of subelement to distribution variables
     *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);

     distributeLDA(m_k1Plus,m_phi);

      residual[1] += m_theta1*m_phiN1[0] + (1.0 - m_theta1)*m_phi[0];
      residual[5] += m_theta1*m_phiN1[1] + (1.0 - m_theta1)*m_phi[1];
      residual[4] += m_theta1*m_phiN1[2] + (1.0 - m_theta1)*m_phi[2];
      residual[7] += m_theta1*m_phiN1[3] + (1.0 - m_theta1)*m_phi[3];



      /*****         tetra 1:          *****/
      
    
      substates[0] = states[5];
      substates[1] = states[2];
      substates[2] = states[6];
      substates[3] = states[8];
 


      nb_tetra = 1;
  
      // NOTE: compute the residual and the upwind parameters k in this cell
      distdata.tStates = computeConsistentStates(&substates);
      distdata.subStates = &substates;
 
     //...............................................................
     *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);
     currBetaMatrix = &distdata.betaMats[1];
     // transform fluxes of subelement to distribution variables
     distributeLDA(m_k2Plus,m_phi);

      residual[5] += m_theta2*m_phiN2[0] + (1.0 - m_theta2)*m_phi[0];
      residual[2] += m_theta2*m_phiN2[1] + (1.0 - m_theta2)*m_phi[1];
      residual[6] += m_theta2*m_phiN2[2] + (1.0 - m_theta2)*m_phi[2];
      residual[8] += m_theta2*m_phiN2[3] + (1.0 - m_theta2)*m_phi[3];

      /*****         tetra 2:           *****/


      substates[0] = states[4];
      substates[1] = states[5];
      substates[2] = states[0];
      substates[3] = states[7];

      distdata.tStates = computeConsistentStates(&substates);
      distdata.subStates = &substates;
  
      // transform fluxes of subelement to distribution variables
      *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);
       currBetaMatrix = &distdata.betaMats[2];
  
       distributeLDA(m_k3Plus,m_phi);

      residual[4] += m_theta3*m_phiN3[0] + (1.0 - m_theta3)*m_phi[0];
      residual[5] += m_theta3*m_phiN3[1] + (1.0 - m_theta3)*m_phi[1];
      residual[0] += m_theta3*m_phiN3[2] + (1.0 - m_theta3)*m_phi[2];
      residual[7] += m_theta3*m_phiN3[3] + (1.0 - m_theta3)*m_phi[3];

      /*****         tetra 3:            *****/


      substates[0] = states[0];
      substates[1] = states[5];
      substates[2] = states[6];
      substates[3] = states[8];

      distdata.tStates = computeConsistentStates(&substates);
      distdata.subStates = &substates;
  
      // transform fluxes of subelement to distribution variables
      *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
       currBetaMatrix = &distdata.betaMats[3];
       distributeLDA(m_k4Plus,m_phi);
  
      residual[0] += m_theta4*m_phiN4[0] + (1.0 - m_theta4)*m_phi[0];
      residual[5] += m_theta4*m_phiN4[1] + (1.0 - m_theta4)*m_phi[1];
      residual[6] += m_theta4*m_phiN4[2] + (1.0 - m_theta4)*m_phi[2];
      residual[8] += m_theta4*m_phiN4[3] + (1.0 - m_theta4)*m_phi[3];

      /*****         tetra 4:          *****/



      substates[0] = states[5];
      substates[1] = states[7];
      substates[2] = states[8];
      substates[3] = states[0];

      // NOTE: compute the residual and the upwind parameters k in this cell
     distdata.tStates = computeConsistentStates(&substates);
     distdata.subStates = &substates;
  
     *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[4]);
     currBetaMatrix = &distdata.betaMats[4];
     distributeLDA(m_k5Plus,m_phi);


      residual[5] += m_theta5*m_phiN5[0] + (1.0 - m_theta5)*m_phi[0];
      residual[7] += m_theta5*m_phiN5[1] + (1.0 - m_theta5)*m_phi[1];
      residual[8] += m_theta5*m_phiN5[2] + (1.0 - m_theta5)*m_phi[2];
      residual[0] += m_theta5*m_phiN5[3] + (1.0 - m_theta5)*m_phi[3];


      /*****         tetra 5:            *****/


      substates[0] = states[7];
      substates[1] = states[8];
      substates[2] = states[9];
      substates[3] = states[3];

      // NOTE: compute the residual and the upwind parameters k in this cell
      distdata.tStates = computeConsistentStates(&substates);
      distdata.subStates = &substates;
  
     // transform fluxes of subelement to distribution variables
     *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[5]);
     currBetaMatrix = &distdata.betaMats[5];
 
     distributeLDA(m_k6Plus,m_phi);

      residual[7] += m_theta6*m_phiN6[0] + (1.0 - m_theta6)*m_phi[0];
      residual[8] += m_theta6*m_phiN6[1] + (1.0 - m_theta6)*m_phi[1];
      residual[9] += m_theta6*m_phiN6[2] + (1.0 - m_theta6)*m_phi[2];
      residual[3] += m_theta6*m_phiN6[3] + (1.0 - m_theta6)*m_phi[3];


       /*****         tetra 6:          *****/


       substates[0] = states[0];
       substates[1] = states[8];
       substates[2] = states[9];
       substates[3] = states[7];

       distdata.tStates = computeConsistentStates(&substates);
       distdata.subStates = &substates;


       *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[6]);
       currBetaMatrix = &distdata.betaMats[6];

       distributeLDA(m_k7Plus,m_phi);


      residual[0] += m_theta7*m_phiN7[0] + (1.0 - m_theta7)*m_phi[0];
      residual[8] += m_theta7*m_phiN7[1] + (1.0 - m_theta7)*m_phi[1];
      residual[9] += m_theta7*m_phiN7[2] + (1.0 - m_theta7)*m_phi[2];
      residual[7] += m_theta7*m_phiN7[3] + (1.0 - m_theta7)*m_phi[3];

 
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HOCRD3D_BT_SysSplitStrategy::computeHOCurvedFluctuation(CFuint &i1)
{

  vector<State*>& states = *getMethodData().getDistributionData().states;
 DistributionData& ddata = getMethodData().getDistributionData();
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[CellID]);

   // The normal of the big elements are used to compute the
   //  shape functions

  const CFreal nx0 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx1 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(2,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(3,XX);

   const CFreal ny0 = cellnormals.getNodalNormComp(0,YY);
   const CFreal ny1 = cellnormals.getNodalNormComp(1,YY);
   const CFreal ny2 = cellnormals.getNodalNormComp(2,YY);
   const CFreal ny3 = cellnormals.getNodalNormComp(3,YY);
   
   const CFreal nz0 = cellnormals.getNodalNormComp(0,ZZ);
   const CFreal nz1 = cellnormals.getNodalNormComp(1,ZZ);
   const CFreal nz2 = cellnormals.getNodalNormComp(2,ZZ);
   const CFreal nz3 = cellnormals.getNodalNormComp(3,ZZ);




  CFreal inv_sixvolume = 1.0/(6.0*ddata.cell->computeVolume());


    State& state0 = *(states[0]);
   State& state1 = *(states[1]);
   State& state2 = *(states[2]);
   State& state3 = *(states[3]);
   State& state4 = *(states[4]);
   State& state5 = *(states[5]);
   State& state6 = *(states[6]);
   State& state7 = *(states[7]);
   State& state8 = *(states[8]);
   State& state9 = *(states[9]);

   const CFreal x0 = states[0]->getCoordinates()[XX];
   const CFreal x1 = states[1]->getCoordinates()[XX];
   const CFreal x2 = states[2]->getCoordinates()[XX];
   const CFreal x3 = states[3]->getCoordinates()[XX];

   const CFreal y0 = states[0]->getCoordinates()[YY];
   const CFreal y1 = states[1]->getCoordinates()[YY];
   const CFreal y2 = states[2]->getCoordinates()[YY];
   const CFreal y3 = states[3]->getCoordinates()[YY];

   const CFreal z0 = states[0]->getCoordinates()[ZZ];
   const CFreal z1 = states[1]->getCoordinates()[ZZ];
   const CFreal z2 = states[2]->getCoordinates()[ZZ];
   const CFreal z3 = states[3]->getCoordinates()[ZZ];


 // Computation  of the flux of each faces of the cell
//and then, computation of the subcell residual..
  for (CFuint i = 0; i < 4; ++i){
      CFuint iFace = subelemtable(i1,i);
      Facenormal0[XX] = FacesMat(iFace,XX);
      Facenormal0[YY] = FacesMat(iFace,YY);
      Facenormal0[ZZ] = FacesMat(iFace,ZZ);

         for (CFuint iQd = 0; iQd < nbQdPts; ++iQd)
	   {	

//       // Values of shape functions in reference space:
	     const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
	     const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
	     const Node& node2 = states[subfacetable(iFace,2)]->getCoordinates();
	     
	     const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX]+qd2[iQd]*node2[XX];
	     const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY]+qd2[iQd]*node2[YY];
	     const CFreal z = qd0[iQd] * node0[ZZ] + qd1[iQd] * node1[ZZ]+qd2[iQd]*node2[ZZ];

	     CFreal L1 = 1.0 + ( ( x - x0 )*nx0 + ( y - y0 )*ny0 + (z - z0)*nz0 )*inv_sixvolume ;
	     CFreal L2 = 1.0 + ( ( x - x1 )*nx1 + ( y - y1 )*ny1 + (z - z1)*nz1 )*inv_sixvolume ;
	     CFreal L3 = 1.0 + ( ( x - x2 )*nx2 + ( y - y2 )*ny2 + (z - z2)*nz2 )*inv_sixvolume ;
	     CFreal L4 = 1.0 + ( ( x - x3 )*nx3 + ( y - y3 )*ny3 + (z - z3)*nz3 )*inv_sixvolume ;

	      const CFreal SF0 = - L1 * (1.0 - 2.0 * L1);
	      const CFreal SF1 = - L2 * (1.0 - 2.0 * L2);
	      const CFreal SF2 = - L3 * (1.0 - 2.0 * L3);
	      const CFreal SF3 = - L4 * (1.0 - 2.0 * L4);
	      const CFreal SF4 = 4.0 * L2 * L1;
	     const CFreal SF5 = 4.0 * L2 * L3;
	      const CFreal SF6 = 4.0 * L3 * L1;
	      const CFreal SF7 = 4.0 * L2 * L4;
	      const CFreal SF8 = 4.0 * L3 * L4;
	      const CFreal SF9 = 4.0 * L4 * L1;

	      (*qdstates[iQd]) = SF0 * state0 + SF1 * state1 + SF2 * state2 + SF3 * state3 + SF4 * state4 +
	        SF5 * state5 + SF6 * state6 + SF7 * state7 + SF8 * state8 + SF9 * state9;


	     
	   }
	     
	   computeStatesData(4, m_updateVar, qdstates, m_pdata, m_qdExtraVars); // three quadrature points per face
	   faceflux[iFace] = 0.;
	   for (CFuint iQd = 0; iQd < 4; ++iQd)
	     {
	       
	       faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()(m_pdata[iQd],Facenormal0);
	     }
  }


 
    RealVector& phi = (*m_phisubT[i1]);
    phi = 0.;
    for (CFuint jCol = 0; jCol < subelemtable.nbCols(); ++jCol)
    {

      phi += subelemfacedir(i1,jCol) * faceflux[subelemtable(i1,jCol)];
    }
  
}
//////////////////////////////////////////////////////////////////////////////////////
void HOCRD3D_BT_SysSplitStrategy::computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus){
  
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFreal m_invDim = 1.0/3.0;
  // The transformation of the normal is needed if the coordinate system is rotated
  for (CFuint iState = 0; iState < 4; ++iState) {
    // The transformation of the normal is needed if the coordinate system
    // is rotated. The normal is adimensionalized, so there is need to multiply
    // by the nodeArea when computing the k parameter
  
    for (CFuint iDim = 0; iDim < 3; ++iDim) {
      
        m_adimNormal[iDim] =  matrix_node_norms(iState, iDim);
     
    }		
    m_adimNormal *= 1. / vector_node_areas[iState];
 
     
     getMethodData().getDistribVar()->splitJacobian(*m_kPlus[iState],
     					       *m_kMin[iState],
     					       *m_eValues[iState],
     					       m_adimNormal);
CFuint nbEqs =  (*m_eValues[iState]).size();
  // First we check if we are at a stagnation point which coorespond 
   // to some eigen value that are null
   bool istagnpoint = false;
   for (CFuint iEq = 0; iEq < nbEqs; ++iEq){
     if ((*m_eValues[iState])[iEq] < 1.0e-8 )
      istagnpoint = true;
   }

   
     CFreal m_nodeArea = vector_node_areas[iState];
     if (!istagnpoint){    
     *m_kPlus[iState] *= m_invDim * m_nodeArea;
     *m_kMin[iState]  *= m_invDim * m_nodeArea;
      }
     else {
     getMethodData().getDistribVar()->computeEigenValuesVectors(m_rightEv, m_leftEv, *m_eValues[iState],m_adimNormal);


     // If we are at a stagnation point we add an epsilon to all the eigenvalue to be sure that
     // the K will keep invertible
     for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            (*m_eValues[iState])[iEq] += 1.0e-8;
              m_eValuesP[iEq] = m_invDim * m_nodeArea*max(0.,(*m_eValues[iState])[iEq]);
             m_eValuesM[iEq] = m_invDim * m_nodeArea*min(0.,(*m_eValues[iState])[iEq]);
     }
     
     // compute jacobian + and -
     *m_kPlus[iState] = m_rightEv*(m_eValuesP*m_leftEv);
     *m_kMin[iState] = m_rightEv*(m_eValuesM*m_leftEv);
     }
    
     
    if (!getMethodData().getDistributionData().isPerturb) {
       const CFreal maxEigenValue = std::max(0.0, m_eValues[iState]->max());
      
       updateCoeff[states[iState]->getLocalID()] += m_invDim*m_nodeArea*maxEigenValue;
      
     }
   }



}



//////////////////////////////////////////////////////////////////////////////
void HOCRD3D_BT_SysSplitStrategy::distributeLDA(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiLDA){

 const RealVector& phiT = getMethodData().getDistributionData().phi;
  DistributionData& ddata = getMethodData().getDistributionData();
   m_sumKplus = *m_kPlus[0];

  for (CFuint iState = 1; iState < 4; ++iState) {
     m_sumKplus  += *m_kPlus[iState];
  }


   m_inverter->invert(m_sumKplus, m_invK);

   m_uTemp = m_invK*phiT;




  for (CFuint iState = 0; iState < 4; ++iState) {
     phiLDA[iState] = (*m_kPlus[iState])*m_uTemp;
    
      (*ddata.currBetaMat)[iState] = (*m_kPlus[iState]) * m_invK;
   
    }


}
//////////////////////////////////////////////////////////////////////////////
void HOCRD3D_BT_SysSplitStrategy::distributeN(std::vector<RealMatrix*> & m_kPlus,vector<RealVector> & phiN){

 const RealVector& phiT = getMethodData().getDistributionData().phi;
   DistributionData& ddata = getMethodData().getDistributionData();
const vector<State*>& tStates = *getMethodData().getDistributionData().tStates;
   m_sumKplus = *m_kPlus[0];
   m_sumKplusU = *m_kPlus[0]*(*tStates[0]);

  for (CFuint iState = 1; iState < 4; ++iState) {
     m_sumKplus  += *m_kPlus[iState];
     m_sumKplusU += *m_kPlus[iState]*(*tStates[iState]);

  }


   m_inverter->invert(m_sumKplus, m_invK);


   m_uTemp = m_invK*(m_sumKplusU -phiT);


  for (CFuint iState = 0; iState < 4; ++iState) {
    phiN[iState] = *m_kPlus[iState]*(*tStates[iState]-m_uTemp);

    }


}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRD3D_BT_SysSplitStrategy::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
   result.push_back(&socket_thetas);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
void HOCRD3D_BT_SysSplitStrategy::computeBlendingCoeff(bool tetra, CFreal & result)
{
 vector<State*> states = substates;

 CFreal vol = getMethodData().getDistributionData().cell->computeVolume();

 if (!tetra) vol = vol/8.0;
 else vol = vol/4.0;

 m_updateVar   = getMethodData().getUpdateVar();
 const RealVector& lData =  _cterm->getPhysicalData();

 const CFuint nbStates = states.size();
 const CFuint dim = PhysicalModelStack::getActive()->getDim();

 _grad = 0.0;
 for (CFuint i = 0; i < nbStates; ++i) {
      m_updateVar->computePhysicalData(*states[i], _pData);
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
           _grad[iDim] += _pData[_varID]*matrix_node_norms(i,iDim);
                                                 }
                                       }
  const CFreal h = 2.0*std::pow(3.0*vol/(4.0*MathTools::MathConsts::CFrealPi()),0.33333);

  CFreal sc;
  if (_sh_detector == "Jirka")
        {
          sc = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY]+_grad[ZZ]*lData[EulerTerm::VZ];
          sc *= _length/(vol*dim*_deltaP*_speed);
          m_sc = sc;
         }
   else if (_sh_detector == "Anton")
            {
              sc = _grad[XX]*lData[EulerTerm::VX] + _grad[YY]*lData[EulerTerm::VY]+_grad[ZZ]*lData[EulerTerm::VZ];
              sc *= std::pow(    (std::pow(lData[EulerTerm::VX],2) + std::pow(lData[EulerTerm::VY],2)+std::pow(lData[EulerTerm::VZ],2)), 0.5   );
              sc *= _length/(vol*dim*_deltaP*_speed * _speed);
              m_sc = sc;
             }

    //second test
    result = min(1.0,(max(0.0, m_sc))*(max(0.0, m_sc))*h);


}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
