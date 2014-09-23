#include "FluctSplit/HOCRDP3_SplitStrategy.hh"

#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/FluctSplitHO.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<HOCRDP3_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
hocrp3dFluctSplitStrategyProvider("HOCRDP3");

//////////////////////////////////////////////////////////////////////////////

HOCRDP3_SplitStrategy::HOCRDP3_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals()
{
}

//////////////////////////////////////////////////////////////////////////////

HOCRDP3_SplitStrategy::~HOCRDP3_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_SplitStrategy::unsetup()
{
  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_SplitStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  if (getMethodData().isMultipleSplitter())
    throw BadValueException (FromHere(),"Cannot use HOCRD with multiple splitters");

  m_splitter = getMethodData().getSplitter();

  m_solutionVar = getMethodData().getSolutionVar();
  m_updateVar   = getMethodData().getUpdateVar();

  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {

    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

   // number of quadrature point used to compute the fluctuatuion
  const CFuint nbQdPts = 5;

 // physical data evaluated in the quadrature points
  m_pdata.resize(nbQdPts);
  for (CFuint  i = 0; i < nbQdPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  

  // sub elemt table
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);

  m_phisubT.resize(9); // P3 triangle has 9 sub triangles

  for (CFuint isub_elem = 0; isub_elem < m_phisubT.size(); ++ isub_elem)
    m_phisubT[isub_elem] = new RealVector(nbEqs);

  // sub elemt table : contain the faces of each sub-element
  subelemtable.resize(9,3); // 9 sub elems with 3 faces each
  subelemfacedir.resize(9,3);

  subelemtable(0,0) = 0;
  subelemtable(0,1) = 1;
  subelemtable(0,2) = 2;

  subelemtable(1,0) = 3;
  subelemtable(1,1) = 4;
  subelemtable(1,2) = 5;

  subelemtable(2,0) = 6;
  subelemtable(2,1) = 7;
  subelemtable(2,2) = 8;

  subelemtable(3,0) = 9;
  subelemtable(3,1) = 10;
  subelemtable(3,2) = 11;

  subelemtable(4,0) = 12;
  subelemtable(4,1) = 13;
  subelemtable(4,2) = 14;

  subelemtable(5,0) = 15;
  subelemtable(5,1) = 16;
  subelemtable(5,2) = 17;

  subelemtable(6,0) = 0;  // faces have negative orientation
  subelemtable(6,1) = 16;
  subelemtable(6,2) = 14;

  subelemtable(7,0) = 15; // faces have negative orientation
  subelemtable(7,1) = 4;
  subelemtable(7,2) = 11;

  subelemtable(8,0) = 12; // faces have negative orientation
  subelemtable(8,1) = 10;
  subelemtable(8,2) = 8;


// sub elemt table : contain the orientation of the faces of each sub-element
  subelemfacedir(0,0) = 1.;
  subelemfacedir(0,1) = 1.;
  subelemfacedir(0,2) = 1.;

  subelemfacedir(1,0) = 1.;
  subelemfacedir(1,1) = 1.;
  subelemfacedir(1,2) = 1.;

  subelemfacedir(2,0) = 1.;
  subelemfacedir(2,1) = 1.;
  subelemfacedir(2,2) = 1.;

  subelemfacedir(3,0) = 1.;
  subelemfacedir(3,1) = 1.;
  subelemfacedir(3,2) = 1.;

  subelemfacedir(4,0) = 1.;
  subelemfacedir(4,1) = 1.;
  subelemfacedir(4,2) = 1.;

  subelemfacedir(5,0) = 1.;
  subelemfacedir(5,1) = 1.;
  subelemfacedir(5,2) = 1.;

  subelemfacedir(6,0) = -1.;
  subelemfacedir(6,1) = -1.;
  subelemfacedir(6,2) = -1.;

  subelemfacedir(7,0) = -1.;
  subelemfacedir(7,1) = -1.;
  subelemfacedir(7,2) = -1.;

  subelemfacedir(8,0) = -1.;
  subelemfacedir(8,1) = -1.;
  subelemfacedir(8,2) = -1.;

  // sub face table : contain the node that contain each face
  subfacetable.resize(18,2); // 18 sub faces with 2 states each

  subfacetable(0,0) = 3;
  subfacetable(0,1) = 8;

  subfacetable(1,0) = 8;
  subfacetable(1,1) = 0;

  subfacetable(2,0) = 0;
  subfacetable(2,1) = 3;

  subfacetable(3,0) = 1;
  subfacetable(3,1) = 5;

  subfacetable(4,0) = 5;
  subfacetable(4,1) = 4;

  subfacetable(5,0) = 4;
  subfacetable(5,1) = 1;

  subfacetable(6,0) = 6;
  subfacetable(6,1) = 2;

  subfacetable(7,0) = 2;
  subfacetable(7,1) = 7;

  subfacetable(8,0) = 7;
  subfacetable(8,1) = 6;

  subfacetable(9,0) = 5;
  subfacetable(9,1) = 6;

  subfacetable(10,0) = 6;
  subfacetable(10,1) = 9;

  subfacetable(11,0) = 9;
  subfacetable(11,1) = 5;

  subfacetable(12,0) = 9;
  subfacetable(12,1) = 7;

  subfacetable(13,0) = 7;
  subfacetable(13,1) = 8;

  subfacetable(14,0) = 8;
  subfacetable(14,1) = 9;

  subfacetable(15,0) = 4;
  subfacetable(15,1) = 9;

  subfacetable(16,0) = 9;
  subfacetable(16,1) = 3;

  subfacetable(17,0) = 3;
  subfacetable(17,1) = 4;

  faceflux.resize(subfacetable.nbRows()); // one flux per sub face
  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(nbEqs);

  const CFreal s1  = (1.0/3.0)*sqrt(5.0-2.0*sqrt(10.0/7.0));
  const CFreal a01 = ( 1.0 - s1 )*0.5;
  const CFreal a11 = ( 1.0 + s1 )*0.5;

  const CFreal s2  = (1.0/3.0)*sqrt(5.0+2.0*sqrt(10.0/7.0));
  const CFreal a02 = ( 1.0 - s2 )*0.5;
  const CFreal a12 = ( 1.0 + s2 )*0.5;

  qd0.resize(nbQdPts); // quadrature points per face
  qd1.resize(nbQdPts); // quadrature points per face

  qd0[0] = a01;  qd1[0] = a11;
  qd0[1] = a11;  qd1[1] = a01;
  qd0[2] = a02;  qd1[2] = a12;
  qd0[3] = a12;  qd1[3] = a02;
  qd0[4] = .5;   qd1[4] = .5;

  wqd.resize(nbQdPts); // 3 quadrature points per face
  wqd[0] = (322.0+13.0*sqrt(70.0))/1800.0;
  wqd[1] = (322.0+13.0*sqrt(70.0))/1800.0;
  wqd[2] = (322.0-13.0*sqrt(70.0))/1800.0;
  wqd[3] = (322.0-13.0*sqrt(70.0))/1800.0;
  wqd[4] = 128.0/450.0;

  qdstates.resize(nbQdPts); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();
  qdstates[3] = new State();
  qdstates[4] = new State();


  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(nbQdPts); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);
  qdnodes[3] = new Node();
  qdstates[3]->setSpaceCoordinates(qdnodes[3]);
  qdnodes[4] = new Node();
  qdstates[4]->setSpaceCoordinates(qdnodes[4]);

  facenormal.resize(2); // only 2D

  // resize the beta's storages in the MethodData
  const CFuint nbSubTriangles = 9;
  const CFuint nbNodesInTriangle = 3;
  DistributionData::BetaMatrices& betaMats =
    getMethodData().getDistributionData().betaMats;
  betaMats.resize(nbSubTriangles);

  for (CFuint t = 0; t < nbSubTriangles; ++t) {
    betaMats[t].resize(nbNodesInTriangle);

    for (CFuint i = 0; i < nbNodesInTriangle; ++i) {

      betaMats[t][i].resize(nbEqs, nbEqs);
    }
  }

  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = true;
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  vector<State*>& states = *getMethodData().getDistributionData().states;
  cf_assert(states.size() == 10); // P3 triangles for solution space
  cf_assert(residual.size() >= states.size()); // state residual
  DistributionData& ddata = getMethodData().getDistributionData();
  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i)
  {
    residual[i] = 0.0;
  }
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);

  cf_assert(cellnormals.nbFaces()  == 3); // triangles have 3 faces

  CFreal onethird = 1.0/3.0;

  // normals are half scale
  cellnormals.scale(onethird);
  computeHOP3Fluctuation();

  /*****         Triangle 0-3-8          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[8];

  getMethodData().getDistributionData().tStates =
    computeConsistentStates(&substates);
  
  ddata.subStates = &substates;
  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[0]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  SafePtr<RealVector> phi = &getMethodData().getDistributionData().phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
 
  

  SafePtr<vector<RealMatrix> >& currBetaMatrix =
    getMethodData().getDistributionData().currBetaMat;


  // cout <<"getMethodData().getDistributionData().betaMats.size() = " <<
  //  getMethodData().getDistributionData().betaMats.size()  << endl;

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[0];

  // cout << "*currBetaMatrix0 = " << (*currBetaMatrix)[0] << endl;
  // cout << "*currBetaMatrix1 = " << (*currBetaMatrix)[1] << endl;
  // cout << "*currBetaMatrix2 = " << (*currBetaMatrix)[2] << endl;

  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[8] += subresidual[2];


  /*****         Triangle 4-1-5          *****/
  substates[0] = states[4];
  substates[1] = states[1];
  substates[2] = states[5];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);
  
  ddata.subStates = &substates;
  cellnormals.scale(onethird);
  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[1]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[1];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[4] += subresidual[0];
  residual[1] += subresidual[1];
  residual[5] += subresidual[2];

  /*****         Triangle 7-6-2          *****/
  substates[0] = states[7];
  substates[1] = states[6];
  substates[2] = states[2];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(onethird);
  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[2]) -= ddata.phiS;
  }
  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);

	currBetaMatrix = &getMethodData().getDistributionData().betaMats[2];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[7] += subresidual[0];
  residual[6] += subresidual[1];
  residual[2] += subresidual[2];

  /*****         Triangle 9-5-6          *****/
  substates[0] = states[9];
  substates[1] = states[5];
  substates[2] = states[6];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(onethird);
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[3]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[3];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[9] += subresidual[0];
  residual[5] += subresidual[1];
  residual[6] += subresidual[2];

 /*****         Triangle 8-9-7          *****/
  substates[0] = states[8];
  substates[1] = states[9];
  substates[2] = states[7];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(onethird);
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[4]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[4]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[4];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[8] += subresidual[0];
  residual[9] += subresidual[1];
  residual[7] += subresidual[2];


/*****         Triangle 3-4-9          *****/
  substates[0] = states[3];
  substates[1] = states[4];
  substates[2] = states[9];

   getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(onethird);
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[5]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[5]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[5];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[3] += subresidual[0];
  residual[4] += subresidual[1];
  residual[9] += subresidual[2];


// compute the residual and the upwind parameters k in these cells
  // the oriantation of the normal in this sub-element is oposite to the one of the element
  cellnormals.scale(-onethird);

/*****         Triangle 9-8-3          *****/
  substates[0] = states[9];
  substates[1] = states[8];
  substates[2] = states[3];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;

  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[6]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[6]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[6];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(-onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[9] += subresidual[0];
  residual[8] += subresidual[1];
  residual[3] += subresidual[2];

/*****         Triangle 5-9-4          *****/
  substates[0] = states[5];
  substates[1] = states[9];
  substates[2] = states[4];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(-onethird);
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[7]) -= ddata.phiS;
  }

  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[7]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[7];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(-onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[5] += subresidual[0];
  residual[9] += subresidual[1];
  residual[4] += subresidual[2];

/*****         Triangle 6-7-9          *****/
  substates[0] = states[6];
  substates[1] = states[7];
  substates[2] = states[9];

  getMethodData().getDistributionData().tStates = computeConsistentStates(&substates);

  ddata.subStates = &substates;
  cellnormals.scale(-onethird);
  m_splitter->computeK(substates,&cellnormals);
  cellnormals.unscale();
  ddata.sourceTermID = 0;

  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[8]) -= ddata.phiS;
  }
  // transform fluxes of subelement to distribution variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[8]);

  currBetaMatrix = &getMethodData().getDistributionData().betaMats[8];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.scale(-onethird);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[6] += subresidual[0];
  residual[7] += subresidual[1];
  residual[9] += subresidual[2];

  cellnormals.unscale();
  //CF_DEBUG_EXIT;
}

//////////////////////////////////////////////////////////////////////////////

void HOCRDP3_SplitStrategy::computeHOP3Fluctuation()
{

  cf_assert(qdstates.size() == 5); // only triags so three quadrature points per face

  vector<State*>& states = *getMethodData().getDistributionData().states;

  const State& state0 = *(states[0]);
  const State& state1 = *(states[1]);
  const State& state2 = *(states[2]);
  const State& state3 = *(states[3]);
  const State& state4 = *(states[4]);
  const State& state5 = *(states[5]);
  const State& state6 = *(states[6]);
  const State& state7 = *(states[7]);
  const State& state8 = *(states[8]);
  const State& state9 = *(states[9]);


  vector<Node*>& nodes = *getMethodData().getDistributionData().cell->getNodes();

  cf_assert(nodes.size()  == 3); // P1 triangles for geometry space

  const CFreal x1 = (*nodes[0])[XX];
  const CFreal x2 = (*nodes[1])[XX];
  const CFreal x3 = (*nodes[2])[XX];

  const CFreal y1 = (*nodes[0])[YY];
  const CFreal y2 = (*nodes[1])[YY];
  const CFreal y3 = (*nodes[2])[YY];

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [getMethodData().getDistributionData().cellID]);
  // we multiply by two because the normals scale has been changed in computeFluctuationAndUpdateCoeff
  // and is 0.5
  const CFreal nx1 = 3. * cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = 3. * cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = 3. * cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = 3. * cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = 3. * cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = 3. * cellnormals.getNodalNormComp(2,YY);


  const CFreal inv_volume = 1.0 / getMethodData().getDistributionData().cell->computeVolume();
// CF_DEBUG_POINT;

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {
    facenormal[XX] = cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = cellnormals.getNodalNormComp(iFace%3,YY);


    for (CFuint iQd = 0; iQd < 5; ++iQd)
    {

      const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
      const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
      const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX];
      const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY];
      const CFreal L1 = 1.0 + 0.5*( ( x - x1 )*nx1 + ( y - y1 )*ny1 ) * inv_volume ;
      const CFreal L2 = 1.0 + 0.5*( ( x - x2 )*nx2 + ( y - y2 )*ny2 ) * inv_volume ;
      const CFreal L3 = 1.0 + 0.5*( ( x - x3 )*nx3 + ( y - y3 )*ny3 ) * inv_volume ;

      const CFreal a0 = 0.5*( 3.0*L1 - 1.0 )*( 3.0*L1 - 2.0 )*L1 ;
      const CFreal a1 = 0.5*( 3.0*L2 - 1.0 )*( 3.0*L2 - 2.0 )*L2 ;
      const CFreal a2 = 0.5*( 3.0*L3 - 1.0 )*( 3.0*L3 - 2.0 )*L3 ;

      const CFreal a3 = 4.5*L1*L2*( 3.0*L1 - 1.0 ) ;
      const CFreal a4 = 4.5*L1*L2*( 3.0*L2 - 1.0 ) ;

      const CFreal a5 = 4.5*L2*L3*( 3.0*L2 - 1.0 ) ;
      const CFreal a6 = 4.5*L2*L3*( 3.0*L3 - 1.0 ) ;

      const CFreal a7 = 4.5*L3*L1*( 3.0*L3 - 1.0 ) ;
      const CFreal a8 = 4.5*L3*L1*( 3.0*L1 - 1.0 ) ;

      const CFreal a9 = 27.0*L1*L2*L3;

      (*qdstates[iQd]) = a0*state0 + a1*state1 + a2*state2 + a3*state3 +
                         a4*state4 + a5*state5 + a6*state6 + a7*state7 +
                         a8*state8 + a9*state9;

      (*qdnodes[iQd])[XX] = x;
      (*qdnodes[iQd])[YY] = y;
    }
   
    computeStatesData(5, m_updateVar, qdstates, m_pdata, m_qdExtraVars); // three quadrature points per face

    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 5; ++iQd)
    {

      faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()((*qdstates[iQd]),facenormal);
    }
  }

  cf_assert(m_phisubT.size() == subelemtable.nbRows());

 // recomposition of the fluctuation of each sub cell (taking care of the orientation)
  for (CFuint iRow = 0; iRow < subelemtable.nbRows(); ++iRow)
  {
    RealVector& phi = (*m_phisubT[iRow]);
    phi = 0.;
    for (CFuint jCol = 0; jCol < subelemtable.nbCols(); ++jCol)
    {

      phi -= subelemfacedir(iRow,jCol) * faceflux[subelemtable(iRow,jCol)];
    }
  }


}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
HOCRDP3_SplitStrategy::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////
    }//end namespace FluctSplit

}// end namespace COOLFluiD
