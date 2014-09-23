#include "FluctSplit/HOCRD_SplitStrategy.hh"

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

MethodStrategyProvider<HOCRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitHOModule>
hocrdFluctSplitStrategyProvider("HOCRD");

//////////////////////////////////////////////////////////////////////////////

HOCRD_SplitStrategy::HOCRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  m_solutionVar(CFNULL),
  m_updateVar(CFNULL),
  m_unitFaceNormals(),
  m_phiT()
{
}

//////////////////////////////////////////////////////////////////////////////

HOCRD_SplitStrategy::~HOCRD_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategy::unsetup()
{
  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
   deletePtr(m_qdExtraVars[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategy::setup()
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
  const CFuint nbQdPts = 3;

  // physical data evaluated in the quadrature points
  m_pdata.resize(nbQdPts);
  for (CFuint  i = 0; i < nbQdPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
  // sub element table: states, residual
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);

  m_phisubT.resize(4); // P2 triangle has 4 sub triangles

  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);


  // sub elemt table : contain the faces of each sub-element
  // See the manual for the numbering
  //
  subelemtable.resize(4,3); // 4 sub elems with 3 faces each
  subelemfacedir.resize(4,3); //this table will set the sign of the face residual when the assembling will be done

  subelemtable(0,0) = 0;
  subelemtable(0,1) = 1;
  subelemtable(0,2) = 2;

  subelemtable(1,0) = 3;
  subelemtable(1,1) = 4;
  subelemtable(1,2) = 5;

  subelemtable(2,0) = 6;
  subelemtable(2,1) = 7;
  subelemtable(2,2) = 8;

  subelemtable(3,0) = 0; // faces have negative orientation
  subelemtable(3,1) = 4;
  subelemtable(3,2) = 8;

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

  subelemfacedir(3,0) = -1.; // faces have negative orientation
  subelemfacedir(3,1) = -1.;
  subelemfacedir(3,2) = -1.;


  // sub face table : contain the node that contain each face
  subfacetable.resize(9,2); // 9 sub faces with 2 states each

  subfacetable(0,0) = 3;
  subfacetable(0,1) = 5;

  subfacetable(1,0) = 5;
  subfacetable(1,1) = 0;

  subfacetable(2,0) = 0;
  subfacetable(2,1) = 3;

  subfacetable(3,0) = 1;
  subfacetable(3,1) = 4;

  subfacetable(4,0) = 4;
  subfacetable(4,1) = 3;

  subfacetable(5,0) = 3;
  subfacetable(5,1) = 1;

  subfacetable(6,0) = 4;
  subfacetable(6,1) = 2;

  subfacetable(7,0) = 2;
  subfacetable(7,1) = 5;

  subfacetable(8,0) = 5;
  subfacetable(8,1) = 4;

  faceflux.resize(subfacetable.nbRows()); // one flux per sub face
  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(nbEqs);


// Definition of the quadrature rule
  const CFreal s  = std::sqrt( 0.6 );
  const CFreal a0 = ( 1.0 - s )*0.5;
  const CFreal a1 = ( 1.0 + s )*0.5;

  qd0.resize(nbQdPts); // quadrature points per face
  qd1.resize(nbQdPts); // quadrature points per face

  qd0[0] = a0;  qd1[0] = a1;
  qd0[1] = a1;  qd1[1] = a0;
  qd0[2] = .5;  qd1[2] = .5;

  wqd.resize(nbQdPts); // 3 quadrature points per face
  wqd[0] = 5.0/18.0;
  wqd[1] = 5.0/18.0;
  wqd[2] = 8.0/18.0;

  qdstates.resize(3); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();

  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(qdstates.size());
  
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  qdnodes.resize(3); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);


  facenormal.resize(2); // only 2D

  // The beta storage is used when a viscous term is also discretized
  // resize the beta's storages in the MethodData
  const CFuint nbSubTriangles = 4;
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

 m_phiT.resize(PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  vector<State*>& states = *getMethodData().getDistributionData().states;
  DistributionData& ddata = getMethodData().getDistributionData();
  cf_assert(states.size() == 6); // P2 triangles for solution space
  cf_assert(residual.size() >= states.size()); // state residual


  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i)
  {
    residual[i] = 0.0;
  }

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[ddata.cellID]);


  cf_assert(cellnormals.nbFaces()  == 3); // triangles have 3 faces

  // normals of the sub-element are half of the normals of the big element
  // 
  cellnormals.scale(0.5);
  computeHOFluctuation();

  /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  ddata.tStates = computeConsistentStates(&substates);

   ddata.subStates = &substates;

  // compute the residual and the upwind parameters k in this cell
  m_splitter->computeK(substates,&cellnormals);


  // Before distributing the residual, we point to the good beta matrix
  SafePtr<vector<RealMatrix> >& currBetaMatrix =
    ddata.currBetaMat;

  currBetaMatrix = &ddata.betaMats[0];

  cf_assert(currBetaMatrix.isNotNull());


  cellnormals.unscale();
    ddata.sourceTermID = 0;
  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[0]) -= ddata.phiS;
  }

   // transform fluxes + source term to distribute variables
  SafePtr<RealVector> phi = &ddata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
  cellnormals.scale(0.5);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[5] += subresidual[2];

  /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  ddata.tStates = computeConsistentStates(&substates);
  *ddata.subStates = substates;
  // compute the residual and the upwind parameters k in this cell
  cellnormals.scale(0.5);
  m_splitter->computeK(substates,&cellnormals);

  currBetaMatrix = &ddata.betaMats[1];
  cf_assert(currBetaMatrix.isNotNull());
  cellnormals.unscale();

    ddata.sourceTermID = 0;
  getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[1]) -= ddata.phiS;
  }

   // transform fluxes + source term to distribute variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);
  cellnormals.scale(0.5);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);



  residual[3] += subresidual[0];
  residual[1] += subresidual[1];
  residual[4] += subresidual[2];

  /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  ddata.tStates = computeConsistentStates(&substates);
  *ddata.subStates = substates;
  // compute the residual and the upwind parameters k in this cell
cellnormals.scale(0.5);
  m_splitter->computeK(substates,&cellnormals);

  currBetaMatrix = &ddata.betaMats[2];
  cf_assert(currBetaMatrix.isNotNull());
 cellnormals.unscale();

    ddata.sourceTermID = 0;

getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[2]) -= ddata.phiS;
  }

   // transform fluxes + source term to distribute variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);
  cellnormals.scale(0.5);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);


  residual[5] += subresidual[0];
  residual[4] += subresidual[1];
  residual[2] += subresidual[2];

  /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  ddata.tStates = computeConsistentStates(&substates);
  *ddata.subStates = substates;
  // compute the residual and the upwind parameters k in this cell
  // the oriantation of the normal in this sub-element is oposite to the one of the element
  cellnormals.scale(-0.5);
  m_splitter->computeK(substates,&cellnormals);

  currBetaMatrix = &ddata.betaMats[3];

  cf_assert(currBetaMatrix.isNotNull());
 cellnormals.unscale();


    ddata.sourceTermID = 0;
getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

  if (getMethodData().includeSourceInFlux()) {
    // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[3]) -= ddata.phiS;
  }

   // transform fluxes + source term to distribute variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
  cellnormals.scale(-0.5);
  m_splitter->distribute(subresidual);
  cellnormals.unscale();
  getMethodData().getSourceTermSplitter(0)->distribute(subresidual);



  residual[4] += subresidual[0];
  residual[5] += subresidual[1];
  residual[3] += subresidual[2];
}

//////////////////////////////////////////////////////////////////////////////

void HOCRD_SplitStrategy::computeHOFluctuation()
{
  DistributionData& ddata = getMethodData().getDistributionData();
  cf_assert(qdstates.size() == 3); // only triags so three quadrature points per face

  vector<State*>& states = *ddata.states;

  const State& state0 = *(states[0]);
  const State& state1 = *(states[1]);
  const State& state2 = *(states[2]);
  const State& state3 = *(states[3]);
  const State& state4 = *(states[4]);
  const State& state5 = *(states[5]);

  vector<Node*>& nodes = *ddata.cell->getNodes();

  cf_assert(nodes.size()  == 3); // P1 triangles for geometry space

  const CFreal x1 = (*nodes[0])[XX];
  const CFreal x2 = (*nodes[1])[XX];
  const CFreal x3 = (*nodes[2])[XX];

  const CFreal y1 = (*nodes[0])[YY];
  const CFreal y2 = (*nodes[1])[YY];
  const CFreal y3 = (*nodes[2])[YY];

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  // we multiply by two because the normals scale has been changed in computeFluctuationAndUpdateCoeff
  // and is 0.5
  const CFreal nx1 = 2. * cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = 2. * cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = 2. * cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = 2. * cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = 2. * cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = 2. * cellnormals.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0 / ddata.cell->computeVolume();

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {
    facenormal[XX] = cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = cellnormals.getNodalNormComp(iFace%3,YY);


    for (CFuint iQd = 0; iQd < 3; ++iQd)
    {

      const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
      const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
      const CFreal x = qd0[iQd] * node0[XX] + qd1[iQd] * node1[XX];
      const CFreal y = qd0[iQd] * node0[YY] + qd1[iQd] * node1[YY];
      const CFreal L1 = 1.0 + 0.5*( ( x - x1 )*nx1 + ( y - y1 )*ny1 ) * inv_volume ;
      const CFreal L2 = 1.0 + 0.5*( ( x - x2 )*nx2 + ( y - y2 )*ny2 ) * inv_volume ;
      const CFreal L3 = 1.0 + 0.5*( ( x - x3 )*nx3 + ( y - y3 )*ny3 ) * inv_volume ;

      (*qdstates[iQd]) = (L1*( 2.0*L1 - 1.0 ) * state0) +
                         (L2*( 2.0*L2 - 1.0 ) * state1) +
                         (L3*( 2.0*L3 - 1.0 ) * state2) +
                         (4.0*L1*L2           * state3) +
                         (4.0*L3*L2           * state4) +
                         (4.0*L1*L3           * state5);

	(*qdnodes[iQd])[XX] = x;
	(*qdnodes[iQd])[YY] = y;

    }

    computeStatesData(3, m_updateVar, qdstates, m_pdata, m_qdExtraVars);

    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 3; ++iQd)
    {

      faceflux[iFace] += wqd[iQd] * m_updateVar->getFlux()(m_pdata[iQd],facenormal);
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
HOCRD_SplitStrategy::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result
    = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
