#include "FluctSplit/STM_HOCRD_SplitStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STM_HOCRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
spaceTimemhocrdFluctSplitStrategyProvider("STM_HOCRD");

//////////////////////////////////////////////////////////////////////////////

STM_HOCRD_SplitStrategy::STM_HOCRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  socket_updateCoeff("updateCoeff"),
  socket_pastNormals("pastNormals",false),
  socket_pastStates("pastStates"),
  socket_pastCellVolume("pastCellVolume",false),
  socket_cellSpeed("cellSpeed",false),
  m_splitter(CFNULL),
  _pastStates(0),
  _null(),
  m_flux_past_time(),
  m_flux_time(),
  m_flux_past_space(),
  m_flux_space(),
  CellVolume(),
  m_qdstates(0),
  m_qdExtraVars(0),
  m_contourIntegrator(CFNULL),
  m_nbQPointsInCell(),
  m_unitFaceNormals(),
  m_statesBkp(0),
  m_phi(),
  _pastResiduals(0),
  _pastResiduals_order1(0),
  temp_residual(0),
  m_phisubT(0),
  subelemfacedir(0,0),
  subelemtable(0,0),
  subfacetable(0,0),
  substates(0),
  subresidual(0),
  faceflux(0),
  qd0(0),
  qd1(0),
  wqd(0),
  qdstates(0),
  qdExtraVars(0),
  qdnodes(0),
  facenormal(0),
  m_updateVar(CFNULL),
  kappa1(0,0),
  kappa2(0,0),
  kappa3(0,0),
  kappa4(0,0),
  kappab1(0,0),
  kappab2(0,0),
  kappab3(0,0),
  kappab4(0,0),
  temp_vect0(0),
  temp_vect1(0),
  temp_vect2(0),
  temp_vect3(0),
  temp_vect4(0),
  temp_vect5(0)
{
}

//////////////////////////////////////////////////////////////////////////////

STM_HOCRD_SplitStrategy::~STM_HOCRD_SplitStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////
void STM_HOCRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  ///@todo Nv: There is no moving mesh implemented
  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  // Since there is no moving mesh here, we set the past cell volume to the same as
  // the non moving one
  CFreal volume = 2.0*volumes[ddata.cellID];
  m_splitter->setCellVolume(volumes[ddata.cellID]);
  m_splitter->setPastCellVolume(volumes[ddata.cellID]);
  m_splitter->setCellSpeed(_null);

  setCurrentCell();

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < residual.size(); ++i)
    {
    residual[i] = 0.0;
    }

  for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
    for (CFuint iEq = 0; iEq <  nbEqs; ++iEq){
      temp_residual[iState][iEq]=0.0;
    }
  }

  cf_assert(_pastStates.size() == 6); // P2 triangles for solution space
  cf_assert(residual.size() >= _pastStates.size()); // state residual
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				      [ddata.cellID]);

  cf_assert(cellnormals.nbFaces()  == 3); // triangles have 3 faces

  // If it is the first step, then first
  // compute the residuals due to the past
  if (SubSystemStatusStack::getActive()->isFirstStep()){
    DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();

    // get the paststates in this cell
    for (CFuint i = 0; i < nbStatesInCell; ++i) {
      const CFuint stateID = (*ddata.states)[i]->getLocalID();
      *_pastStates[i] = *pastStatesStorage[stateID];
    }

    // back up the flag telling if you are perturbing and set it to true
    // so that the update coefficient is not computed
    bool backUpPerturb = ddata.isPerturb;
    ddata.isPerturb = true;


    //We point the past_residual from the DistributeData to the past residual of cellID
    ddata.past_residuals = &_pastResiduals[ddata.cellID];

    //We point the past_residual of order1 from the DistributeData to the past residual of order 1 of cellID
    // If we are not unsing blending scheme these are vector of dimension 0
    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];

   //The time is the present time (this is used for the source term)
   ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();

 computeHOFluctuation_past();

    //We point the past_residual from the DistributeData to the past residual of cellID
    ddata.past_residuals = &_pastResiduals[ddata.cellID];

    /*****         Triangle 0-3-5          *****/

    substates[0] = _pastStates[0];
    substates[1] = _pastStates[3];
    substates[2] = _pastStates[5];

    ddata.tStates = computeConsistentStates(&substates);
    ddata.subStates = &substates;
    // normals are half scale
    cellnormals.scale(0.5);

    // compute the upwind parameters k in this cell
    m_splitter->computeK(substates,&cellnormals);
    //We unscale the normals to compute the source term
    cellnormals.unscale();
    ddata.sourceTermID = 0;

    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, _pastStates[iState]);
      }

    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
      }

    if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
      (*m_phisubT[0]) -= ddata.phiS;
    }

    // We point the current beta to beta of the first triangle
    SafePtr<vector<RealMatrix> >& currBetaMatrix = ddata.currBetaMat;
    currBetaMatrix = &ddata.betaMats[0];

    cf_assert(currBetaMatrix.isNotNull());

    // transform fluxes + source term to distribute variables and distribute
    SafePtr<RealVector> phi = &ddata.phi;
    *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
    m_splitter->distributePast(subresidual);
    getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

    // call beta of the sub-triangles, betas are always the one of the LDA schemes
    vector<RealMatrix> betasInTriag;
    betasInTriag.resize(3);
    for (CFuint i = 0 ; i< 3; ++i)
      betasInTriag[i].resize(4,4);

    betasInTriag = ddata.betaMats[0];

    temp_vect0 = betasInTriag[0]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[0]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[0]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[0]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[0]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[0]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[0][iEq] += 0.5*subresidual[0][iEq] -
                                  volume*(kappa1(0,0)*(*_pastStates[0])[iEq] + kappab1(0,0)*temp_vect0[iEq]  +
                                          kappa1(0,1)*(*_pastStates[1])[iEq] + kappab1(0,1)*temp_vect1[iEq]  +
                                          kappa1(0,2)*(*_pastStates[2])[iEq] + kappab1(0,2)*temp_vect2[iEq]  +
                                          kappa1(0,3)*(*_pastStates[3])[iEq] + kappab1(0,3)*temp_vect3[iEq]  +
                                          kappa1(0,4)*(*_pastStates[4])[iEq] + kappab1(0,4)*temp_vect4[iEq]  +
                                          kappa1(0,5)*(*_pastStates[5])[iEq] + kappab1(0,5)*temp_vect5[iEq]) ;


      temp_residual[1][iEq] -= volume*(kappa1(1,0)*(*_pastStates[0])[iEq] + kappa1(1,1)*(*_pastStates[1])[iEq] +
                                       kappa1(1,2)*(*_pastStates[2])[iEq] + kappa1(1,3)*(*_pastStates[3])[iEq] +
                                       kappa1(1,4)*(*_pastStates[4])[iEq] + kappa1(1,5)*(*_pastStates[5])[iEq]) ;

      temp_residual[2][iEq] -= volume*(kappa1(2,0)*(*_pastStates[0])[iEq] + kappa1(2,1)*(*_pastStates[1])[iEq] +
                                       kappa1(2,2)*(*_pastStates[2])[iEq] + kappa1(2,3)*(*_pastStates[3])[iEq] +
                                       kappa1(2,4)*(*_pastStates[4])[iEq] + kappa1(2,5)*(*_pastStates[5])[iEq]) ;
    }

    temp_vect0 = betasInTriag[1]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[1]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[1]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[1]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[1]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[1]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[3][iEq] += 0.5*subresidual[1][iEq] -
                                  volume*(kappa1(3,0)*(*_pastStates[0])[iEq] +  kappab1(1,0)*temp_vect0[iEq] +
                                          kappa1(3,1)*(*_pastStates[1])[iEq] +  kappab1(1,1)*temp_vect1[iEq] +
                                          kappa1(3,2)*(*_pastStates[2])[iEq] +  kappab1(1,2)*temp_vect2[iEq] +
                                          kappa1(3,3)*(*_pastStates[3])[iEq] +  kappab1(1,3)*temp_vect3[iEq] +
                                          kappa1(3,4)*(*_pastStates[4])[iEq] +  kappab1(1,4)*temp_vect4[iEq] +
                                          kappa1(3,5)*(*_pastStates[5])[iEq] +  kappab1(1,5)*temp_vect5[iEq] );

      temp_residual[4][iEq] -= volume*(kappa1(4,0)*(*_pastStates[0])[iEq] + kappa1(4,1)*(*_pastStates[1])[iEq] +
                                       kappa1(4,2)*(*_pastStates[2])[iEq] + kappa1(4,3)*(*_pastStates[3])[iEq] +
                                       kappa1(4,4)*(*_pastStates[4])[iEq] + kappa1(4,5)*(*_pastStates[5])[iEq]);
    }

    temp_vect0 = betasInTriag[2]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[2]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[2]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[2]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[2]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[2]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[5][iEq] += 0.5*subresidual[2][iEq] -
                                  volume*(kappa1(5,0)*(*_pastStates[0])[iEq] + kappab1(2,0)*temp_vect0[iEq] +
                                          kappa1(5,1)*(*_pastStates[1])[iEq] + kappab1(2,1)*temp_vect1[iEq] +
                                          kappa1(5,2)*(*_pastStates[2])[iEq] + kappab1(2,2)*temp_vect2[iEq] +
                                          kappa1(5,3)*(*_pastStates[3])[iEq] + kappab1(2,3)*temp_vect3[iEq] +
                                          kappa1(5,4)*(*_pastStates[4])[iEq] + kappab1(2,4)*temp_vect4[iEq] +
                                          kappa1(5,5)*(*_pastStates[5])[iEq] + kappab1(2,5)*temp_vect5[iEq]);

    }


    /*****         Triangle 3-1-4          *****/

    substates[0] = _pastStates[3];
    substates[1] = _pastStates[1];
    substates[2] = _pastStates[4];

    ddata.tStates = computeConsistentStates(&substates);
    ddata.subStates = &substates;
    // normals are half scale
    cellnormals.scale(0.5);

    // compute the residual and the upwind parameters k in this cell
    m_splitter->computeK(substates,&cellnormals);

    cellnormals.unscale();
    //We unscale the normals to compute the source term
    ddata.sourceTermID = 0;
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, _pastStates[iState]);
    }
    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
    }
    if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
      (*m_phisubT[1]) -= ddata.phiS;
    }

    // We point the current beta to beta of the second triangle
    currBetaMatrix = &ddata.betaMats[1];
    cf_assert(currBetaMatrix.isNotNull());

    // transform fluxes + source term to distribute variables and distribute
    *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);
    m_splitter->distributePast(subresidual);
    getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

    betasInTriag = ddata.betaMats[1];

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[0][iEq] -= volume*(kappa2(0,0)*(*_pastStates[0])[iEq] + kappa2(0,1)*(*_pastStates[1])[iEq] +
                                       kappa2(0,2)*(*_pastStates[2])[iEq] + kappa2(0,3)*(*_pastStates[3])[iEq] +
                                       kappa2(0,4)*(*_pastStates[4])[iEq] + kappa2(0,5)*(*_pastStates[5])[iEq]);
    }

    temp_vect0 = betasInTriag[1]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[1]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[1]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[1]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[1]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[1]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[1][iEq] += 0.5*subresidual[1][iEq] -
                                  volume*(kappa2(1,0)*(*_pastStates[0])[iEq] + kappab2(1,0)*temp_vect0[iEq] +
                                          kappa2(1,1)*(*_pastStates[1])[iEq] + kappab2(1,1)*temp_vect1[iEq] +
                                          kappa2(1,2)*(*_pastStates[2])[iEq] + kappab2(1,2)*temp_vect2[iEq] +
                                          kappa2(1,3)*(*_pastStates[3])[iEq] + kappab2(1,3)*temp_vect3[iEq] +
                                          kappa2(1,4)*(*_pastStates[4])[iEq] + kappab2(1,4)*temp_vect4[iEq] +
                                          kappa2(1,5)*(*_pastStates[5])[iEq] + kappab2(1,5)*temp_vect5[iEq] );

      temp_residual[2][iEq] -= volume*(kappa2(2,0)*(*_pastStates[0])[iEq] + kappa2(2,1)*(*_pastStates[1])[iEq] +
                                       kappa2(2,2)*(*_pastStates[2])[iEq] + kappa2(2,3)*(*_pastStates[3])[iEq] +
                                       kappa2(2,4)*(*_pastStates[4])[iEq] + kappa2(2,5)*(*_pastStates[5])[iEq]) ;
    }

    temp_vect0 = betasInTriag[0]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[0]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[0]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[0]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[0]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[0]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[3][iEq] += 0.5*subresidual[0][iEq] -
                                 volume*(kappa2(3,0)*(*_pastStates[0])[iEq] + kappab2(0,0)*temp_vect0[iEq] +
                                         kappa2(3,1)*(*_pastStates[1])[iEq] + kappab2(0,1)*temp_vect1[iEq] +
                                         kappa2(3,2)*(*_pastStates[2])[iEq] + kappab2(0,2)*temp_vect2[iEq] +
                                         kappa2(3,3)*(*_pastStates[3])[iEq] + kappab2(0,3)*temp_vect3[iEq] +
                                         kappa2(3,4)*(*_pastStates[4])[iEq] + kappab2(0,4)*temp_vect4[iEq] +
                                         kappa2(3,5)*(*_pastStates[5])[iEq] + kappab2(0,5)*temp_vect5[iEq] );
    }

    temp_vect0 = betasInTriag[2]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[2]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[2]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[2]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[2]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[2]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[4][iEq] += 0.5*subresidual[2][iEq] -
                                  volume*(kappa2(4,0)*(*_pastStates[0])[iEq] + kappab2(2,0)*temp_vect0[iEq] +
                                          kappa2(4,1)*(*_pastStates[1])[iEq] + kappab2(2,1)*temp_vect1[iEq] +
                                          kappa2(4,2)*(*_pastStates[2])[iEq] + kappab2(2,2)*temp_vect2[iEq] +
                                          kappa2(4,3)*(*_pastStates[3])[iEq] + kappab2(2,3)*temp_vect3[iEq] +
                                          kappa2(4,4)*(*_pastStates[4])[iEq] + kappab2(2,4)*temp_vect4[iEq] +
                                          kappa2(4,5)*(*_pastStates[5])[iEq] + kappab2(2,5)*temp_vect5[iEq]);

      temp_residual[5][iEq] -= volume*(kappa2(5,0)*(*_pastStates[0])[iEq] + kappa2(5,1)*(*_pastStates[1])[iEq] +
                                       kappa2(5,2)*(*_pastStates[2])[iEq] + kappa2(5,3)*(*_pastStates[3])[iEq] +
                                       kappa2(5,4)*(*_pastStates[4])[iEq] + kappa2(5,5)*(*_pastStates[5])[iEq]) ;

    }


    /*****         Triangle 5-4-2          *****/

    substates[0] = _pastStates[5];
    substates[1] = _pastStates[4];
    substates[2] = _pastStates[2];

    ddata.tStates = computeConsistentStates(&substates);

    // normals are half scale
    cellnormals.scale(0.5);

    // compute the residual and the upwind parameters k in this cell
    m_splitter->computeK(substates,&cellnormals);

    //We unscale the normals to compute the source term
    cellnormals.unscale();
    ddata.sourceTermID = 0;
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, _pastStates[iState]);
      }
    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
      }
    if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
     (*m_phisubT[2]) -= ddata.phiS;
    }

    // We point the current beta to beta of the third triangle
    currBetaMatrix = &ddata.betaMats[2];

    // transform fluxes + source term to distribute variables and distribute
    *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);
    m_splitter->distributePast(subresidual);
    getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

    betasInTriag = ddata.betaMats[2];

    temp_vect0 = betasInTriag[2]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[2]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[2]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[2]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[2]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[2]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    temp_residual[0][iEq] -= volume*(kappa3(0,0)*(*_pastStates[0])[iEq] + kappa3(0,1)*(*_pastStates[1])[iEq] +
                                     kappa3(0,2)*(*_pastStates[2])[iEq] + kappa3(0,3)*(*_pastStates[3])[iEq] +
                                     kappa3(0,4)*(*_pastStates[4])[iEq] + kappa3(0,5)*(*_pastStates[5])[iEq]);

    temp_residual[1][iEq] -= volume*(kappa3(1,0)*(*_pastStates[0])[iEq] + kappa3(1,1)*(*_pastStates[1])[iEq] +
                                     kappa3(1,2)*(*_pastStates[2])[iEq] + kappa3(1,3)*(*_pastStates[3])[iEq] +
                                     kappa3(1,4)*(*_pastStates[4])[iEq] + kappa3(1,5)*(*_pastStates[5])[iEq]);

    temp_residual[2][iEq] += 0.5*subresidual[2][iEq] -
                                  volume*(kappa3(2,0)*(*_pastStates[0])[iEq] + kappab3(2,0)*temp_vect0[iEq] +
                                          kappa3(2,1)*(*_pastStates[1])[iEq] + kappab3(2,1)*temp_vect1[iEq] +
                                          kappa3(2,2)*(*_pastStates[2])[iEq] + kappab3(2,2)*temp_vect2[iEq] +
                                          kappa3(2,3)*(*_pastStates[3])[iEq] + kappab3(2,3)*temp_vect3[iEq] +
                                          kappa3(2,4)*(*_pastStates[4])[iEq] + kappab3(2,4)*temp_vect4[iEq] +
                                          kappa3(2,5)*(*_pastStates[5])[iEq] + kappab3(2,5)*temp_vect5[iEq]);


    temp_residual[3][iEq] -= volume*(kappa3(3,0)*(*_pastStates[0])[iEq] + kappa3(3,1)*(*_pastStates[1])[iEq] +
                                     kappa3(3,2)*(*_pastStates[2])[iEq] + kappa3(3,3)*(*_pastStates[3])[iEq] +
                                     kappa3(3,4)*(*_pastStates[4])[iEq] + kappa3(3,5)*(*_pastStates[5])[iEq]);
    }


    temp_vect0 = betasInTriag[1]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[1]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[1]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[1]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[1]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[1]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

      temp_residual[4][iEq] += 0.5*subresidual[1][iEq] -
                                  volume*(kappa3(4,0)*(*_pastStates[0])[iEq] + kappab3(1,0)*temp_vect0[iEq] +
                                          kappa3(4,1)*(*_pastStates[1])[iEq] + kappab3(1,1)*temp_vect1[iEq] +
                                          kappa3(4,2)*(*_pastStates[2])[iEq] + kappab3(1,2)*temp_vect2[iEq] +
                                          kappa3(4,3)*(*_pastStates[3])[iEq] + kappab3(1,3)*temp_vect3[iEq] +
                                          kappa3(4,4)*(*_pastStates[4])[iEq] + kappab3(1,4)*temp_vect4[iEq] +
                                          kappa3(4,5)*(*_pastStates[5])[iEq] + kappab3(1,5)*temp_vect5[iEq]) ;

    }

    temp_vect0 = betasInTriag[0]*(*_pastStates[0]);
    temp_vect1 = betasInTriag[0]*(*_pastStates[1]);
    temp_vect2 = betasInTriag[0]*(*_pastStates[2]);
    temp_vect3 = betasInTriag[0]*(*_pastStates[3]);
    temp_vect4 = betasInTriag[0]*(*_pastStates[4]);
    temp_vect5 = betasInTriag[0]*(*_pastStates[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      temp_residual[5][iEq] += 0.5*subresidual[0][iEq] -
                                  volume*(kappa3(5,0)*(*_pastStates[0])[iEq] + kappab3(0,0)*temp_vect0[iEq] +
                                          kappa3(5,1)*(*_pastStates[1])[iEq] + kappab3(0,1)*temp_vect1[iEq] +
                                          kappa3(5,2)*(*_pastStates[2])[iEq] + kappab3(0,2)*temp_vect2[iEq] +
                                          kappa3(5,3)*(*_pastStates[3])[iEq] + kappab3(0,3)*temp_vect3[iEq] +
                                          kappa3(5,4)*(*_pastStates[4])[iEq] + kappab3(0,4)*temp_vect4[iEq] +
                                          kappa3(5,5)*(*_pastStates[5])[iEq] + kappab3(0,5)*temp_vect5[iEq]) ;
    }


   /*****         Triangle 4-5-3          *****/

   substates[0] = _pastStates[4];
   substates[1] = _pastStates[5];
   substates[2] = _pastStates[3];

   ddata.tStates = computeConsistentStates(&substates);
   ddata.subStates = &substates;
   // compute the upwind parameters k in this cell
   // the oriantation of the normal in this sub-element is oposite to the one of the element
   cellnormals.scale(-0.5);
   m_splitter->computeK(substates,&cellnormals);

   //   We unscale the normals to compute the source term
   cellnormals.unscale();
   ddata.sourceTermID = 0;
   for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, _pastStates[iState]);
   }
   getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);
   for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        ddata.cell->setState(iState, m_statesBkp[iState]);
   }
   if (getMethodData().includeSourceInFlux()) {
     // in this case converctive and source fluctuations will be distributed together
     (*m_phisubT[3]) -= ddata.phiS;
   }

   // We point the current beta to beta of the third triangle
   currBetaMatrix = &ddata.betaMats[3];

   // transform fluxes + source term to distribute variables and distribute
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
   m_splitter->distributePast(subresidual);
   getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

   betasInTriag = ddata.betaMats[3];

   temp_vect0 = betasInTriag[2]*(*_pastStates[0]);
   temp_vect1 = betasInTriag[2]*(*_pastStates[1]);
   temp_vect2 = betasInTriag[2]*(*_pastStates[2]);
   temp_vect3 = betasInTriag[2]*(*_pastStates[3]);
   temp_vect4 = betasInTriag[2]*(*_pastStates[4]);
   temp_vect5 = betasInTriag[2]*(*_pastStates[5]);

   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     temp_residual[0][iEq] -= volume*(kappa4(0,0)*(*_pastStates[0])[iEq] + kappa4(0,1)*(*_pastStates[1])[iEq] +
                                      kappa4(0,2)*(*_pastStates[2])[iEq] + kappa4(0,3)*(*_pastStates[3])[iEq] +
                                      kappa4(0,4)*(*_pastStates[4])[iEq] + kappa4(0,5)*(*_pastStates[5])[iEq]);

     temp_residual[1][iEq] -= volume*(kappa4(1,0)*(*_pastStates[0])[iEq] + kappa4(1,1)*(*_pastStates[1])[iEq] +
                                      kappa4(1,2)*(*_pastStates[2])[iEq] + kappa4(1,3)*(*_pastStates[3])[iEq] +
                                      kappa4(1,4)*(*_pastStates[4])[iEq] + kappa4(1,5)*(*_pastStates[5])[iEq]);

     temp_residual[2][iEq] -= volume*(kappa4(2,0)*(*_pastStates[0])[iEq] + kappa4(2,1)*(*_pastStates[1])[iEq] +
                                      kappa4(2,2)*(*_pastStates[2])[iEq] + kappa4(2,3)*(*_pastStates[3])[iEq] +
                                      kappa4(2,4)*(*_pastStates[4])[iEq] + kappa4(2,5)*(*_pastStates[5])[iEq]);

     temp_residual[3][iEq] += 0.5*subresidual[2][iEq] -
                                  volume*(kappa4(3,0)*(*_pastStates[0])[iEq] + kappab4(2,0)*temp_vect0[iEq] +
                                          kappa4(3,1)*(*_pastStates[1])[iEq] + kappab4(2,1)*temp_vect1[iEq] +
                                          kappa4(3,2)*(*_pastStates[2])[iEq] + kappab4(2,2)*temp_vect2[iEq] +
                                          kappa4(3,3)*(*_pastStates[3])[iEq] + kappab4(2,3)*temp_vect3[iEq] +
                                          kappa4(3,4)*(*_pastStates[4])[iEq] + kappab4(2,4)*temp_vect4[iEq] +
                                          kappa4(3,5)*(*_pastStates[5])[iEq] + kappab4(2,5)*temp_vect5[iEq]) ;
   }

   temp_vect0 = betasInTriag[0]*(*_pastStates[0]);
   temp_vect1 = betasInTriag[0]*(*_pastStates[1]);
   temp_vect2 = betasInTriag[0]*(*_pastStates[2]);
   temp_vect3 = betasInTriag[0]*(*_pastStates[3]);
   temp_vect4 = betasInTriag[0]*(*_pastStates[4]);
   temp_vect5 = betasInTriag[0]*(*_pastStates[5]);

   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     temp_residual[4][iEq] += 0.5*subresidual[0][iEq] -
                                  volume*(kappa4(4,0)*(*_pastStates[0])[iEq] + kappab4(0,0)*temp_vect0[iEq] +
                                          kappa4(4,1)*(*_pastStates[1])[iEq] + kappab4(0,1)*temp_vect1[iEq] +
                                          kappa4(4,2)*(*_pastStates[2])[iEq] + kappab4(0,2)*temp_vect2[iEq] +
                                          kappa4(4,3)*(*_pastStates[3])[iEq] + kappab4(0,3)*temp_vect3[iEq] +
                                          kappa4(4,4)*(*_pastStates[4])[iEq] + kappab4(0,4)*temp_vect4[iEq] +
                                          kappa4(4,5)*(*_pastStates[5])[iEq] + kappab4(0,5)*temp_vect5[iEq]) ;
   }

   temp_vect0 = betasInTriag[1]*(*_pastStates[0]);
   temp_vect1 = betasInTriag[1]*(*_pastStates[1]);
   temp_vect2 = betasInTriag[1]*(*_pastStates[2]);
   temp_vect3 = betasInTriag[1]*(*_pastStates[3]);
   temp_vect4 = betasInTriag[1]*(*_pastStates[4]);
   temp_vect5 = betasInTriag[1]*(*_pastStates[5]);

   for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
     temp_residual[5][iEq] += 0.5*subresidual[1][iEq] -
                                 volume*( kappa4(5,0)*(*_pastStates[0])[iEq] + kappab4(1,0)*temp_vect0[iEq] +
                                          kappa4(5,1)*(*_pastStates[1])[iEq] + kappab4(1,1)*temp_vect1[iEq] +
                                          kappa4(5,2)*(*_pastStates[2])[iEq] + kappab4(1,2)*temp_vect2[iEq] +
                                          kappa4(5,3)*(*_pastStates[3])[iEq] + kappab4(1,3)*temp_vect3[iEq] +
                                          kappa4(5,4)*(*_pastStates[4])[iEq] + kappab4(1,4)*temp_vect4[iEq] +
                                          kappa4(5,5)*(*_pastStates[5])[iEq] + kappab4(1,5)*temp_vect5[iEq]) ;

  }

   //We unscale at the end!
   cellnormals.unscale();

   // restore the backed up flag
   ddata.isPerturb = backUpPerturb;

   RealVector& past_residuals = *ddata.past_residuals;

   for (CFuint iState = 0; iState < 6; ++ iState)
     {
      for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
        {
         past_residuals[iState*nbEqs+iEq] = temp_residual[iState][iEq];
        }
     }

  }

  if (!SubSystemStatusStack::getActive()->isFirstStep())
   {
    // If this is not the first pseudo-time iteration we point to the past_residual of the cell
    ddata.past_residuals = &_pastResiduals[ddata.cellID];
    //    ddata.past_residuals_order1 = &_pastResiduals_order1[ddata.cellID];
   }


  const RealVector& past_residuals = *ddata.past_residuals;

  for (CFuint iState = 0; iState < 6; ++ iState)
    {
     for (CFuint iEq = 0; iEq < nbEqs; ++ iEq)
      {
        residual[iState][iEq] = past_residuals[iState*nbEqs+ iEq];
      }
    }

  vector<State*>& states = *ddata.states;


  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  //The time is the present time (this is used for the source term)
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime()+dt;

  computeHOFluctuation_present();

  /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  // normals are half scale
           cellnormals.scale(0.5);

  // compute the upwind parameters k in this cell
  m_splitter->computeK(substates,&cellnormals);

   //We unscale the normals to compute the source term
   cellnormals.unscale();
   ddata.sourceTermID = 0;
   getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

   if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
      (*m_phisubT[0]) -= ddata.phiS;
   }

   // We point the current beta to beta of the first triangle
   SafePtr<vector<RealMatrix> >& currBetaMatrix = ddata.currBetaMat;
   currBetaMatrix = &ddata.betaMats[0];

   cf_assert(currBetaMatrix.isNotNull());

   // transform fluxes of subelement to distribution variables
   SafePtr<RealVector> phi = &ddata.phi;
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[0]);
   m_splitter->distribute(subresidual);
   getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

  // call beta of the sub-triangles, betas are always the one of the LDA schemes
  vector<RealMatrix> betasInTriag;
  betasInTriag.resize(3);
  for (CFuint i = 0 ; i< 3; ++i)
    betasInTriag[i].resize(4,4);


  betasInTriag = ddata.betaMats[0];

    temp_vect0 = betasInTriag[0]*(*states[0]);
    temp_vect1 = betasInTriag[0]*(*states[1]);
    temp_vect2 = betasInTriag[0]*(*states[2]);
    temp_vect3 = betasInTriag[0]*(*states[3]);
    temp_vect4 = betasInTriag[0]*(*states[4]);
    temp_vect5 = betasInTriag[0]*(*states[5]);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
  residual[0][iEq] += 0.5*subresidual[0][iEq] +
                                  volume*(kappa1(0,0)*(*states[0])[iEq] + kappab1(0,0)*temp_vect0[iEq] +
                                          kappa1(0,1)*(*states[1])[iEq] + kappab1(0,1)*temp_vect1[iEq] +
                                          kappa1(0,2)*(*states[2])[iEq] + kappab1(0,2)*temp_vect2[iEq] +
                                          kappa1(0,3)*(*states[3])[iEq] + kappab1(0,3)*temp_vect3[iEq] +
                                          kappa1(0,4)*(*states[4])[iEq] + kappab1(0,4)*temp_vect4[iEq] +
                                          kappa1(0,5)*(*states[5])[iEq] + kappab1(0,5)*temp_vect5[iEq]) ;


   residual[1][iEq] += volume*(kappa1(1,0)*(*states[0])[iEq] + kappa1(1,1)*(*states[1])[iEq] +
                          kappa1(1,2)*(*states[2])[iEq] + kappa1(1,3)*(*states[3])[iEq] +
                          kappa1(1,4)*(*states[4])[iEq] + kappa1(1,5)*(*states[5])[iEq]) ;

   residual[2][iEq] += volume*(kappa1(2,0)*(*states[0])[iEq] + kappa1(2,1)*(*states[1])[iEq] +
                          kappa1(2,2)*(*states[2])[iEq] + kappa1(2,3)*(*states[3])[iEq] +
                          kappa1(2,4)*(*states[4])[iEq] + kappa1(2,5)*(*states[5])[iEq]) ;
}

 temp_vect0 = betasInTriag[1]*(*states[0]);
 temp_vect1 = betasInTriag[1]*(*states[1]);
 temp_vect2 = betasInTriag[1]*(*states[2]);
 temp_vect3 = betasInTriag[1]*(*states[3]);
 temp_vect4 = betasInTriag[1]*(*states[4]);
 temp_vect5 = betasInTriag[1]*(*states[5]);

         for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[3][iEq] += 0.5*subresidual[1][iEq] +
                                  volume*(kappa1(3,0)*(*states[0])[iEq] + kappab1(1,0)*temp_vect0[iEq] +
                                          kappa1(3,1)*(*states[1])[iEq] + kappab1(1,1)*temp_vect1[iEq] +
                                          kappa1(3,2)*(*states[2])[iEq] + kappab1(1,2)*temp_vect2[iEq] +
                                          kappa1(3,3)*(*states[3])[iEq] + kappab1(1,3)*temp_vect3[iEq] +
                                          kappa1(3,4)*(*states[4])[iEq] + kappab1(1,4)*temp_vect4[iEq] +
                                          kappa1(3,5)*(*states[5])[iEq] + kappab1(1,5)*temp_vect5[iEq]) ;



   residual[4][iEq] += volume*(kappa1(4,0)*(*states[0])[iEq] + kappa1(4,1)*(*states[1])[iEq] +
                          kappa1(4,2)*(*states[2])[iEq] + kappa1(4,3)*(*states[3])[iEq] +
                          kappa1(4,4)*(*states[4])[iEq] + kappa1(4,5)*(*states[5])[iEq]) ;
               }

temp_vect0 = betasInTriag[2]*(*states[0]);
temp_vect1 = betasInTriag[2]*(*states[1]);
temp_vect2 = betasInTriag[2]*(*states[2]);
temp_vect3 = betasInTriag[2]*(*states[3]);
temp_vect4 = betasInTriag[2]*(*states[4]);
temp_vect5 = betasInTriag[2]*(*states[5]);
 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[5][iEq] += 0.5*subresidual[2][iEq] +
                                  volume*(kappa1(5,0)*(*states[0])[iEq] + kappab1(2,0)*temp_vect0[iEq] +
                                          kappa1(5,1)*(*states[1])[iEq] + kappab1(2,1)*temp_vect1[iEq] +
                                          kappa1(5,2)*(*states[2])[iEq] + kappab1(2,2)*temp_vect2[iEq] +
                                          kappa1(5,3)*(*states[3])[iEq] + kappab1(2,3)*temp_vect3[iEq] +
                                          kappa1(5,4)*(*states[4])[iEq] + kappab1(2,4)*temp_vect4[iEq] +
                                          kappa1(5,5)*(*states[5])[iEq] + kappab1(2,5)*temp_vect5[iEq]) ;
  }


   /*****         Triangle 3-1-4          *****/

   substates[0] = states[3];
   substates[1] = states[1];
   substates[2] = states[4];

   ddata.tStates = computeConsistentStates(&substates);
    ddata.subStates = &substates;
   // normals are half scale
 cellnormals.scale(0.5);

   // compute  the upwind parameters k in this cell
   m_splitter->computeK(substates,&cellnormals);

   cellnormals.unscale();
   //We unscale the normals to compute the source term
   ddata.sourceTermID = 0;
   getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

   if (getMethodData().includeSourceInFlux()) {
   // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[1]) -= ddata.phiS;
   }

   // We point the current beta to beta of the second triangle
   currBetaMatrix = &ddata.betaMats[1];
   cf_assert(currBetaMatrix.isNotNull());

   // transform fluxes + source term to distribute variables and distribute
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[1]);
   m_splitter->distribute(subresidual);
   getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

   betasInTriag = ddata.betaMats[1];
temp_vect0 = betasInTriag[1]*(*states[0]);
temp_vect1 = betasInTriag[1]*(*states[1]);
temp_vect2 = betasInTriag[1]*(*states[2]);
temp_vect3 = betasInTriag[1]*(*states[3]);
temp_vect4 = betasInTriag[1]*(*states[4]);
temp_vect5 = betasInTriag[1]*(*states[5]);
 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[0][iEq] += volume*(kappa2(0,0)*(*states[0])[iEq] + kappa2(0,1)*(*states[1])[iEq] +
                          kappa2(0,2)*(*states[2])[iEq] + kappa2(0,3)*(*states[3])[iEq] +
                          kappa2(0,4)*(*states[4])[iEq] + kappa2(0,5)*(*states[5])[iEq]);

  residual[1][iEq]+= 0.5*subresidual[1][iEq] +
                                  volume*(kappa2(1,0)*(*states[0])[iEq] + kappab2(1,0)*temp_vect0[iEq] +
                                          kappa2(1,1)*(*states[1])[iEq] + kappab2(1,1)*temp_vect1[iEq] +
                                          kappa2(1,2)*(*states[2])[iEq] + kappab2(1,2)*temp_vect2[iEq] +
                                          kappa2(1,3)*(*states[3])[iEq] + kappab2(1,3)*temp_vect3[iEq] +
                                          kappa2(1,4)*(*states[4])[iEq] + kappab2(1,4)*temp_vect4[iEq] +
                                          kappa2(1,5)*(*states[5])[iEq] + kappab2(1,5)*temp_vect5[iEq]) ;



   residual[2][iEq] += volume*(kappa2(2,0)*(*states[0])[iEq] + kappa2(2,1)*(*states[1])[iEq] +
                          kappa2(2,2)*(*states[2])[iEq] + kappa2(2,3)*(*states[3])[iEq] +
                          kappa2(2,4)*(*states[4])[iEq] + kappa2(2,5)*(*states[5])[iEq]) ;
}

temp_vect0 = betasInTriag[0]*(*states[0]);
temp_vect1 = betasInTriag[0]*(*states[1]);
temp_vect2 = betasInTriag[0]*(*states[2]);
temp_vect3 = betasInTriag[0]*(*states[3]);
temp_vect4 = betasInTriag[0]*(*states[4]);
temp_vect5 = betasInTriag[0]*(*states[5]);
 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[3][iEq] += 0.5*subresidual[0][iEq] +
                                 volume*(kappa2(3,0)*(*states[0])[iEq] + kappab2(0,0)*temp_vect0[iEq] +
                                         kappa2(3,1)*(*states[1])[iEq] + kappab2(0,1)*temp_vect1[iEq] +
                                         kappa2(3,2)*(*states[2])[iEq] + kappab2(0,2)*temp_vect2[iEq] +
                                         kappa2(3,3)*(*states[3])[iEq] + kappab2(0,3)*temp_vect3[iEq] +
                                         kappa2(3,4)*(*states[4])[iEq] + kappab2(0,4)*temp_vect4[iEq] +
                                         kappa2(3,5)*(*states[5])[iEq] + kappab2(0,5)*temp_vect5[iEq] );
}


temp_vect0 = betasInTriag[2]*(*states[0]);
temp_vect1 = betasInTriag[2]*(*states[1]);
temp_vect2 = betasInTriag[2]*(*states[2]);
temp_vect3 = betasInTriag[2]*(*states[3]);
temp_vect4 = betasInTriag[2]*(*states[4]);
temp_vect5 = betasInTriag[2]*(*states[5]);
 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[4][iEq] += 0.5*subresidual[2][iEq] +
                                  volume*(kappa2(4,0)*(*states[0])[iEq] + kappab2(2,0)*temp_vect0[iEq] +
                                          kappa2(4,1)*(*states[1])[iEq] + kappab2(2,1)*temp_vect1[iEq] +
                                          kappa2(4,2)*(*states[2])[iEq] + kappab2(2,2)*temp_vect2[iEq] +
                                          kappa2(4,3)*(*states[3])[iEq] + kappab2(2,3)*temp_vect3[iEq] +
                                          kappa2(4,4)*(*states[4])[iEq] + kappab2(2,4)*temp_vect4[iEq] +
                                          kappa2(4,5)*(*states[5])[iEq] + kappab2(2,5)*temp_vect5[iEq] );


   residual[5][iEq] += volume*(kappa2(5,0)*(*states[0])[iEq] + kappa2(5,1)*(*states[1])[iEq] +
                          kappa2(5,2)*(*states[2])[iEq] + kappa2(5,3)*(*states[3])[iEq] +
                          kappa2(5,4)*(*states[4])[iEq] + kappa2(5,5)*(*states[5])[iEq]) ;
}


   /*****         Triangle 5-4-2          *****/

   substates[0] = states[5];
   substates[1] = states[4];
   substates[2] = states[2];

   ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;
   // normals are half scale
   cellnormals.scale(0.5);

   // compute the upwind parameters k in this cell
   m_splitter->computeK(substates,&cellnormals);

   //We unscale the normals to compute the source term
   cellnormals.unscale();
   ddata.sourceTermID = 0;
   getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

   if (getMethodData().includeSourceInFlux()) {
   // in this case converctive and source fluctuations will be distributed together
    (*m_phisubT[2]) -= ddata.phiS;
   }

   // We point the current beta to beta of the third triangle
   currBetaMatrix = &ddata.betaMats[2];


   // transform fluxes of subelement to distribution variables
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[2]);
   m_splitter->distribute(subresidual);
   getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

   betasInTriag = ddata.betaMats[2];

temp_vect0 = betasInTriag[2]*(*states[0]);
temp_vect1 = betasInTriag[2]*(*states[1]);
temp_vect2 = betasInTriag[2]*(*states[2]);
temp_vect3 = betasInTriag[2]*(*states[3]);
temp_vect4 = betasInTriag[2]*(*states[4]);
temp_vect5 = betasInTriag[2]*(*states[5]);
 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[0][iEq] += volume*(kappa3(0,0)*(*states[0])[iEq] + kappa3(0,1)*(*states[1])[iEq] +
                          kappa3(0,2)*(*states[2])[iEq] + kappa3(0,3)*(*states[3])[iEq] +
                          kappa3(0,4)*(*states[4])[iEq] + kappa3(0,5)*(*states[5])[iEq]);

   residual[1][iEq] += volume*(kappa3(1,0)*(*states[0])[iEq] + kappa3(1,1)*(*states[1])[iEq] +
                          kappa3(1,2)*(*states[2])[iEq] + kappa3(1,3)*(*states[3])[iEq] +
                          kappa3(1,4)*(*states[4])[iEq] + kappa3(1,5)*(*states[5])[iEq]);

   residual[2][iEq] += 0.5*subresidual[2][iEq] +
                                  volume*(kappa3(2,0)*(*states[0])[iEq] + kappab3(2,0)*temp_vect0[iEq] +
                                          kappa3(2,1)*(*states[1])[iEq] + kappab3(2,1)*temp_vect1[iEq] +
                                          kappa3(2,2)*(*states[2])[iEq] + kappab3(2,2)*temp_vect2[iEq] +
                                          kappa3(2,3)*(*states[3])[iEq] + kappab3(2,3)*temp_vect3[iEq] +
                                          kappa3(2,4)*(*states[4])[iEq] + kappab3(2,4)*temp_vect4[iEq] +
                                          kappa3(2,5)*(*states[5])[iEq] + kappab3(2,5)*temp_vect5[iEq] ) ;

   residual[3][iEq] += volume*(kappa3(3,0)*(*states[0])[iEq] + kappa3(3,1)*(*states[1])[iEq] +
                          kappa3(3,2)*(*states[2])[iEq] + kappa3(3,3)*(*states[3])[iEq] +
                          kappa3(3,4)*(*states[4])[iEq] + kappa3(3,5)*(*states[5])[iEq]);
    }

temp_vect0 = betasInTriag[1]*(*states[0]);
temp_vect1 = betasInTriag[1]*(*states[1]);
temp_vect2 = betasInTriag[1]*(*states[2]);
temp_vect3 = betasInTriag[1]*(*states[3]);
temp_vect4 = betasInTriag[1]*(*states[4]);
temp_vect5 = betasInTriag[1]*(*states[5]);

 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[4][iEq] += 0.5*subresidual[1][iEq]+
                                  volume*(kappa3(4,0)*(*states[0])[iEq] + kappab3(1,0)*temp_vect0[iEq] +
                                          kappa3(4,1)*(*states[1])[iEq] + kappab3(1,1)*temp_vect1[iEq] +
                                          kappa3(4,2)*(*states[2])[iEq] + kappab3(1,2)*temp_vect2[iEq] +
                                          kappa3(4,3)*(*states[3])[iEq] + kappab3(1,3)*temp_vect3[iEq] +
                                          kappa3(4,4)*(*states[4])[iEq] + kappab3(1,4)*temp_vect4[iEq] +
                                          kappa3(4,5)*(*states[5])[iEq] + kappab3(1,5)*temp_vect5[iEq] ) ;
}


temp_vect0 = betasInTriag[0]*(*states[0]);
temp_vect1 = betasInTriag[0]*(*states[1]);
temp_vect2 = betasInTriag[0]*(*states[2]);
temp_vect3 = betasInTriag[0]*(*states[3]);
temp_vect4 = betasInTriag[0]*(*states[4]);
temp_vect5 = betasInTriag[0]*(*states[5]);

 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {

   residual[5][iEq] += 0.5*subresidual[0][iEq]+
                                  volume*(kappa3(5,0)*(*states[0])[iEq] + kappab3(0,0)*temp_vect0[iEq] +
                                          kappa3(5,1)*(*states[1])[iEq] + kappab3(0,1)*temp_vect1[iEq] +
                                          kappa3(5,2)*(*states[2])[iEq] + kappab3(0,2)*temp_vect2[iEq] +
                                          kappa3(5,3)*(*states[3])[iEq] + kappab3(0,3)*temp_vect3[iEq] +
                                          kappa3(5,4)*(*states[4])[iEq] + kappab3(0,4)*temp_vect4[iEq] +
                                          kappa3(5,5)*(*states[5])[iEq] + kappab3(0,5)*temp_vect5[iEq]) ;
}

   /*****         Triangle 4-5-3          *****/

   substates[0] = states[4];
   substates[1] = states[5];
   substates[2] = states[3];

   ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;
   // compute the residual and the upwind parameters k in this cell
   // the oriantation of the normal in this sub-element is oposite to the one of the element
   cellnormals.scale(-0.5);
   m_splitter->computeK(substates,&cellnormals);

   //We unscale the normals to compute the source term
   cellnormals.unscale();
   ddata.sourceTermID = 0;
   getMethodData().getSourceTermSplitter(0)->computeSourceTerm(cellnormals);

   if (getMethodData().includeSourceInFlux()) {
     // in this case converctive and source fluctuations will be distributed together
     (*m_phisubT[3]) -= ddata.phiS;
   }

   // We point the current beta to beta of the third triangle
   currBetaMatrix = &ddata.betaMats[3];

   // transform fluxes of subelement to distribution variables
   *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(m_phisubT[3]);
   m_splitter->distribute(subresidual);
   getMethodData().getSourceTermSplitter(0)->distribute(subresidual);

   betasInTriag = ddata.betaMats[3];

temp_vect0 = betasInTriag[2]*(*states[0]);
temp_vect1 = betasInTriag[2]*(*states[1]);
temp_vect2 = betasInTriag[2]*(*states[2]);
temp_vect3 = betasInTriag[2]*(*states[3]);
temp_vect4 = betasInTriag[2]*(*states[4]);
temp_vect5 = betasInTriag[2]*(*states[5]);

 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[0][iEq] += volume*(kappa4(0,0)*(*states[0])[iEq] + kappa4(0,1)*(*states[1])[iEq] +
                          kappa4(0,2)*(*states[2])[iEq] + kappa4(0,3)*(*states[3])[iEq] +
                          kappa4(0,4)*(*states[4])[iEq] + kappa4(0,5)*(*states[5])[iEq]);

   residual[1][iEq] += volume*(kappa4(1,0)*(*states[0])[iEq] + kappa4(1,1)*(*states[1])[iEq] +
                          kappa4(1,2)*(*states[2])[iEq] + kappa4(1,3)*(*states[3])[iEq] +
                          kappa4(1,4)*(*states[4])[iEq] + kappa4(1,5)*(*states[5])[iEq]);

   residual[2][iEq] += volume*(kappa4(2,0)*(*states[0])[iEq] + kappa4(2,1)*(*states[1])[iEq] +
                          kappa4(2,2)*(*states[2])[iEq] + kappa4(2,3)*(*states[3])[iEq] +
                          kappa4(2,4)*(*states[4])[iEq] + kappa4(2,5)*(*states[5])[iEq]);

   residual[3][iEq] += 0.5*subresidual[2][iEq] +
                                  volume*(kappa4(3,0)*(*states[0])[iEq] + kappab4(2,0)*temp_vect0[iEq] +
                                          kappa4(3,1)*(*states[1])[iEq] + kappab4(2,1)*temp_vect1[iEq] +
                                          kappa4(3,2)*(*states[2])[iEq] + kappab4(2,2)*temp_vect2[iEq] +
                                          kappa4(3,3)*(*states[3])[iEq] + kappab4(2,3)*temp_vect3[iEq] +
                                          kappa4(3,4)*(*states[4])[iEq] + kappab4(2,4)*temp_vect4[iEq] +
                                          kappa4(3,5)*(*states[5])[iEq] + kappab4(2,5)*temp_vect5[iEq] ) ;
}




temp_vect0 = betasInTriag[0]*(*states[0]);
temp_vect1 = betasInTriag[0]*(*states[1]);
temp_vect2 = betasInTriag[0]*(*states[2]);
temp_vect3 = betasInTriag[0]*(*states[3]);
temp_vect4 = betasInTriag[0]*(*states[4]);
temp_vect5 = betasInTriag[0]*(*states[5]);

 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[4][iEq] += 0.5*subresidual[0][iEq] +
                                  volume*(kappa4(4,0)*(*states[0])[iEq] + kappab4(0,0)*temp_vect0[iEq]  +
                                          kappa4(4,1)*(*states[1])[iEq] + kappab4(0,1)*temp_vect1[iEq]  +
                                          kappa4(4,2)*(*states[2])[iEq] + kappab4(0,2)*temp_vect2[iEq]  +
                                          kappa4(4,3)*(*states[3])[iEq] + kappab4(0,3)*temp_vect3[iEq]  +
                                          kappa4(4,4)*(*states[4])[iEq] + kappab4(0,4)*temp_vect4[iEq]  +
                                          kappa4(4,5)*(*states[5])[iEq] + kappab4(0,5)*temp_vect5[iEq] ) ;
}




temp_vect0 = betasInTriag[1]*(*states[0]);
temp_vect1 = betasInTriag[1]*(*states[1]);
temp_vect2 = betasInTriag[1]*(*states[2]);
temp_vect3 = betasInTriag[1]*(*states[3]);
temp_vect4 = betasInTriag[1]*(*states[4]);
temp_vect5 = betasInTriag[1]*(*states[5]);

 for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
   residual[5][iEq] += 0.5*subresidual[1][iEq] +
                                  volume*(kappa4(5,0)*(*states[0])[iEq] + kappab4(1,0)*temp_vect0[iEq] +
                                          kappa4(5,1)*(*states[1])[iEq] + kappab4(1,1)*temp_vect1[iEq] +
                                          kappa4(5,2)*(*states[2])[iEq] + kappab4(1,2)*temp_vect2[iEq] +
                                          kappa4(5,3)*(*states[3])[iEq] + kappab4(1,3)*temp_vect3[iEq] +
                                          kappa4(5,4)*(*states[4])[iEq] + kappab4(1,4)*temp_vect4[iEq] +
                                          kappa4(5,5)*(*states[5])[iEq] + kappab4(1,5)*temp_vect5[iEq]) ;



             }




cellnormals.unscale();
}

//////////////////////////////////////////////////////////////////////////////

void STM_HOCRD_SplitStrategy::computeHOFluctuation_present()
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

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0 / ddata.cell->computeVolume();

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {
    // We devide the face normal by two because we integrate on each sub-faces
    facenormal[XX] = 0.5*cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = 0.5*cellnormals.getNodalNormComp(iFace%3,YY);

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

    computeStatesData(3, m_updateVar, qdstates, m_pdata, qdExtraVars); // three quadrature points per face
      
    
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

void STM_HOCRD_SplitStrategy::computeHOFluctuation_past()
{
DistributionData& ddata = getMethodData().getDistributionData();
  cf_assert(qdstates.size() == 3); // only triags so three quadrature points per face

  const State& state0 = *(_pastStates[0]);
  const State& state1 = *(_pastStates[1]);
  const State& state2 = *(_pastStates[2]);
  const State& state3 = *(_pastStates[3]);
  const State& state4 = *(_pastStates[4]);
  const State& state5 = *(_pastStates[5]);

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

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0 / ddata.cell->computeVolume();

 // Computation  of the fluctuation of each faces of the cell
for (CFuint iFace = 0; iFace < subfacetable.nbRows(); ++iFace)
  {    // We devide the face normal by two because we integrate on each sub-faces
    facenormal[XX] = 0.5*cellnormals.getNodalNormComp(iFace%3,XX);
    facenormal[YY] = 0.5*cellnormals.getNodalNormComp(iFace%3,YY);

    for (CFuint iQd = 0; iQd < 3; ++iQd)
    {

      const Node& node0 = _pastStates[subfacetable(iFace,0)]->getCoordinates();
      const Node& node1 = _pastStates[subfacetable(iFace,1)]->getCoordinates();
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

    computeStatesData(3, m_updateVar, qdstates, m_pdata, qdExtraVars); // three quadrature points per face
      
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

void STM_HOCRD_SplitStrategy::unsetup()
{
  for (CFuint i = 0; i < _pastStates.size(); ++i) {
    deletePtr(_pastStates[i]);
  }

  // AL: this could be changed ...
  // where are deleted the qdstates ???
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    deletePtr(m_qdExtraVars[i]);
  }

  for (CFuint i = 0; i < qdExtraVars.size(); ++i) {
    deletePtr(qdExtraVars[i]);
  }

  FluctuationSplitStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STM_HOCRD_SplitStrategy::setup()
{
  CFAUTOTRACE;

  FluctuationSplitStrategy::setup();

  _null.resize(PhysicalModelStack::getActive()->getDim());
  _null = 0.;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();
  m_splitter->setUpdateCoeff(updateCoeff);

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  _pastStates.resize(maxNbStatesInCell);
  // Resizing pastStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    _pastStates[i] = new State();
   }

  // Resizing m_flux_past, m_flux, m_phi
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_flux_past_space.resize(nbEqs);
  m_flux_space.resize(nbEqs);
  m_flux_past_time.resize(nbEqs);
  m_flux_time.resize(nbEqs);
  m_phi.resize(nbEqs);
 (getMethodData().getDistributionData().phi_time).resize(nbEqs);
 // AL: check if this is effective
  getMethodData().getUpdateVar()->setExtraData(true);

  // now complement it
  m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

  // set up the contour integrator
  m_contourIntegrator = getMethodData().getContourIntegrator();
  m_contourIntegrator->setFaceNormals(&m_unitFaceNormals);

  const CFuint maxNbQPts = m_contourIntegrator->getMaxIntegratorPattern().totalNbPts();
  // physical data evaluated in the quadrature points
  m_pdata.resize(maxNbQPts);
  for (CFuint  i = 0; i < maxNbQPts; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
  m_contourIntegrator->getValues(m_qdstates);

  const CFuint extra_var_size = getMethodData().getUpdateVar()->getExtraPhysicalVarsSize();
  m_qdExtraVars.resize(m_qdstates.size());
  for (CFuint i = 0; i < m_qdExtraVars.size(); ++i) {
    m_qdExtraVars[i] = new RealVector(extra_var_size);
  }

  SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbGeoEnts = innerCells->getLocalNbGeoEnts();
  m_nbQPointsInCell.resize(nbGeoEnts);
   CFuint nbStatesInCell;
  // back up cell state pointers
  m_statesBkp.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell());

  // Resizing the storage for the full past Residuals
  _pastResiduals_order1.resize(nbGeoEnts);

  // Resizing the storage for the full past Residuals
  _pastResiduals.resize(nbGeoEnts);
  _stdTrsGeoBuilder.getDataGE().trs = innerCells;

  // loop over all the cells and set the number of quadrature
  // points in the m_nbQPointsInCell vector
  for (CFuint iGeoEnt = 0; iGeoEnt < nbGeoEnts; ++iGeoEnt) {

    // build the GeometricEntity
    _stdTrsGeoBuilder.getDataGE().idx = iGeoEnt;
    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();
    // CFout << "m_nbQPointsInCell["<<iGeoEnt<<"]: "<< m_nbQPointsInCell[iGeoEnt]<<"\n";
    // CFout << "cell["<<iGeoEnt<<"].nbNodes: "<< cell.nbNodes()<<"\n";
    m_contourIntegrator->setNbSolQuadraturePoints(&cell, m_nbQPointsInCell[iGeoEnt]);

    nbStatesInCell = cell.nbStates();
    _pastResiduals_order1[iGeoEnt].resize(nbStatesInCell*nbEqs);
   _pastResiduals[iGeoEnt].resize(nbStatesInCell*nbEqs);
    // release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
  }

  temp_residual.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
    temp_residual[iState].resize(nbEqs);

 m_unitFaceNormals.resize(MeshDataStack::getActive()->Statistics().getMaxNbFacesInCell());
  for (CFuint  i = 0; i < m_unitFaceNormals.size(); ++i) {
    m_unitFaceNormals[i].resize(PhysicalModelStack::getActive()->getDim());
  }

   // number of quadrature point used to compute the fluctuatuion
  const CFuint nbQdPts = 3;
  
  // sub elemt table
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element

  subresidual[0].resize(nbEqs);
  subresidual[1].resize(nbEqs);
  subresidual[2].resize(nbEqs);

  m_phisubT.resize(4); // P2 triangle has 4 sub triangles

  m_phisubT[0] = new RealVector(nbEqs);
  m_phisubT[1] = new RealVector(nbEqs);
  m_phisubT[2] = new RealVector(nbEqs);
  m_phisubT[3] = new RealVector(nbEqs);


  // sub elemt table : contain the faces of each sub-element
  subelemtable.resize(4,3); // 4 sub elems with 3 faces each
  subelemfacedir.resize(4,3);

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

  qdExtraVars.resize(3); // 3 quadrature points per face
  qdExtraVars[0] = new RealVector();
  qdExtraVars[1] = new RealVector();
  qdExtraVars[2] = new RealVector();

  qdnodes.resize(3); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);


  facenormal.resize(2); // only 2D

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

  kappa1.resize(6,6);
  kappa2.resize(6,6);
  kappa3.resize(6,6);
  kappa4.resize(6,6);

  CFreal coeff1 = 29./40320.;
  CFreal coeff2 = 1./20160.;
  CFreal coeff3 = 1./4480.;
  CFreal coeff4 = -1./8064.;
  CFreal coeff5 = (-1./1344.);
  CFreal coeff6 = (1./10080.);
  CFreal coeff7 = 1./2520.;
  CFreal coeff8 = (-1./504.);
  CFreal coeff9 = (-1./4032.);
  CFreal coeff10 = (41./6720.);
  CFreal coeff11 = (-1./336.);
  CFreal coeff12 = (-1./630.);
  CFreal coeff13 = 1./3360.;
  CFreal coeff14 = 1./1260.;


  kappa1(0,0) = coeff10;
  kappa1(0,1) = coeff1;
  kappa1(0,2) = coeff1;
  kappa1(0,3) = coeff11;
  kappa1(0,4) = coeff12;
  kappa1(0,5) = coeff11;

  kappa1(1,0) = coeff3;
  kappa1(1,1) = coeff2;
  kappa1(1,2) = coeff4;
  kappa1(1,3) = coeff5;
  kappa1(1,4) = -coeff2;
  kappa1(1,5) = 13.*coeff2;

  kappa1(2,0) = coeff3;
  kappa1(2,1) = coeff4;
  kappa1(2,2) = coeff2;
  kappa1(2,3) = 13.0*coeff2;
  kappa1(2,4) = -coeff2;
  kappa1(2,5) = coeff5;

  kappa1(3,0) = coeff8;
  kappa1(3,1) = coeff9;
  kappa1(3,2) = 23.*coeff2;
  kappa1(3,3) = 83.*coeff6;
  kappa1(3,4) = -coeff6;
  kappa1(3,5) = -71.*coeff6;

  kappa1(4,0) = -11.*coeff6;
  kappa1(4,1) = -coeff2;
  kappa1(4,2) = -coeff2;
  kappa1(4,3) = coeff7;
  kappa1(4,4) = coeff7;
  kappa1(4,5) = coeff7;

  kappa1(5,0) = coeff8;
  kappa1(5,1) = 23.*coeff2;
  kappa1(5,2) = coeff9;
  kappa1(5,3) = -71.*coeff6;
  kappa1(5,4) = -coeff6;
  kappa1(5,5) = 83.*coeff6;



  kappa2(0,0) = coeff2;
  kappa2(0,1) = coeff3;
  kappa2(0,2) = coeff4;
  kappa2(0,3) = coeff5;
  kappa2(0,4) = 13.*coeff2;
  kappa2(0,5) = -coeff2;

  kappa2(1,0) = coeff1;
  kappa2(1,1) = coeff10;
  kappa2(1,2) = coeff1;
  kappa2(1,3) = coeff11;
  kappa2(1,4) = coeff11;
  kappa2(1,5) = coeff12;

  kappa2(2,0) = coeff4;
  kappa2(2,1) = coeff3;
  kappa2(2,2) = coeff2;
  kappa2(2,3) = 13.*coeff2;
  kappa2(2,4) = coeff5;
  kappa2(2,5) = -coeff2;

  kappa2(3,0) = coeff9;
  kappa2(3,1) = coeff8;
  kappa2(3,2) = 23.0*coeff2;
  kappa2(3,3) = 83.0*coeff6;
  kappa2(3,4) = -71.0*coeff6;
  kappa2(3,5) = -coeff6;

  kappa2(4,0) = 23.0*coeff2;
  kappa2(4,1) = coeff8;
  kappa2(4,2) = coeff9;
  kappa2(4,3) = -71.0*coeff6;
  kappa2(4,4) = 83.0*coeff6;
  kappa2(4,5) = -coeff6;

  kappa2(5,0) = -coeff2;
  kappa2(5,1) = -11.0*coeff6;
  kappa2(5,2) = -coeff2;
  kappa2(5,3) = coeff7;
  kappa2(5,4) = coeff7;
  kappa2(5,5) = coeff7;



  kappa3(0,0) = coeff2;
  kappa3(0,1) = coeff4;
  kappa3(0,2) = coeff3;
  kappa3(0,3) = -coeff2;
  kappa3(0,4) = 13.0*coeff2;
  kappa3(0,5) = coeff5;

  kappa3(1,0) = coeff4;
  kappa3(1,1) = coeff2;
  kappa3(1,2) = coeff3;
  kappa3(1,3) = -coeff2;
  kappa3(1,4) = coeff5;
  kappa3(1,5) = 13.0*coeff2;

  kappa3(2,0) = coeff1;
  kappa3(2,1) = coeff1;
  kappa3(2,2) = coeff10;
  kappa3(2,3) = coeff12;
  kappa3(2,4) = coeff11;
  kappa3(2,5) = coeff11;

  kappa3(3,0) = -coeff2;
  kappa3(3,1) = -coeff2;
  kappa3(3,2) = -11.*coeff6;
  kappa3(3,3) = coeff7;
  kappa3(3,4) = coeff7;
  kappa3(3,5) = coeff7;

  kappa3(4,0) = 23.0*coeff2;
  kappa3(4,1) = coeff9;
  kappa3(4,2) = coeff8;
  kappa3(4,3) = -coeff6;
  kappa3(4,4) = 83.*coeff6;
  kappa3(4,5) = -71.*coeff6;

  kappa3(5,0) = coeff9;
  kappa3(5,1) = 23.*coeff2;
  kappa3(5,2) = coeff8;
  kappa3(5,3) = -coeff6;
  kappa3(5,4) = -71.*coeff6;
  kappa3(5,5) = 83.*coeff6;



  kappa4(0,0) = coeff2;
  kappa4(0,1) = coeff4;
  kappa4(0,2) = coeff4;
  kappa4(0,3) = coeff13;
  kappa4(0,4) = -coeff7;
  kappa4(0,5) = coeff13;

  kappa4(1,0) = coeff4;
  kappa4(1,1) = coeff2;
  kappa4(1,2) = coeff4;
  kappa4(1,3) = coeff13;
  kappa4(1,4) = coeff13;
  kappa4(1,5) = -coeff7;

  kappa4(2,0) = coeff4;
  kappa4(2,1) = coeff4;
  kappa4(2,2) = coeff2;
  kappa4(2,3) = -coeff7;
  kappa4(2,4) = coeff13;
  kappa4(2,5) = coeff13;

  kappa4(3,0) = coeff14;
  kappa4(3,1) = coeff14;
  kappa4(3,2) = coeff6;
  kappa4(3,3) = 41.*coeff6;
  kappa4(3,4) = -29.0*coeff6;
  kappa4(3,5) = -29.0*coeff6;

  kappa4(4,0) = coeff6;
  kappa4(4,1) = coeff14;
  kappa4(4,2) = coeff14;
  kappa4(4,3) = -29.0*coeff6;
  kappa4(4,4) = 41.0*coeff6;
  kappa4(4,5) = -29.0*coeff6;

  kappa4(5,0) = coeff14;
  kappa4(5,1) = coeff6;
  kappa4(5,2) = coeff14;
  kappa4(5,3) = -29.*coeff6;
  kappa4(5,4) = -29.0*coeff6;
  kappa4(5,5) = 41.0*coeff6;

  kappab1.resize(3,6);
  kappab2.resize(3,6);
  kappab3.resize(3,6);
  kappab4.resize(3,6);

  coeff1 = 1./84.;
  coeff2 = 3./56.;
  coeff3 = 5./168.;

  kappab1(0,0) = coeff3;
  kappab1(0,1) = -coeff1;
  kappab1(0,2) = -coeff1;
  kappab1(0,3) = coeff2;
  kappab1(0,4) = coeff1;
  kappab1(0,5) = coeff2;

  kappab1(1,0) = coeff3;
  kappab1(1,1) = -coeff1;
  kappab1(1,2) = -coeff1;
  kappab1(1,3) = coeff2;
  kappab1(1,4) = coeff1;
  kappab1(1,5) = coeff2;

  kappab1(2,0) = coeff3;
  kappab1(2,1) = -coeff1;
  kappab1(2,2) = -coeff1;
  kappab1(2,3) = coeff2;
  kappab1(2,4) = coeff1;
  kappab1(2,5) = coeff2;



  kappab2(0,0) = -coeff1;
  kappab2(0,1) = coeff3;
  kappab2(0,2) = -coeff1;
  kappab2(0,3) = coeff2;
  kappab2(0,4) = coeff2;
  kappab2(0,5) = coeff1;

  kappab2(1,0) = -coeff1;
  kappab2(1,1) = coeff3;
  kappab2(1,2) = -coeff1;
  kappab2(1,3) = coeff2;
  kappab2(1,4) = coeff2;
  kappab2(1,5) = coeff1;

  kappab2(2,0) = -coeff1;
  kappab2(2,1) = coeff3;
  kappab2(2,2) = -coeff1;
  kappab2(2,3) = coeff2;
  kappab2(2,4) = coeff2;
  kappab2(2,5) = coeff1;



  kappab3(0,0) = -coeff1;
  kappab3(0,1) = -coeff1;
  kappab3(0,2) = coeff3;
  kappab3(0,3) = coeff1;
  kappab3(0,4) = coeff2;
  kappab3(0,5) = coeff2;

  kappab3(1,0) = -coeff1;
  kappab3(1,1) = -coeff1;
  kappab3(1,2) = coeff3;
  kappab3(1,3) = coeff1;
  kappab3(1,4) = coeff2;
  kappab3(1,5) = coeff2;

  kappab3(2,0) = -coeff1;
  kappab3(2,1) = -coeff1;
  kappab3(2,2) = coeff3;
  kappab3(2,3) = coeff1;
  kappab3(2,4) = coeff2;
  kappab3(2,5) = coeff2;



  kappab4(0,0) = -coeff1;
  kappab4(0,1) = -coeff1;
  kappab4(0,2) = -coeff1;
  kappab4(0,3) = coeff2;
  kappab4(0,4) = coeff2;
  kappab4(0,5) = coeff2;

  kappab4(1,0) = -coeff1;
  kappab4(1,1) = -coeff1;
  kappab4(1,2) = -coeff1;
  kappab4(1,3) = coeff2;
  kappab4(1,4) = coeff2;
  kappab4(1,5) = coeff2;

  kappab4(2,0) = -coeff1;
  kappab4(2,1) = -coeff1;
  kappab4(2,2) = -coeff1;
  kappab4(2,3) = coeff2;
  kappab4(2,4) = coeff2;
  kappab4(2,5) = coeff2;


  // tell the splitters to compute the betas
  getMethodData().getDistributionData().computeBetas = true;

  m_updateVar   = getMethodData().getUpdateVar();

temp_residual.resize(maxNbStatesInCell);

  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
    temp_residual[iState].resize(nbEqs);

  getMethodData().getDistributionData().isHO = true;

temp_vect0.resize(nbEqs);
    temp_vect1.resize(nbEqs);
    temp_vect2.resize(nbEqs);
    temp_vect3.resize(nbEqs);
    temp_vect4.resize(nbEqs);
    temp_vect5.resize(nbEqs);


}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
STM_HOCRD_SplitStrategy::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();

   result.push_back(&socket_pastNormals);
   result.push_back(&socket_pastStates);
   result.push_back(&socket_pastCellVolume);
   result.push_back(&socket_cellSpeed);
   result.push_back(&socket_updateCoeff);

   return result;
}
//////////////////////////////////////////////////////////////////////////////

void STM_HOCRD_SplitStrategy::setCurrentCell()
{
  DistributionData& ddata = getMethodData().getDistributionData();

  // back up the update states and sets them in the linearizer
  // AL: don't touch this! update states have to be backed up
  // because linearizing states are temporarily stored
  // in the cell and used to contour integrate
  const CFuint nbCellStates = ddata.states->size();

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {

    m_statesBkp[iState] = (*ddata.states)[iState];
  }
  getMethodData().getLinearizer()->setUpdateStates(&m_statesBkp);

  ddata.tStates = computeConsistentStates(ddata.states);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
