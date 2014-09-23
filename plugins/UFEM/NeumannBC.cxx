#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/NeumannBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< NeumannBC,UFEMSolverData,UFEMPlugin > NeumannBCProvider("NeumannBC");

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::defineConfigOptions(Config::OptionList& options)
{
   CFAUTOTRACE;

   options.addConfigOption< bool >( "Implicit", "Apply the BC implicitly. Can be false (dU/dn*k_U=Vars_U) or true (dU/dn*k_U=U*Vars_U). (default is false)" );
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::NeumannBC(const std::string& name) :
  BaseBC(name),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors")
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_isImplicit = false;

  setParameter( "Implicit",      &m_isImplicit);
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setup()
{
  CFAUTOTRACE;
  BaseBC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  BaseBC::configure(args);

  CFLog(INFO, getClassName() << ": Implicit: "       << m_isImplicit  << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::executeOnTrs()
{
  CFAUTOTRACE;
  BaseBC::executeOnTrs();

#define kronecker(i,j) ((i)==(j)?1.:0.)

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle< bool > isUpdated = socket_isUpdated.getDataHandle();
  DataHandle< State* > interStates = socket_interStates.getDataHandle();
  DataHandle< State*, GLOBAL > states = socket_states.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "NeumannBC::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // PhysicalModel properties and auxiliary variables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  const CFuint dim   = PhysicalModelStack::getActive()->getDim();
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  // loop over the faces
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity *const face = geoBuilder->buildGE();

    vector<State*>* faceStates = face->getStates();

    // getting number of states
    const CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();

    // getting normal vector
    CFreal vol=(face->computeAvgCellNormal()).norm2();
    vol/=((CFreal)((dim)*(dim+1)));

    // evaluate boundary condition at space, time and states,
    // these are the dFI/dn values, needs a preliminary loop over nbStatesInFace
    RealVector states(nbStatesInFace*nbEqs);
    RealMatrix applyvars(nbStatesInFace,nbEqs);
    for (CFuint iState = 0; iState < nbStatesInFace; ++iState) {
      State *const currInterState = interStates[(*faceStates)[iState]->getLocalID()];
      State *const currState = (*faceStates)[iState];
      computeStateValuesNeumannBC(currInterState);
      for (CFuint i=0; i<nbEqs; ++i) {
        states[iState*nbEqs+i] = (*currState)[i];
        applyvars(iState,i)=m_applyVars[i];
      }
    }

    // the element local matrix
    RealMatrix A(nbStatesInFace*nbEqs,nbStatesInFace*nbEqs);
    RealVector b(nbStatesInFace*nbEqs);
    RealMatrix eK(nbStatesInFace*nbEqs,nbStatesInFace*nbEqs);
    RealVector eR(nbStatesInFace*nbEqs);
    A=0.;
    b=0.;

    if (m_isImplicit) {

      // assembly
      for(CFuint iEq=0; iEq<nbEqs; ++iEq) if (m_applyFlags[0][iEq]){
        for (CFuint i=0; i<nbStatesInFace; i++) {
          State *const currState = (*faceStates)[i];
            if (!isUpdated[currState->getLocalID()*nbEqs+iEq]) for (CFuint j=0; j<nbStatesInFace; j++){
              A(i*nbEqs+iEq,j*nbEqs+iEq)-=vol*(1.+kronecker(i,j))*applyvars(j,iEq);
            }  // i
        } // j
      } // iEq

    } else {

      // assembly
      for(CFuint iEq=0; iEq<nbEqs; ++iEq) if (m_applyFlags[0][iEq]) {
        for (CFuint i=0; i<nbStatesInFace; i++) {
          State *const currState = (*faceStates)[i];
          if (!isUpdated[currState->getLocalID()*nbEqs+iEq]) for (CFuint j=0; j<nbStatesInFace; j++){
            b[i*nbEqs+iEq]-=vol*(1.+kronecker(i,j))*applyvars(j,iEq);
          }  // i
        } // j
      } // iEq

    }

    // time split
    eK=A;
    eR=A*states+b;

    // blockacvcumulator
    BlockAccumulator *acc = getMethodData().getLinearSystemSolver()[0]->createBlockAccumulator(nbStatesInFace,nbStatesInFace,nbEqs);

    // and system matrix
    SafePtr<LSSMatrix> jacobMatrix = getMethodData().getLinearSystemSolver()[0]->getMatrix();

    // adding to system matrix
    for (CFuint iState=0; iState<nbStatesInFace; ++iState)
    {
      CFuint nLocalID=((*faceStates)[iState])->getLocalID();
      // add the contribution to the RHS vector
      for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
        rhs(nLocalID, iEq, nbEqs) -= eR[iState*nbEqs+iEq];
      }

      // set the index of the block corresponding to the current
      // state in the jacobian matrix
      acc->setRowColIndex(iState, nLocalID);
    }

    // add the contribution to the jacobian matrix
    acc->setValuesM(eK);
    jacobMatrix->addValues( *acc );

    //release the GeometricEntity
    geoBuilder->releaseGE();

  } // iGeoEnt loop

/*
  // loop over the faces, fixing isUpdated, now its not set because dirichlet carefully separated but neumanns can sum their contributions to corner nodes
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {
    // getting geoent
    geoData.idx = iGeoEnt;
    GeometricEntity *const face = geoBuilder->buildGE();
    // getting 
    vector<State*>* faceStates = face->getStates();
    // getting number of states
    const CFuint nbStatesInFace = face->getNbNodesSolutionShapeFunction();
    // loop over equations
    for(CFuint eqi=0; eqi<nbApplyEqs; ++eqi) {
      // no isupdated fix, because neumanns can be mixed for a corner node
      for (CFuint i=0; i<nbStatesInFace; i++) {
        isUpdated[((*faceStates)[i])->getLocalID()*nbEqs+m_applyEqs[eqi]]=true;
      }
    } // equation loop 
    //release the GeometricEntity
    geoBuilder->releaseGE();
  } // iGeoEnt loop
*/

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::computeStateValuesNeumannBC(const Framework::State* currState)
{
  computeStateValuesBaseBC(currState);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > NeumannBC::needsSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<BaseDataSocketSink> > result=BaseBC::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

