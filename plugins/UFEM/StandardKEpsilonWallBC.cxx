#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/StandardKEpsilonWallBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StandardKEpsilonWallBC,UFEMSolverData,UFEMPlugin > StandardKEpsilonWallBCProvider("StandardKEpsilonWallBC");

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilonWallBC::defineConfigOptions(Config::OptionList& options)
{
   CFAUTOTRACE;
   DirichletBC::defineConfigOptions(options);

   options.addConfigOption< CFreal >("Cmu",      "Cmu");
   options.addConfigOption< CFreal >("K",        "Karman constant");
   options.addConfigOption< CFreal >("RhoElm",   "Density");
   options.addConfigOption< CFreal >("MuLam",    "Laminar viscosity");
}

//////////////////////////////////////////////////////////////////////////////

StandardKEpsilonWallBC::StandardKEpsilonWallBC(const std::string& name) :
  DirichletBC(name),
  socket_connBndFace2InnerCell("connBndFace2InnerCell"),
  socket_connBndFace2BndState("connBndFace2BndState")
{
  CFAUTOTRACE;

  m_Cmu = 90.; //0.09
  m_K = 90.; //0.41
  m_RhoElm = 90.; //1
  m_MuLam = 90.; //0.0001

  setParameter( "Cmu",                   &m_Cmu);
  setParameter( "K",                     &m_K);
  setParameter( "RhoElm",                &m_RhoElm);
  setParameter( "MuLam",                 &m_MuLam);

}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilonWallBC::setup()
{
  CFAUTOTRACE;
  DirichletBC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilonWallBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  DirichletBC::configure(args);

  CFLog(INFO, getClassName() << ": Cmu: "            << m_Cmu         << "\n");
  CFLog(INFO, getClassName() << ": K: "              << m_K           << "\n");
  CFLog(INFO, getClassName() << ": RhoElm: "         << m_RhoElm      << "\n");
  CFLog(INFO, getClassName() << ": MuLam: "          << m_MuLam       << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilonWallBC::executeOnTrs()
{
  CFAUTOTRACE;
  DirichletBC::executeOnTrs();
/*
  // finding both trs
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<TopologicalRegionSet> innerTrs=0;
  CFuint innerTrsID=0;
  Common::SafePtr<TopologicalRegionSet> bndTrs = 0;
  // THIS DOES NOT WORK PROPERLY, MAYBE WRONG USE? //const CFuint bndTrsID = getCurrentTrsID();
  CFuint bndTrsID=0;
  std::string bndTrsName=getCurrentTRS()->getName();
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
      innerTrsID = i;
    }
    if (trs[i]->getName() == bndTrsName) {
      bndTrs = trs[i];
      bndTrsID = i;
     }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs not found with name: InnerCells");
  if (bndTrs==0) Common::NoSuchValueException(FromHere(),"Trs not found with name: "+bndTrsName);

  // trs and geo props
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > bndGeoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& bndGeoData = bndGeoBuilder->getDataGE();
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> > innerGeoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& innerGeoData = innerGeoBuilder->getDataGE();

  // trs-adjacent innercell and face 2 state connectivity
  DataHandle< std::vector<CFuint> > connBndFace2InnerCell = socket_connBndFace2InnerCell.getDataHandle();
  DataHandle< Common::ConnectivityTable<CFuint> > connBndFace2BndState = socket_connBndFace2BndState.getDataHandle();

  // interstates
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();

  // get/allocate common used variables
  const CFuint nbStates = bndTrs->getNbStatesInTrs();
  const CFuint nbGeos   = bndTrs->getLocalNbGeoEnts();
  const CFuint nbEqs    = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim    = PhysicalModelStack::getActive()->getDim();

  // data from innercells, encoding of wall adjacent cell average values: 0..nbEqs-1,sumarea,walldist
  RealMatrix bndProps(nbStates,nbEqs+2);
  RealVector innerAvgs(nbEqs+1+nbDim);

  // data to store state-wise setting
  MathTools::CFVec<bool> applyEqs(nbEqs);
  MathTools::CFVec<CFreal> applyVars(nbEqs);

  // calculating innercell averages and summing to nodes
  for (CFuint igeo=0; igeo<nbGeos; igeo++){

    // inner's geoent props
    innerGeoData.trs = innerTrs;
    innerGeoData.idx = connBndFace2InnerCell[bndTrsID][igeo];
    GeometricEntity *const innerGeo = bndGeoBuilder->buildGE();

    // innercell averages
    std::vector<State*>* innerStates = innerGeo->getStates();
    const CFuint innerNbGeoStates=innerStates->size();
    innerAvgs=0.;
    CFreal innerVol=innerGeo->computeVolume();
    CFreal invNbStates= 1./((CFreal)innerNbGeoStates);
    for (CFuint iState=0; iState<innerNbGeoStates; ++iState) {
      State *const state = (*innerStates)[iState];
      State& interState = *interStates[state->getLocalID()];
      for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
        innerAvgs[iEq] += interState[iEq];
      }
      Node innerNod=state->getCoordinates();
      for (CFuint idim=0; idim<nbDim; ++idim)
        innerAvgs[nbEqs+1+idim]+=innerNod[idim];
    }
    innerAvgs *= invNbStates;
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) innerAvgs[iEq] *= innerVol;
    innerAvgs[nbEqs]=innerVol;

    // release geometric entity
    innerGeoBuilder->releaseGE();

    // switching to actual trs
    bndGeoData.trs = bndTrs;
    bndGeoData.idx = igeo;
    GeometricEntity *const bndGeo = bndGeoBuilder->buildGE();

    // calculating wall distance
    RealVector bndNormal=bndGeo->computeAvgCellNormal();
    bndNormal.normalize();
    Node* bndZeroNode=bndGeo->getNode(0);
    CFreal wdist=0.;
    for (CFuint idim=0; idim<nbDim; ++idim)
      wdist+=(bndNormal[idim]*((*bndZeroNode)[idim]-innerAvgs[nbEqs+1+idim]))*innerVol;

    // looping on current trs faces and summing to states
    //CFuint bndNbGeoStates=bndGeo->getStates();
    std::vector< State* >* bndGeoStates=bndGeo->getStates();
    const CFuint bndNbGeoStates= bndGeoStates->size();
    for (CFuint istate=0; istate<bndNbGeoStates; istate++){
      for (CFuint ieq=0; ieq<(const CFuint)(nbEqs+1); ieq++)
        bndProps(connBndFace2BndState[bndTrsID](igeo,istate),ieq)+=innerAvgs[ieq];
      bndProps(connBndFace2BndState[bndTrsID](igeo,istate),nbEqs+1)+=wdist;
    }

    // release geometric entity
    bndGeoBuilder->releaseGE();

  }

  // now comes the loop over the states
  for (CFuint istate=0; istate<nbStates; istate++){

    const CFuint nLocalID = (*bndTrs->getStatesInTrs())[istate] ;
    const State *bndState = states[nLocalID];

    if (bndState->isParUpdatable()) {

      // first dividing by sum volume to get the average
      CFreal invSumVol=1./bndProps(istate,nbEqs);
      for (CFuint ieq=0; ieq<nbEqs; ieq++)
        bndProps(istate,ieq)*=invSumVol;
      const CFreal walldist=bndProps(istate,nbEqs+1)*=invSumVol;

      // init
      applyEqs=false;
      applyVars=0.;

      // Defining yplus
      double yplus = (m_RhoElm*pow(m_Cmu,0.25)*walldist*sqrt(bndProps(istate,nbDim+1)))/m_MuLam;

      // velocity always zero
      for(CFuint iDim=0; iDim<nbDim; ++iDim) {
         applyEqs[iDim+1] = true;
         applyVars[iDim+1] = 0.;
      }

      // Laminar region
      if (yplus<=11.225) {

const CFreal klimit= 1.e-10;
const CFreal epslimit= m_Cmu*klimit*klimit/(m_MuLam*1.e-6);

        applyEqs[nbDim+1] = true;
        applyVars[nbDim+1] = klimit; //1.e-6;
        applyEqs[nbDim+2] = true;
        applyVars[nbDim+2] = epslimit; //pow(m_Cmu,3./4.)*pow(bndProps(istate,nbDim+1),3./2.)*yplus/walldist;
      } else {
        applyEqs[nbDim+2] = true;
        applyVars[nbDim+2] = pow(m_Cmu,3./4.)*pow(bndProps(istate,nbDim+1),3./2.)*(1./m_K*log(yplus))/walldist;
      }

    }
  }
*/
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilonWallBC::computeStateValuesStandardKEpsilonWallBC(const Framework::State* currState)
{
  computeStateValuesDirichletBC(currState);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > StandardKEpsilonWallBC::needsSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<BaseDataSocketSink> > result=DirichletBC::needsSockets();

  result.push_back(&socket_connBndFace2InnerCell);
  result.push_back(&socket_connBndFace2BndState);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

