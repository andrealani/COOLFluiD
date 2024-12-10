#include "Framework/MethodStrategyProvider.hh"

#include "HyperPoisson/HyperPoisson3DVarSet.hh"

#include "FluxReconstructionHyperPoisson/FluxReconstructionHyperPoisson.hh"
#include "FluxReconstructionHyperPoisson/BCOutletHyperPoisson.hh"

#include "Common/NotImplementedException.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::HyperPoisson;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCOutletHyperPoisson,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionHyperPoissonModule >
  BCOutletHyperPoissonProvider("OutletHyperPoisson");

//////////////////////////////////////////////////////////////////////////////

BCOutletHyperPoisson::BCOutletHyperPoisson(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_refPhi = 0.0;
  setParameter("refPhi",&m_refPhi);
  
  m_dipole = false;
  setParameter("isDipole",&m_dipole);
}

//////////////////////////////////////////////////////////////////////////////

BCOutletHyperPoisson::~BCOutletHyperPoisson()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCOutletHyperPoisson::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the outlet.");
  options.addConfigOption< bool >("isDipole","Prescribe analytical dipole phi.");
}

//////////////////////////////////////////////////////////////////////////////

void BCOutletHyperPoisson::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  //const CFuint nbrStates = ghostStates.size();
  //cf_assert(nbrStates == intStates.size());
  //cf_assert(nbrStates == normals.size());

  CFuint nbrStates = m_nbrFaceFlxPnts;
  
  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {              
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);
    
    for (CFuint iEq = 1; iEq < (intState.size()); ++iEq)
    {
      ghostState[iEq] = intState[iEq];
    }
    const CFreal x = coords[iState][XX];
    const CFreal y = coords[iState][YY];
    const CFreal z = coords[iState][ZZ];
    const CFreal r = (coords[iState]).norm2();
    if (m_dipole)
    {
      m_refPhi=-1.0/3.0*0.666*z/(r*r*r);
    }
    //phi
    ghostState[0] = (2.*m_refPhi)-intState[0];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCOutletHyperPoisson::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  //const CFuint nbrStateGrads = intGrads.size();
  //cf_assert(nbrStateGrads == ghostGrads.size());
  //cf_assert(nbrStateGrads == normals.size());

  
}

//////////////////////////////////////////////////////////////////////////////

void BCOutletHyperPoisson::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = true;
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  m_nbrFaceFlxPnts = frLocalData[0]->getFaceFlxPntsFaceLocalCoords()->size();
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  CFuint nbrFaceFlxPntsMax= (order+1)*(order+1);

  //m_tempStates.resize(nbrFaceFlxPnts);

  m_tempStates.resize(nbrFaceFlxPntsMax);
  
  // number of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  //for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPntsMax; ++iFlx)
  {
    m_tempStates[iFlx].resize(nbEqs); 
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

