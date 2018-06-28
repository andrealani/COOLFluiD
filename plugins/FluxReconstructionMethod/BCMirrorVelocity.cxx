#include "Framework/MethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCMirrorVelocity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCMirrorVelocity,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionModule >
  BCMirrovVelocityProvider("MirrorVelocity");

//////////////////////////////////////////////////////////////////////////////

BCMirrorVelocity::BCMirrorVelocity(const std::string& name) :
  BCStateComputer(name),
  m_isVelocityComp(),
  m_tangent(),
  m_velocityNGradI(),
  m_velocityTGradI(),
  m_velocityNGradG(),
  m_velocityTGradG()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCMirrorVelocity::~BCMirrorVelocity()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorVelocity::computeGhostStates(const vector< State* >& intStates,
					  vector< State* >& ghostStates,
					  const std::vector< RealVector >& normals,
					  const std::vector< RealVector >& coords)
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  CFreal vn = 0.;
  CFreal area2 = 0.;
  
  cf_assert(m_velocityIDs.size() == dim);
  for (CFuint iState = 0; iState < intStates.size() ; ++iState){
   for (CFuint i = 0; i < m_velocityIDs.size(); ++i) { 
     const CFreal nxComp = normals[iState][i];
     vn += nxComp*(*intStates[iState])[m_velocityIDs[i]];
     area2 += nxComp*nxComp;
   }
   
   CFuint jxx = 0;
   for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
     if (!m_isVelocityComp[i]) {
       (*ghostStates[iState])[i] = (*intStates[iState])[i];
     }
     else {
       (*ghostStates[iState])[i]= (*intStates[iState])[i] - 2.0*vn*normals[iState][jxx]/area2;
       jxx++;
     }
     // if (i == 2 &&  (*ghostStates[iState])[i] < -0.000001) CFLog(INFO, "intState: " << *intStates[iState] << ", ghost: " << *ghostStates[iState] << "\n");
   }
   
   //CFLog(DEBUG_MAX, "MirrorVelocity::setGhostState() => ghostState = " << *ghostState << "\n"); 
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorVelocity::computeGhostGradients
(const std::vector< std::vector< RealVector* > >& intGrads,
 std::vector< std::vector< RealVector* > >& ghostGrads,
 const std::vector< RealVector >& normals,
 const std::vector< RealVector >& coords)
{
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());
  
  // number of gradient variables
  cf_assert(nbrStateGrads > 0);
  
  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState) {
    // normal
    const RealVector& normal = normals[iState];
    
    // tangential unit vector
    m_tangent[XX] = -normal[YY];
    m_tangent[YY] =  normal[XX];
    
    const vector< RealVector* >& velocityGradI = intGrads[iState];
    m_velocityNGradI = 0.;
    m_velocityTGradI = 0.;
    CFuint jxx = 0;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (!m_isVelocityComp[i]) {
	const RealVector& varGradI =  *intGrads[iState][i];
	RealVector& varGradG =  *ghostGrads[iState][i];
	const CFreal nVarGrad = MathTools::MathFunctions::innerProd(varGradI, normal);
	varGradG = varGradI - 2.0*nVarGrad*normal;
      }
      else {
	// internal normal and tangential component
	m_velocityNGradI +=  *velocityGradI[i]*normal[jxx];
	m_velocityTGradI +=  *velocityGradI[i]*m_tangent[jxx];
	++jxx;
      }
    }
    
    m_velocityNGradG = m_velocityNGradI;
    
    jxx = 0;
    CFreal nGradUT = 0.;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (m_isVelocityComp[i]){
        // ghost normal and tangential component
        nGradUT += m_velocityTGradI[jxx]*normal[jxx];  
	++jxx;
      }
    }
    
    m_velocityTGradG = m_velocityTGradI - 2.0*nGradUT*normal;
    
    // project onto x- and y-axis
    jxx = 0;
    for (CFuint i = 0; i < m_isVelocityComp.size(); ++i) {
      if (m_isVelocityComp[i]) {
        *ghostGrads[iState][i] = m_velocityNGradG*normal[jxx] + m_velocityTGradG*m_tangent[jxx];
	++jxx;
      }
    }
  }
}



//////////////////////////////////////////////////////////////////////////////

void BCMirrorVelocity::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // set the IDs corresponding to the velocity components
  getMethodData().getUpdateVar()->setStateVelocityIDs(m_velocityIDs);
  
  if(m_velocityIDs.size() == 0) {
    CFLog(NOTICE, "MirrorVelocity::setup() => choosing default\n");
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_velocityIDs.resize(dim);
    for (CFuint i = 0 ; i < dim; ++i) {
      m_velocityIDs[i] = 1 + i; 
    }
  }
  
  m_isVelocityComp.resize(PhysicalModelStack::getActive()->getNbEq());
  m_isVelocityComp = false;
  for (CFuint i = 0 ; i < m_velocityIDs.size(); ++i) {
    m_isVelocityComp[m_velocityIDs[i]] = true;
  }

}

//////////////////////////////////////////////////////////////////////////////

void BCMirrorVelocity::unsetup()
{
  CFAUTOTRACE;

  // unsetup of the parent class
  BCStateComputer::unsetup();

  // allocate arrays to store temporary data
  const CFuint nbDims = PhysicalModelStack::getActive()->getDim();
  m_tangent.resize(nbDims);
  m_velocityNGradI.resize(nbDims);
  m_velocityTGradI.resize(nbDims);
  m_velocityNGradG.resize(nbDims);
  m_velocityTGradG.resize(nbDims);
}
    
//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

