#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

#include "FluxReconstructionTurb/FluxReconstructionSA.hh"
#include "FluxReconstructionTurb/SA2DSourceTerm.hh"

#include "SA/NavierStokesSAVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::SA;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SA2DSourceTerm, FluxReconstructionSolverData, FluxReconstructionSAModule>
SA2DSourceTermProvider("SA2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("CompressibilityCorrectionTerm","Add the extra destruction term (Default = False)");
}

//////////////////////////////////////////////////////////////////////////////

SA2DSourceTerm::SA2DSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    socket_gradients("gradients"),
    socket_wallDistance("wallDistance"),
    socket_wallShearStressVelocity("wallShearStressVelocity"),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_diffVarSet(CFNULL),
    m_solPhysData(),
    m_dummyGradients(),
    m_cellGrads(),
    m_currWallDist()
{
  addConfigOptionsTo(this);
  
  m_compTerm = false;
  setParameter("CompressibilityCorrectionTerm",&m_compTerm);
}

//////////////////////////////////////////////////////////////////////////////

SA2DSourceTerm::~SA2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  SA2DSourceTerm::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_wallShearStressVelocity);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{     
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      *(m_cellGrads[iState][iEq]) = gradients[stateID][iEq];
    }
    
    // Get the wall distance
    DataHandle< CFreal > wallDist = socket_wallDistance.getDataHandle();

    m_currWallDist[iState] = wallDist[stateID];
  }
  
  const CFreal R = m_eulerVarSet->getModel()->getR();
  
  const EquationSubSysDescriptor& eqData = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq(); 
  
  for (CFuint iSol = 0; iSol < nbrStates; ++iSol)
  {
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
    
    const CFreal dUdX = (*(m_cellGrads[iSol][1]))[XX];
    const CFreal dUdY = (*(m_cellGrads[iSol][1]))[YY];
    const CFreal dVdX = (*(m_cellGrads[iSol][2]))[XX];
    const CFreal dVdY = (*(m_cellGrads[iSol][2]))[YY];
    const CFreal dKdX = (*(m_cellGrads[iSol][4]))[XX];
    const CFreal dKdY = (*(m_cellGrads[iSol][4]))[YY];
    
    const CFreal avRhoR = m_solPhysData[EulerTerm::RHO]*R;
    const CFreal RT = R*m_solPhysData[EulerTerm::T];
    
    const CFreal dRhodX = ((*(m_cellGrads[iSol][0]))[XX] - avRhoR*(*(m_cellGrads[iSol][3]))[XX])/RT;
    const CFreal dRhodY = ((*(m_cellGrads[iSol][0]))[YY] - avRhoR*(*(m_cellGrads[iSol][3]))[YY])/RT;
    
    const CFuint iK = m_eulerVarSet->getModel()->getFirstScalarVar(0);
    
    const CFreal mu = m_diffVarSet->getLaminarDynViscosityFromGradientVars(*((*m_cellStates)[iSol]));
    const CFreal rho = m_diffVarSet->getDensity(*((*m_cellStates)[iSol]));
    const CFreal NIU = mu / rho;
    CFreal NIUtilda = (*((*m_cellStates)[iSol]))[4];
    
    //make sure we don't have negative values of niutilda
    NIUtilda = max(0.,NIUtilda);

    // To prevent division by zero (see Ashford, G., An Unstructured Grid
    // Generation and Adaptive Solution Technique for High-Reynolds-Number Compressible Flows, Ph.D. Thesis, University of Michigan 1996.)
    // Page 155
    CFreal Qsi = NIUtilda / NIU;

    //constants of the SA turbulence model:
    const CFreal Cb1 = 0.1355;
    const CFreal Cb2 = 0.622;
    const CFreal sigma = 2.0/3.0;
    const CFreal kappa = 0.41;
    const CFreal Cw1 = ( Cb1 / (kappa*kappa) ) + (( 1. + Cb2)/sigma );
    const CFreal Cw2 = 0.3;
    const CFreal Cw3 = 2.0;
    const CFreal Cv1 = 7.1;
    
    const CFreal Cv2 = 0.7;// these two parameters are used for the modified Stilda
    const CFreal Cv3 = 0.9; ///@see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
    ///of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012. 
    
    const CFreal fv1 = Qsi*Qsi*Qsi / (Qsi*Qsi*Qsi + Cv1*Cv1*Cv1);
    const CFreal fv2 = 1. - ( Qsi /(1. + (Qsi * fv1)));
    
    const CFreal d = max(m_currWallDist[iSol],1.e-10);
	 
    const CFreal Nitiloverkapa2d2 = NIUtilda/( kappa*kappa*d*d);
	 
    const CFreal S = fabs(dVdX - dUdY);// definition of the 2D vorticity magnitude
	 
    const CFreal Soverbar = Nitiloverkapa2d2*fv2;
	 
    CFreal Stilda = 0.; // definition and initialization of Stilda
	 
    //Preventing Negative Values of Modified Vorticity Stilda
    ///@see Allmaras, S. R., Johnson, F. T., and Spalart, P. R., "Modifications and Clarifications for the Implementation 
    ///of the Spalart-Allmaras Turbulence Model," ICCFD7-1902, 7th International Conference on Computational Fluid Dynamics, Big Island, Hawaii, 9-13 July 2012.
    if (Soverbar >= -(Cv2*S))
    {
      Stilda = S + Soverbar;
    }
    else
    {
      Stilda = S + (S*(Cv2*Cv2*S + Cv3*Soverbar))/((Cv3 - 2.0*Cv2)*S - Soverbar);
    }
    
    const CFreal rlim = 10.0; // definition of the first SA model
    
    CFreal r = Nitiloverkapa2d2/Stilda;
    
    r = min(r, rlim);
    
    const CFreal g = r + Cw2 * ( pow(r,6) - r);
    const CFreal g6 = g*g*g*g*g*g;
    const CFreal Cw3_6 = Cw3*Cw3*Cw3*Cw3*Cw3*Cw3;
    const CFreal sixth = 1./6.;
    const CFreal fw = g * pow( (1. + Cw3_6) / (g6 + Cw3_6) ,sixth);
    
    ///// In fully turbulent flow, ft2 doesn't need to be taken into account
    const CFreal ct3 = 1.2;
    const CFreal ct4 = 0.5;
    const CFreal ft2 = 0.0;//ct3*exp(-ct4*Qsi*Qsi);
    
    const CFreal adimCoef = m_diffVarSet->getModel().getCoeffTau();
    
    const CFreal nonConsDiffTerm = adimCoef * ( Cb2 / sigma ) * (dKdX*dKdX + dKdY*dKdY);
    const CFreal P = Cb1 * Stilda * NIUtilda * (1-ft2); // production term
    
    // correction for the introduction of rho in the convective flux : G
    const CFreal G = ( 1. / sigma ) * (NIU + NIUtilda) * ( dKdX + dKdY ) * ( dRhodX + dRhodY );
    
    CFreal D =  adimCoef * ( Cw1 * fw - Cb1 * ft2/(kappa*kappa)) * ((NIUtilda*NIUtilda) / (m_currWallDist[iSol]*m_currWallDist[iSol]));
    
    // tranform the model into SA - noft2 - comp by adding the extra destruction term
    // It improves the performance of the model in compressible mixing layers
    ///@see the site http://turbmodels.larc.nasa.gov/spalart.html#qcr2000
    if (m_compTerm)
    {
      const CFreal C5 = 3.5;
      
      const CFreal spSound = m_solPhysData[EulerTerm::A];
      
      const CFreal sum = (dUdX*dUdX + dVdX*dVdX + dUdY*dUdY + dVdY*dVdY);
      
      const CFreal extraDest = (C5*NIUtilda*NIUtilda/(spSound*spSound))*sum;
      
      D += extraDest;
    }
    
    //Compute Reynolds stress tensor 
    const CFreal positivePart = (P + nonConsDiffTerm)*rho;
    const CFreal negativePart = - (D * rho + G);
      
    /// Compute the rhs contribution
    resUpdates[m_nbrEqs*iSol + 4] = positivePart + negativePart;
    
    if (!m_isPerturbed)
    {
      DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
    
      const CFreal niuTot = NIU + NIUtilda*fv1;
      
      // take the absolute value of dUdY to avoid nan which causes tecplot to be unable to load the file
      wallShearStressVelocity[(((*m_cellStates)[iSol]))->getLocalID()] = sqrt(niuTot*fabs(dUdY));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in SA2DSourceTerm!");
  }
  cf_assert(m_nbrEqs == 5);
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
  
  m_diffVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokes2DSA>();
  cf_assert(m_diffVarSet.isNotNull());
  
  m_currWallDist.resize(m_nbrSolPnts);
  
  // size cell gradients vector
  m_cellGrads.resize(m_nbrSolPnts);
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    m_cellGrads[iState].resize(m_nbrEqs);
    
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_cellGrads[iState][iEq] = new RealVector(m_dim);
    }
  }
  
  DataHandle< CFreal > wallShearStressVelocity = socket_wallShearStressVelocity.getDataHandle();
  
  DataHandle< CFreal > wallDistance = socket_wallDistance.getDataHandle();
  
  // resize socket
  wallShearStressVelocity.resize(wallDistance.size());

}

//////////////////////////////////////////////////////////////////////////////

void SA2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
  
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      deletePtr(m_cellGrads[iState][iVar]); 
    }
    
    m_cellGrads[iState].clear();
  }
  
  m_cellGrads.clear();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    SA2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
