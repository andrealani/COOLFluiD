#include "MathTools/MathConsts.hh"
#include "MathTools/JacobiEigenSolver.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/SUPG_ArtDiffStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<SUPG_ArtDiffStrategy,
                       FluctuationSplitData,
                       ArtificialDiffusionStrategy,
                       FluctSplitModule>
                       supgArtDiffStrategyProvider("SUPG");

//////////////////////////////////////////////////////////////////////////////

void SUPG_ArtDiffStrategy::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("WithShockDetection","Add shock detection to the SUPG artificial diffusion.");
   options.addConfigOption< CFreal >("Viscosity","Artifical viscosity.");
   options.addConfigOption< bool >("WithHughes","Add Hughes diffusion factor");
}

//////////////////////////////////////////////////////////////////////////////

SUPG_ArtDiffStrategy::SUPG_ArtDiffStrategy(const std::string& name) :
  ArtificialDiffusionStrategy(name),
  m_adv(0),
  m_tmpAD(0)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_with_shock_detect = true;
  setParameter("WithShockDetection",&m_with_shock_detect);
  
  m_userViscosity = 1;
  setParameter("Viscosity",&m_userViscosity);
  
  m_with_Hughes = false;
  setParameter("WithHughes", &m_with_Hughes);
  
}

//////////////////////////////////////////////////////////////////////////////

SUPG_ArtDiffStrategy::~SUPG_ArtDiffStrategy()
{
}

//////////////////////////////////////////////////////////////////////////////

void SUPG_ArtDiffStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ArtificialDiffusionStrategy::setup();

  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  CFuint dim   = PhysicalModelStack::getActive()->getDim();

  m_min_states.resize(nbEqs);
  m_max_states.resize(nbEqs);
  m_adv.resize(dim);
  m_result.resize(nbEqs);
  m_k.resize(nbEqs,nbEqs);
  m_v.resize(nbEqs);
  m_tmpAD.resize(nbEqs);
  theta.resize(nbEqs);
  tau.resize(nbEqs);
  fluxdot.resize(nbEqs);
  unitVector.resize(nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void SUPG_ArtDiffStrategy::addArtificialDiff(std::vector<RealVector>& residual)
{
  //compute the average jacobians on the cell
  _distribVar->computeJacobians();

  const CFuint nbEqs = m_min_states.size();
  const CFreal volume = m_cell->computeVolume();
  const CFreal h   = std::sqrt(4.0*MathTools::MathConsts::CFrealPi()*volume); // must get this from the m_cell

  // Use average cell jacobians to compute Hughes diffusion factor
  // Diffusion factor is set to 1 by default, i.e. no effect
  
  tau = 1;
  if (m_with_Hughes)
  {
    const std::vector<RealMatrix>& jacobs = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();    
    
   fluxdot = 0;
   unitVector = 1; 
    // Loop over the jacobian matrices and compute the sum of each of their rows squared and store the sums in
    // a vector (this is a sort of approximate multivector dot product)
    
    CFuint dim   = PhysicalModelStack::getActive()->getDim();
    
    for (CFuint i = 0; i<dim; i++){
      fluxdot +=  jacobs[i]*jacobs[i]*unitVector;
    }
    
    // Now compute the Hughes diffusion factor
    // The artificial diffusion residual is then later multiplied by this factor
    tau = 1/(2*sqrt(fluxdot));
  }
  
  // by default is one, meaning no shock detection
  theta = 1;
  if (m_with_shock_detect)
  {

    m_cellStates = m_cell->getStates();// store the pointers to state in another array (of RealVector*)

    m_min_states = MathTools::MathConsts::CFrealMax();
    m_max_states = -MathTools::MathConsts::CFrealMax();

    std::vector<State*>::iterator itr = m_cellStates->begin();

    for (;itr!=m_cellStates->end(); ++itr)
    {
      for (CFuint i = 0; i<nbEqs; ++i)
      {

        m_min_states[i] = std::min((**itr)[i],m_min_states[i]);
        m_max_states[i] = std::max((**itr)[i],m_max_states[i]);
      }

    }



    for (CFuint i = 0; i< nbEqs ; ++i)
    {
      theta[i] -= std::abs(m_max_states[i] - m_min_states[i])/(std::abs(m_max_states[i]) + std::abs(m_min_states[i])+10.0E-10);

    }
  }

  CFreal theta_real = theta.min();
  

  // only works for advection ?
//   _distribVar->getAverageAdvectionVector(m_adv,0);



//   static RealVector phi (PhysicalModelStack::getActive()->getNbEq());
//   phi = 0.0;
//   for (CFuint node = 0; node < residual.size(); ++node)
//   {
//     phi += residual[node];
//   }

//   cout << "Residual : " << phi << endl;

//   phi.abs(phi);

//   cout << "Abs residual : " << phi << endl;

#if 0
  const CFreal eps = 1E-10;
  const CFreal factor = ( phi[0] / ( m_adv.norm2() * volume ) ) + eps ;
  const CFreal theta = std::min ( 1.0 , 1.0 / factor );
#else

#endif

//   cout << "Theta : " << theta << endl;
//   cout << "Theta*h : " << theta*h << endl;


  // compute and add the artificial diffusion contribution to the nodal residual
  for (CFuint node = 0; node < residual.size(); ++node)
  {
    setIdx(node);
    getMethodData().getVolumeIntegrator()->integrateGeneralFunctorOnGeoEnt<SUPG_ArtDiffStrategy>(m_cell, *this, m_tmpAD);

    for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      residual[node][iEq] = m_userViscosity*tau[iEq]*theta_real*h*m_tmpAD[iEq];

  }
}

//////////////////////////////////////////////////////////////////////////////

 RealVector& SUPG_ArtDiffStrategy::operator()(const std::vector<Framework::State*>& vars,
                            const RealVector& shapeF,
                            const RealMatrix& grad,
                            Framework::GeometricEntity* const geo){
  // using average jacobians per cell
      // as an alternative we could calculate them per quadrature point
      const std::vector<RealMatrix>& jacob = *PhysicalModelStack::getActive()->getImplementor()->getJacobians();

//       CF_DEBUG_OBJ(grad);
//  CF_DEBUG_OBJ(shapeF);
      std::vector<State*>& states = *(geo->getStates());


//       CF_DEBUG_OBJ(states[0]->getCoordinates());
//       CF_DEBUG_OBJ(states[1]->getCoordinates());
//       CF_DEBUG_OBJ(states[2]->getCoordinates());
//       CF_DEBUG_OBJ(states[3]->getCoordinates());

//       CF_DEBUG_OBJ(*(states[0]));
//       CF_DEBUG_OBJ(*(states[1]));
//       CF_DEBUG_OBJ(*(states[2]));
//       CF_DEBUG_OBJ(*(states[3]));

//        CF_DEBUG_OBJ(jacob[XX]);
//        CF_DEBUG_OBJ(jacob[YY]);

CFuint m_dim   = PhysicalModelStack::getActive()->getDim();
      m_k = 0.0;
      m_v = 0.0;

      for (CFuint d = 0; d < m_dim; ++d)
      {
        m_k += grad(m_idx,d) * jacob[d];
        for (CFuint s = 0; s < states.size(); ++s)
        {
          m_v += jacob[d] * ( *(states[s]) * grad(s,d));
        }

      }
      m_result = m_k*m_v;
      return m_result;

#if 0
      CFreal A = jacob[XX](0,0);
      CFreal B = jacob[YY](0,0);
      CFreal u_x = 0.;
      CFreal u_y = 0.;
      for (CFuint s = 0; s < states.size(); ++s)
      {
        u_x += ( (*(states[s]))[0] * grad(s,XX));
        u_y += ( (*(states[s]))[0] * grad(s,YY));
      }
      cout << "Jacobians : " << A << " " << B << endl;
      cout << "Solution Grads : " << u_x << " " << u_y << endl;
      CFreal res = A * u_x + B * u_y;
      cout << "Recomputed residual : " << res << endl;

      m_v = res;
      m_result = res;
      CF_DEBUG_OBJ(m_result);
      return m_result;
#endif
}

    }// End namespace FluctSplit

}// End namespace COOLFluiD
