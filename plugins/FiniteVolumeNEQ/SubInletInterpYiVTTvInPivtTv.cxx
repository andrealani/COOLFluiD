#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumeNEQ/SubInletInterpYiVTTvInPivtTv.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/OldLookupTable.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletInterpYiVTTvInPivtTv,
                      CellCenterFVMData,
		      FiniteVolumeModule>
subInletInterpYiVTTvInPivtTvFVMCCProvider("SubInletInterpYiVTTvInPivtTv");

//////////////////////////////////////////////////////////////////////////////

SubInletInterpYiVTTvInPivtTv::SubInletInterpYiVTTvInPivtTv(const std::string& name) :
  SubInletInterpYiVTTv(name),
  m_rhoi(),
  m_rhoiB()
{
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletInterpYiVTTvInPivtTv::~SubInletInterpYiVTTvInPivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletInterpYiVTTvInPivtTv::setGhostState(GeometricEntity *const face)
{   
  this->computeGhostPosition(face);
  
  // add lookup table in the nodal extrapolator
  SafePtr<StateInterpolator> interp = getStateInterpolator();
  getMethodData().getNodalStatesExtrapolator()->addLookupState
    (getCurrentTRS()->getName(), interp);
  
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  
  // this point [rho_i V T Tv] are set for the boundary state
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = m_library->getNbSpecies();
  const bool readFile = (m_yiVTTv.size() == 0);
  if (readFile) {
    const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] +
			       innerState.getCoordinates()[YY]);
    
    const CFuint nbEqs = innerState.size();
    for (CFuint i = 0; i < nbEqs; ++i) {
      // interpolated state value in input variables
      interp->interpolate(i, yCoord, (*m_tstate)[i]);
    }
    *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  }
  else {
    cf_assert(m_yiVTTv.size() == m_bstate->size());
    for (CFuint i = 0; i < m_yiVTTv.size(); ++i) {
      (*m_tstate)[i] = m_yiVTTv[i];
    }
    *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  } 

  (*m_bstate)[nbSpecies+1] = 0.;
  if (m_blowVelocity > 0.) {
    if ((*m_bstate)[nbSpecies] < m_blowVelocity) {
      (*m_bstate)[nbSpecies]   = m_blowVelocity;  
      (*m_bstate)[nbSpecies+1] = 0.;
    }
  }
  
  // p = R*T*(sum_i rho_i/M_i) + R*T_e*rho_e/M_e
  const CFreal RT = m_library->getRgas()*innerState[nbSpecies + dim];
  const CFreal RTb  = m_library->getRgas()*(*m_bstate)[nbSpecies + dim];
  
  // assumption here: second temperature is electron temperature
  const CFreal RTe = (m_library->presenceElectron()) ? m_library->getRgas()*innerState[nbSpecies + dim + 1] : RT;
  const CFreal RTeb = (m_library->presenceElectron()) ? m_library->getRgas()*(*m_bstate)[nbSpecies + dim + 1] : RTb;
  
  // rho_i = p_i/(R/M_i T_i)
  CFreal ovRhoInner = 0.;
  CFreal ovRhoIB = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal RT_i = (i != 0) ? RT : RTe; 
    m_rhoi[i] = innerState[i]/(RT_i/m_mmasses[i]);
    ovRhoInner += m_rhoi[i];
    
    const CFreal RTb_i = (i != 0) ? RTb : RTeb; 
    m_rhoiB[i] = (*m_bstate)[i]/(RTb_i/m_mmasses[i]);
    ovRhoIB += m_rhoiB[i];
  }
  ovRhoInner = 1./ovRhoInner;
  ovRhoIB = 1./ovRhoIB;
  
  const CFuint tempID = nbSpecies + dim;
  const CFuint nbEqs = innerState.size();
  
  CFuint idx = 0;
  for (CFuint i = 0; i < nbEqs; ++i) {
    cf_assert(idx < m_innerYiTTv.size());
    if (i < nbSpecies) {
      m_innerYiTTv[idx] = m_rhoi[i]*ovRhoInner;
      m_boundYiTTv[idx] = m_rhoiB[i]*ovRhoIB;
      idx++;
    }
    if (i >= tempID) {
      m_innerYiTTv[idx] = innerState[i];
      m_boundYiTTv[idx] = (*m_bstate)[i];
      idx++;
    }
  }
  
  this->repositionNode(m_innerYiTTv, m_ghostYiTTv, m_boundYiTTv, m_minYiTTv); 
    
  // reset the ghost node with the new position
  ghostState.getCoordinates() = this->m_tempGhostNode;
  
  for (CFuint i = 0; i < nbEqs; ++i) {
    // velocity components
    if (i >= nbSpecies && i < tempID) {
      this->linearInterpolate(innerState[i], (*m_bstate)[i], ghostState[i]);
    }
    else if (i < nbSpecies) {
      // p_i_ghost = p_i_in
      cf_assert(m_ghostYiTTv[i] >= 0.);
      ghostState[i] = innerState[i];
    }
    else if (i >= tempID) {
      ghostState[i] = m_ghostYiTTv[i - dim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void SubInletInterpYiVTTvInPivtTv::setup()
{
  SubInletInterpYiVTTv::setup();
  
  m_rhoi.resize(m_library->getNbSpecies()); 
  m_rhoiB.resize(m_library->getNbSpecies());
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
