#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumeNEQ/SubInletInterpYiVTTv.hh"
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

MethodCommandProvider<SubInletInterpYiVTTv,
                      CellCenterFVMData,
		      FiniteVolumeModule>
subInletInterpYiVTTvFVMCCProvider("SubInletInterpYiVTTv");

//////////////////////////////////////////////////////////////////////////////

void SubInletInterpYiVTTv::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >
    ("YiVTTv", "Inlet mass fractions, velocity components and temperatures.");
  
  options.addConfigOption< CFreal >
    ("BlowVelocity", "Blowing velocity to stabilize the flow.");
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletInterpYiVTTv::SubInletInterpYiVTTv(const std::string& name) :
  SuperInletInterp(name),
  m_library(CFNULL),
  m_mmasses(),
  m_ghostYiTTv(),
  m_innerYiTTv(),
  m_boundYiTTv(),
  m_minYiTTv()
{
  addConfigOptionsTo(this);
  
  m_yiVTTv = std::vector<CFreal>();
  setParameter("YiVTTv",&m_yiVTTv);
  
  m_blowVelocity = 0.;
  setParameter("BlowVelocity",&m_blowVelocity);
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletInterpYiVTTv::~SubInletInterpYiVTTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletInterpYiVTTv::setGhostState(GeometricEntity *const face)
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
    const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] + innerState.getCoordinates()[YY]);
    const CFuint nbEqs = innerState.size();
    for (CFuint i = 0; i < nbEqs; ++i) {
      // interpolated state value in input variables
      interp->interpolate(i, yCoord, (*m_tstate)[i]);
    }
    *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  }
  else {
    assert(m_yiVTTv.size() == m_bstate->size());
    for (CFuint i = 0; i < m_yiVTTv.size(); ++i) {
      (*m_tstate)[i] = m_yiVTTv[i];
    }
    *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
  } 
  
  if (readFile) {
    (*m_bstate)[nbSpecies+1] = 0.;
  }
  if (m_blowVelocity > 0.) {
    if ((*m_bstate)[nbSpecies] < m_blowVelocity) {
      (*m_bstate)[nbSpecies]   = m_blowVelocity;  
      (*m_bstate)[nbSpecies+1] = 0.;
    }
  }
  
  // compute internal pressure
  // p = R*T*(sum_i rho_i/M_i) + R*T_e*rho_e/M_e
  CFreal pInner = 0.;
  CFreal rhoInner = 0.;
  CFreal rhoIB = 0.;
  const CFreal RT = m_library->getRgas()*innerState[nbSpecies + dim];
  // assumption here: second temperature is electron temperature
  const CFreal RTe = (m_library->presenceElectron()) ? m_library->getRgas()*innerState[nbSpecies + dim + 1] : RT;
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    rhoInner += innerState[i];
    rhoIB += (*m_bstate)[i];
    
    (i != 0) ?
      pInner += (innerState[i]/m_mmasses[i])*RT :
      pInner += (innerState[0]/m_mmasses[0])*RTe;
  }
  
  const CFuint tempID = nbSpecies + dim;
  const CFuint nbEqs = innerState.size();
  
  CFuint idx = 0;
  for (CFuint i = 0; i < nbEqs; ++i) {
    cf_assert(idx < m_innerYiTTv.size());
    if (i < nbSpecies) {
      m_innerYiTTv[idx] = innerState[i]/rhoInner;
      m_boundYiTTv[idx] = (*m_bstate)[i]/rhoIB;
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
  
  const CFreal RTghost  = m_library->getRgas()*m_ghostYiTTv[nbSpecies];
  // assumption here: second temperature is electron temperature
  const CFreal RTeGhost = (m_library->presenceElectron()) ? m_library->getRgas()*m_ghostYiTTv[nbSpecies+1] : RTghost;
  
  // rho_G = p_G/(R T_G_i sum_i y_i_G/M_i)
  // pInner = pGhost
  CFreal RTovMGhost = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    (i != 0) ? 
      RTovMGhost += m_ghostYiTTv[i]/m_mmasses[i]*RTghost : 
      RTovMGhost += m_ghostYiTTv[0]/m_mmasses[0]*RTeGhost;
  }
  const CFreal rhoG = pInner/RTovMGhost;
  
  for (CFuint i = 0; i < nbEqs; ++i) {
    // velocity components
    if (i >= nbSpecies && i < tempID) {
      this->linearInterpolate(innerState[i], (*m_bstate)[i], ghostState[i]);
    }
    else if (i < nbSpecies) {
      // rho_i_ghost = rho_i_in * T_in / T_ghost 
      cf_assert(m_ghostYiTTv[i] >= 0.);
      ghostState[i] = m_ghostYiTTv[i]*rhoG;
    }
    else if (i >= tempID) {
      ghostState[i] = m_ghostYiTTv[i - dim];
    }
  }
}
      
// //////////////////////////////////////////////////////////////////////////////

// void SubInletInterpYiVTTv::setGhostState(GeometricEntity *const face)
// {
//   // add lookup table in the nodal extrapolator
//   getMethodData().getNodalStatesExtrapolator()->addLookupState(getCurrentTRS()->getName(), &m_lookupState);
  
//   State& innerState = *face->getState(0);
//   State& ghostState = *face->getState(1);
  
//   // this point [rho_i V T Tv] are set for the boundary state
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   const CFuint nbSpecies = m_library->getNbSpecies();

//   const bool readFile = (m_yiVTTv.size() == 0);
  
//   if (readFile) {
//     const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] +
// 			       innerState.getCoordinates()[YY]);
    
//     const CFuint nbEqs = innerState.size();
//     for (CFuint i = 0; i < nbEqs; ++i) {
//       // interpolated state value in input variables
//       (*m_tstate)[i] = m_lookupState[i]->get(yCoord);
//     }
//     *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
//   }
//   else {
//     assert(m_yiVTTv.size() == m_bstate->size());
//     for (CFuint i = 0; i < ghostState.size(); ++i) {
//       (*m_bstate)[i] = m_yiVTTv[i]; // y_i is set for the first nbSpecies components
//     }
//   }
  
//   // if (m_blowVelocity > 0.) {
//   //(*m_bstate)[nbSpecies] = std::max((*m_bstate)[nbSpecies], m_blowVelocity);  
//   (*m_bstate)[nbSpecies] = std::max((*m_bstate)[nbSpecies], 50.);  
//   // }
  
//   // set ghost values for V,T,Tv
//   for (CFuint i = nbSpecies; i < ghostState.size(); ++i) {
//     ghostState[i] = 2.*(*m_bstate)[i] - innerState[i];
    
//     // fix for negative temperatures
//     if (i >= nbSpecies+dim && ghostState[i] < 0.) {
//       ghostState[i] = (*m_bstate)[i];
//       //cout << "fix T" << endl;
//    }
//   }
  
//   // compute internal pressure
//   // p = R*T*(sum_i rho_i/M_i)
//   CFreal pInner = 0.;
//   CFreal rhoInner = 0.;
//   CFreal rhoIB = (readFile) ? 0. : 1.;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     rhoInner += innerState[i];
//     if (readFile) {rhoIB += (*m_bstate)[i];}
//     pInner   += innerState[i]/m_mmasses[i]; // Y_i
//   }
//   pInner *= m_library->getRgas()*innerState[nbSpecies + dim];
  
//   bool negativeYi = false;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     // initially yiGhost is stored into the ghost state instead of rho_i
//     // y_G = 2*y_B - y_I
//     ghostState[i] = 2.*(*m_bstate)[i]/rhoIB - innerState[i]/rhoInner;
//     // assert((*m_bstate)[i]/rhoIB < 1.000001);
//     // assert(innerState[i]/rhoInner < 1.000001);	 

//     if (ghostState[i] < 0.) {
//       // cout << "fix rho_in " << 0.5*(ghostState.getCoordinates()[YY] + innerState.getCoordinates()[YY])  << endl;
//       //  cout << i << " => 2*" << (*m_bstate)[i]/rhoIB << " - " <<  innerState[i]/rhoInner << endl;
//       negativeYi = true;
//     //  break;
//     }
//   }
  
//   // fix for negative mass fractions
//   if (negativeYi) {
//     //cout << "SubInletInterpYiVTTv => Correcting mass fractions < 0" << endl;
//     for (CFuint i = 0; i < nbSpecies; ++i) {
//       ghostState[i] = (*m_bstate)[i];
//     }
//   }
  
//   CFreal sumYiOvMGhost = 0.;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     sumYiOvMGhost += ghostState[i]/m_mmasses[i];
//   }
  
//   // rho_G = p_G/(sum_i y_i_G/M_i)/(R T_G)
//   // pInner = pGhost
//   const CFreal rhoG = pInner/(m_library->getRgas()*sumYiOvMGhost*ghostState[nbSpecies + dim]);
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     ghostState[i] *= rhoG;
//   }
// }

// // //////////////////////////////////////////////////////////////////////////////

// void SubInletInterpYiVTTv::setGhostState(GeometricEntity *const face)
// {
//    this->computeGhostPosition(face);
//  State& innerState = *face->getState(0);
//   State& ghostState = *face->getState(1);
  
//   // this point [rho_i V T Tv] are set for the boundary state
//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEqs();
//   const CFuint nbSpecies = m_library->getNbSpecies();
//   const bool readFile = (m_yiVTTv.size() == 0);
  
//   if (readFile) {
//     const CFreal yCoord = 0.5*(ghostState.getCoordinates()[YY] +
// 			       innerState.getCoordinates()[YY]);
    
//     for (CFuint i = 0; i < nbEqs; ++i) {
//       // interpolated state value in input variables
//       (*m_tstate)[i] = m_lookupState[i]->get(yCoord);
//     }
//     *m_bstate = *m_inputToUpdateVar->transform(m_tstate);
//   }
//   else {
//     assert(m_yiVTTv.size() == m_bstate->size());
//     for (CFuint i = 0; i < nbEqs; ++i) {
//       (*m_bstate)[i] = m_yiVTTv[i]; // y_i is set for the first nbSpecies components
//     }
//   }
  
//    if (m_blowVelocity > 0.) {
//     (*m_bstate)[nbSpecies] = std::max((*m_bstate)[nbSpecies], m_blowVelocity);  
//  }
      
//   // set the ghost values
//   // if needed, reposition ghost state so that the corresponding ghost temperature is > 0
//   const CFuint tempID = nbSpecies + dim;
//   CFuint iTemp = tempID;
//   for (CFuint i = 0; i < m_innerTTvib.size(); ++i, ++iTemp) {
//     m_innerTTvib[i] = innerState[iTemp];
//     m_ghostTTvib[i] = ghostState[iTemp];
//   }
//   this->repositionNode(m_innerTTvib, m_ghostTTvib, (*m_bstate)[tempID], this->m_ghostTempMin);
  
//   // reset the ghost node with the new position
//   ghostState.getCoordinates() = this->m_tempGhostNode;
  
//   // compute internal pressure
//   // p = R*T*(sum_i rho_i/M_i)
//   CFreal pInner = 0.;
//   CFreal rhoInner = 0.;
//   CFreal rhoIB = (readFile) ? 0. : 1.;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     rhoInner += innerState[i];
//     if (readFile) {rhoIB += (*m_bstate)[i];}
//     pInner   += innerState[i]/m_mmasses[i]; // Y_i
//   }
//   pInner *= m_library->getRgas()*innerState[tempID];
  
//   // extrapolate the mass fractions in the ghost in order to respect the imposed boundary values
//   bool negativeYi = false;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     // initially yiGhost is stored into the ghost state instead of rho_i
//     // y_G = linear_interp(y_B, y_I)
//     this->linearInterpolate(innerState[i]/rhoInner, (*m_bstate)[i]/rhoIB, ghostState[i]);
//     if (ghostState[i] < 0.) {
//       negativeYi = true;
//       break;
//     }
//   }
  
//   // if the ghost mass fractions are < 0 reset them all to be equal to the boundary values
//   if (negativeYi) {
//     cout << "SubInletInterpYiVTTv => Correcting mass fractions < 0" << endl;
//     for (CFuint i = 0; i < nbSpecies; ++i) {
//       ghostState[i] = (*m_bstate)[i];
//     }
//   }
  
//   CFreal sumYiOvMGhost = 0.;
//   for (CFuint i = 0; i < nbSpecies; ++i) {
//     sumYiOvMGhost += ghostState[i]/m_mmasses[i];
//   }
  
//   // pInner = pGhost
//   // rho_G = p_G/(sum_i y_i_G/M_i)/(R T_G)
//   const CFreal rhoG = pInner/(m_library->getRgas()*sumYiOvMGhost*m_ghostTTvib[0]);
  
//   for (CFuint i = 0; i < nbEqs; ++i) {
//     // velocity components
//     if (i >= nbSpecies && i < tempID) {
//       this->linearInterpolate(innerState[i], (*m_bstate)[i], ghostState[i]);
//     }
//     else if (i < nbSpecies) {
//       // rho_i_ghost = rho_i_in * T_in / T_ghost
//       ghostState[i] *= rhoG;
//     }
//     else if (i >= tempID) {
//       ghostState[i] = m_ghostTTvib[i - tempID];
//     }
//   }
// }
      
//////////////////////////////////////////////////////////////////////////////

void SubInletInterpYiVTTv::setup()
{
  SuperInletInterp::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = m_library->getNbSpecies();
  m_mmasses.resize(nbSpecies);
  m_library->getMolarMasses(m_mmasses);

  const bool readFile = (m_yiVTTv.size() == 0);
  const CFuint yiTTvSize = (readFile) ?  nbSpecies + m_library->getNbTempVib() + 1 :
    m_yiVTTv.size(); // size can be bigger than case in which file is read 
                     // AL: gory fix, to be improved by allowin g to read more variables  
  
  m_ghostYiTTv.resize(yiTTvSize); // y_i, roto-translational + vibrational temperatures
  m_innerYiTTv.resize(yiTTvSize); // y_i, roto-translational + vibrational temperatures
  m_boundYiTTv.resize(yiTTvSize); // y_i, roto-translational + vibrational temperatures
  
  // set minimum ghost values for Yi, T , Tv
  m_minYiTTv.resize(yiTTvSize, 0.);
  for (CFuint i = nbSpecies; i < yiTTvSize; ++i) {
    m_minYiTTv[i] = 40.; // minimum temperature allowed in the ghost state
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
