#ifndef COOLFluiD_Numerics_FluctSplit_NavierStokesTermHOIsoP2_hh
#define COOLFluiD_Numerics_FluctSplit_NavierStokesTermHOIsoP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DiffusiveVarSet.hh"
#include "FluctSplit/ComputeDiffusiveTerm.hh"
#include "FluctSplit/P2Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
      class EulerVarSet;
    }
  }



    namespace FluctSplit {
            
//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the diffusive flux corresponding to the Navier Stokes
 * physical model for isoparametric P2 elements
 *
 * @author Martin Vymazal
 * @author Andrea Lani
 *
 */
class NavierStokesTermHOIsoP2 : public ComputeDiffusiveTerm {
public:

  /**
   * Constructor
   */
  NavierStokesTermHOIsoP2(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokesTermHOIsoP2();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  void setup();
  
  /**
   * Set the update variable set
   */
  void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar);
   
  /**
   * Compute the diffusive term flux in the current cell
   */
  void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
			    std::vector<RealVector>& result,
			    const bool updateCoeff);
  /**
    * Set the update variable set
    */
  void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>); 
// Integral of diffusive part times basis function

  void fluctuation_diff_galerkin(const CFuint& i0, const CFuint& i1, const CFuint& i2);

// Integral of diffusive part times bubblefunction
  /**
    * i0,i1,i2 are vertices of the subtriangle
    * f0, f1, f2 are nodal faces, i.e. f0 is face of subtriangle i0-i1-i2 which is opposite to node i0
    */
  void fluctuation_diff_bubble(const CFuint& i0, const CFuint& i1, const CFuint& i2);

  void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob){
throw Common::NotImplementedException (FromHere(),"NavierStokesTermHOisoP2::computePicardJacob()");
  
}
protected: // data

  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;

  // acquaintance of the update convective var set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVar;
    // array of cell states
  std::vector<Framework::State*> * m_cellStates;

  // array of cell nodes
  std::vector<Framework::Node*> * m_cellNodes;

  // array of cell states
  std::vector<RealVector*> m_states;

  // array of values (rho, u, v, T)
  RealMatrix m_values;

  // array of gradients
  std::vector<RealVector*> m_gradients;

  // array of average values (rho, u, v, T)
  RealVector m_avValues;

  // unit normal
  RealVector m_normal;

  /// fluctuation from diffusion part bubble part
  RealVector* m_phi_diff_bub;

 /// fluctuation from diffusion part galerkin
  std::vector<RealVector*> m_phi_diff_gal;

/// fluctuation from diffusion part bubble part splitted with the beta of LDA normally
   std::vector<RealVector*> m_phi_diff_bub_split;

/// Pointer to the current cell being processed.
  CFuint m_cellID;

  /// quadrature points per cell
  RealVector m_qd0;
  RealVector m_qd1;
  RealVector m_qd2;
  RealVector m_wqd;
/// states at quadrature points
  RealVector m_qdstates;

  /// quadrature points per face
//   RealVector m_faceQd0;
//   RealVector m_faceQd1;
//   RealVector m_faceWQd;

  ///Coordinates of vertices of current triangle
//   RealVector m_x;
//   RealVector m_y;

  ///Derivatives of shape functions in reference space
//   RealVector m_dNdxi;
//   RealVector m_dNdeta;
  ///Derivatives of shape functions in physical space, Jacobian of transformation physical->reference space
  RealVector m_dNdx;
  RealVector m_dNdy;
  CFreal m_J;

  CFreal m_dxdxi, m_dxdeta, m_dydxi, m_dydeta;

  RealVector m_faceNormal0;
  RealVector m_faceNormal1;
  RealVector m_SFGrad;

  ///Gradients of states:
  RealVector m_gradState_x;
  RealVector m_gradState_y;

  ///Vectors of viscous fluxes:
  RealVector m_F0;
  RealVector m_F1;
  RealVector m_F2;
  RealVector m_F3;
  RealVector m_F4;
  RealVector m_F5;
  RealVector m_F;

  ///Coordinates of nodes in reference space
  static CFreal xi_ref[6];
  static CFreal eta_ref[6];
  ///Names of vertices in reference space
  enum { V0, V1, V2, V3, V4, V5 };

  ///Number of nodes in P2 triangle:
  enum { NNODES = 6 };

  ///Object to compute normals of P2P2 triangle on the fly
  P2Normal m_CP2N;

  ///Number of equations
  CFuint m_nbEqs;

/// Coeficients used to copute the diffusion residual
  MathTools::CFMat<CFreal> kappa;

// volume of the cell
  CFreal m_cellVolume;

  //Helper variables:
  RealVector m_refCoord;
  RealVector m_physCoord;
  CFreal m_L0, m_L1, m_L2;
  RealVector m_gradS;

  //Functions
  //Compute mapped coordinate from given reference coordinate
  void ComputeMappedCoord(const std::vector<Framework::Node*> & nodes, const RealVector & refcoord, RealVector & physcoord);

  //Compute gradients of shape functions at given integration point
  void ComputeSFGradients(const std::vector<Framework::Node*>& nodes, const CFreal xi, const CFreal eta, CFreal J, RealVector & dNdx, RealVector & dNdy);

//Compute gradient of bubble function at given integration point
  void ComputeBubbleGradient(const std::vector<Framework::Node*>& nodes, const CFreal xi, const CFreal eta, CFreal J, RealVector & gradS);


}; // end of class NavierStokesTermHOIsoP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NavierStokesTermHOIsoP2_hh
