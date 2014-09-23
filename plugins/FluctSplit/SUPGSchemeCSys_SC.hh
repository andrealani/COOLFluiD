#ifndef COOLFluiD_Numerics_FluctSplit_SUPGSchemeCSys_SC_hh
#define COOLFluiD_Numerics_FluctSplit_SUPGSchemeCSys_SC_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/NSchemeCSys.hh"
#include "FluctSplit/BSchemeBase.hh"

#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the SUPG scheme for RDS space discretization
 * based on the CRD approach
 *
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 *
 */
class SUPGSchemeCSys_SC : public BSchemeBase<NSchemeCSys> {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor
   */
   explicit SUPGSchemeCSys_SC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SUPGSchemeCSys_SC();

  /// Configure this object with user defined parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

  /// Distribute the residual
  virtual void distribute(std::vector<RealVector>& residual);

  ///  Distribute part of the residual
  virtual void distributePart(std::vector<RealVector>& residual);
 
protected: // functions

  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);

  /// Compute the blending coefficients
  virtual void computeBlendingCoeff();

  /**
   * Add an isotropic dissipative term.
   */
  virtual void addExtraDissipation(std::vector<RealVector>& residual);
  
  /**
   * ...
   */
   // not needed virtual CFreal getDissFactor();
  
  /**
   * ...
   */  
  // not needed virtual CFreal getMa();
  
  /**
   * ...
   */  
  // not needed virtual RealVector getCellAvgSpeed();
  
  /**
   * ...
   */    
  // not needed virtual RealVector getPressureBounds();
  
  // -- SUPG KIRK
  /**
   * ...
   */
  CFreal getNormSq(const RealVector vectorW);
  

protected: // data

  RealVector  _sumKplusU;

  RealMatrix  _sumKplus;

  RealMatrix  _invK;

  RealVector  _uInflow;

  RealVector  _uDiff;

  RealVector  _temp;
  
  RealVector  _tempBkp;
  
  RealMatrix _tempMat;

  RealMatrix _tmp;
  
  RealVector _sumKU;


  /// Auxiliary vectors - stabilization terms
  RealVector  _stab_Adv;
  RealVector  _stab_Shock;

  /// Matrix of dissipative time-scales
  RealMatrix _invTau;

  /// Matrix of dissipative time-scales
  RealMatrix _Tau;

  /// Jacobian matrices of the Advective Term
  RealMatrix _Ax;
  RealMatrix _Ay;

  /// Average state over the cell
  RealVector _Uavg;

  /// Unit vector indicating speed direction
  RealVector _speedNormal;

  /// Vector holding physical data
  RealVector _pData;
  
  // -- Shock capturing term --  
  /// Vectors holding gradXi and gradEta with respect to global x and y coordinates
  RealVector _gradXi;
  RealVector _gradEta;
  
  /// Vectors holding dUdx, dUdy and dUdz (i.e. components of the gradient of conserved variables)
  RealVector _dUdx;
  RealVector _dUdy;
  RealVector _dUdz;

  /// Transformation matrix from Conserved to Symmetrizing variables  
  RealMatrix _dVdU;
  RealMatrix _dUdV;

  /// Transformation matrix needed to compute gradXi and gradEta
  RealMatrix _dUHdU;

  /// Transformation matrix needed to compute gradXi and gradEta
  RealMatrix _M;

//   /// gIJ matrix
//   RealMatrix _gIJ;

  // --------------------------------------------------------------------------
  
  /// flag telling if to run first order
  CFuint     m_firstOrder;

  /// LDA residual ???
  RealVector m_phiLDA;

  // --------------------------------------------------------------------------
  //void distributePart( std::vector<RealVector>& residual );

}; // end of class SUPGSchemeCSys_SC

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SUPGSchemeCSys_SC_hh
