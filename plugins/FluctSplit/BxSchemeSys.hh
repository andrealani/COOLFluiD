#ifndef COOLFluiD_Numerics_FluctSplit_BxSchemeSys_hh
#define COOLFluiD_Numerics_FluctSplit_BxSchemeSys_hh

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents the N scheme for RDS space discretization
   *
   * @author Andrea Lani
   */
template <class BASE, class MODEL>
class BxSchemeSys : public BASE {
public:
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit BxSchemeSys(const std::string& name);

  /// Destructor
  virtual ~BxSchemeSys();

  /// Configure this object with user defined parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Setup this object with data depending on the mesh
  virtual void setup();

protected:
  
  /**
   * Compute the blending coefficients
   */
  virtual void computeBlendingCoeff();
  
  /**
   * Update the scaling values
   */
  virtual void updateScalingValues();


  /**
   * Add an isotropic dissipative term.
   */
   virtual void addExtraDissipation(std::vector<RealVector>& residual);
    
protected:

  /// convective term
  Common::SafePtr<MODEL> _cterm;
  
  /// temporary physical data array
  RealVector _pData;
  
  /// temporary gradient
  RealVector _grad;
  
  /// ID of the variable chosen for blending coefficient
  CFuint _varID;
  
  /// var_Max - var_min
  CFreal _deltaVar;
  
  /// reference length
  CFreal _length;
  
  /// reference speed
  CFreal _speed;

  /// name of the variable chosen for blending coefficient
  std::string _varName;

/// Order of the discretization
  CFuint _order;

  /// Shock detector that we want to use
  std::string _sh_detector;


protected:


  /// Vector holding an average of the Conservative variables
  RealVector _Uavg;

  /// pointer holding the physical model
  Common::SafePtr<COOLFluiD::Framework::MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> > _model;

  /// socket for storage of the cells where the extra dissipation is applied
  //   Framework::DataSocketSink<CFreal> socket_ExtraDiss_active;

}; // end of class BxSchemeSys

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace FluctSplit



} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "BxSchemeSys.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_BxSchemeSys_hh
