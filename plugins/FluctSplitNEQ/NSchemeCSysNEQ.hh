#ifndef COOLFluiD_Numerics_FluctSplitNEQ_NSchemeCSysNEQ_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_NSchemeCSysNEQ_hh

//////////////////////////////////////////////////////////////////////////////

// This class derives from:
#include "FluctSplit/NSchemeCSys.hh"

#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the N scheme for RDS space discretization
 * based on the CRD approach.
 * It provides extra dissipation to face the problem of non-positivity
 * across strong bow shock waves, which plagues the standard Nc scheme. 
 * @author Andrea Lani
 * @author Jesús Garicano Mena
 *
 */

class NSchemeCSysNEQ : public FluctSplit::NSchemeCSys {
public:

  /**
   * Default constructor.
   */
  explicit NSchemeCSysNEQ(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NSchemeCSysNEQ();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Configure
   */

  void configure(Config::ConfigArgs& args);

  /**
   * Define the configuration options
   */

  void static defineConfigOptions(Config::OptionList& options);

  /**
   * Distribute the residual
   */
  virtual void distribute(std::vector<RealVector>& residual);
  /**
   * Distribute the residual
   */
  virtual void distributePart(std::vector<RealVector>& residual);
  
  /**
   * Compute all the contributions for the Picard jacobian
   */
  void computePicardJacob(std::vector<RealMatrix*>& jacob);


  /**
   * Store cells where extra dissipation is active
   */
  void store_ExtraDiss(const CFreal alpha_max);


  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = NSchemeCSys::needsSockets(); 

    result.push_back(&socket_ExtraDiss_active);
    return result;
  }
  
protected:
  /// Current global iteration
  static CFuint _actual_iter;

  /// last iteration in which store_CarbuncleFixActive() and store_ArtViscCoeff() were accessed
  static CFuint _last_accessed_at_iter;

  /// index of the  variable chosen for blending coefficient
  CFuint _varID;

//   /// flag telling if the locations where the extra dissipation is applied have to be stored
//   CFuint _storeExtraD;
  /// flag telling NOT to recompute the additional dissipation anymore
  CFuint _freezeAdditionalDiss;

  /// Internal variable telling if extra dissipation should be added
  CFreal _doAct;

  /// pmax - pmin
  CFreal _deltaP;

  /// reference length
  CFreal _length;

  /// reference speed
  CFreal _speed;

  /// Reference difference on var
  CFreal _deltaVar;

  /// name of the variable chosen for blending coefficient
  std::string _varName;

  /// Vector holding physical data
  RealVector _pData;

  /// Vector holding an average of the Conservative variables
  RealVector _Uavg;
  
  /// gradiant variables
  RealVector _gradVar;
  
  /// pointer holding the physical model
  Common::SafePtr<COOLFluiD::Framework::MultiScalarTerm<COOLFluiD::Physics::NavierStokes::EulerTerm> > _model;  

  /// socket for storage of the cells where the extra dissipation is applied
  Framework::DataSocketSink<CFreal> socket_ExtraDiss_active;

}; // end of class NSchemeCSysNEQ

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NSchemeCSys_hh
