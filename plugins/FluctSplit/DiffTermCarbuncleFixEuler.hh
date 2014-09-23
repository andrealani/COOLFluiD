#ifndef COOLFluiD_Numerics_FluctSplit_DiffTermCarbuncleFixEuler_hh
#define COOLFluiD_Numerics_FluctSplit_DiffTermCarbuncleFixEuler_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeDiffusiveTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes an artificial diffusive term to cure 
 * the carbuncle phenomenon for the Euler equations 
 *
 * @author Andrea Lani
 * @author Jesus Garicano Mena
 *
 */
template <typename UPDATEVAR>
class DiffTermCarbuncleFixEuler : public ComputeDiffusiveTerm {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  DiffTermCarbuncleFixEuler(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~DiffTermCarbuncleFixEuler();
  
  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Set the update variable set
   */
  virtual void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffVar)
  {
  }
 
  /**
   * Compute the diffusive term flux in the current cell
   */
  virtual void computeDiffusiveTerm(Framework::GeometricEntity *const geo,
				    std::vector<RealVector>& result,
				    const bool updateCoeff);
  
  /**
   * Set the update variable set
   */
  virtual void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet>); 
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeDiffusiveTerm::needsSockets();
    
    result.push_back(&socket_ArtViscCoeff);
	result.push_back(&socket_ViscCoeff);
    result.push_back(&socket_fix_active);
    result.push_back(&socket_updateCoeff);

    result.push_back(&socket_uCsi);
    result.push_back(&socket_uEta);

    result.push_back(&socket_duCsidCsi);
    result.push_back(&socket_duEtadCsi);
    result.push_back(&socket_duCsidEta);
    result.push_back(&socket_duEtadEta);
	
	result.push_back(&socket_dpdCsi);
    result.push_back(&socket_dpdEta);
    
    return result;
  }

private:
  
  /**
   * Compute Picard jacobian
   */
  void computePicardDiffJacob(Framework::GeometricEntity *const cell,std::vector<RealMatrix*>& _diffjacob)
  {
    throw Common::NotImplementedException (FromHere(),"DiffTermCarbuncleFixEuler::computePicardJacob()");
  }

protected:

  /**
   * Store the cells where the carbuncle fix is active
   */
  void store_FixActiveCells();

  /**
   * Store the cells where the
   */
  void store_ExtraInfo(const CFreal uCsi, const CFreal uEta, const CFreal duCsidCsi, const CFreal duEtadCsi, const CFreal duCsidEta, const CFreal duEtadEta, const CFreal dpdCsi, const CFreal dpdEta);

  /**
   * Store the cell Peclet number
   */
  void store_PeCell(CFreal thePecletNb);
  
  /// Tells if this term can be added to another derived one
  virtual bool addToDerivedTerm() {return true;}
  
protected: // data
  
  RealVector _normal0;
  
  RealVector _normal1;
  
  RealVector _normal2;
   
  RealVector _gradientM2;
  
  RealVector _nodalM2;
  
  RealVector _gradientT;
  
  RealVector _nodalT;
  
  RealVector _gradientP;
  
  RealVector _nodalP;
  
  RealVector _unitGradientM2;
  
  RealVector _uCsiAtNodes;
  
  RealVector _uEtaAtNodes;
  
  RealVector _nEtaAtNodes;
  
  // dimension
  CFreal _dim;
  
  // acquaintance of the update convective var set
  Common::SafePtr<UPDATEVAR> _updateVar;
  
  /// socket for storage of artificial diffusion coefficient mu_s
  Framework::DataSocketSink<CFreal> socket_ArtViscCoeff;
  
  /// socket for storage of the physical diffusion coefficient mu
  Framework::DataSocketSink<CFreal> socket_ViscCoeff;
  
  /// socket for storage of the cells where the fix is active
  Framework::DataSocketSink<CFreal> socket_fix_active;

  /// socket for storage of uCsi, uEta and their derivatives
  Framework::DataSocketSink<CFreal> socket_uCsi;
  Framework::DataSocketSink<CFreal> socket_uEta;
  Framework::DataSocketSink<CFreal> socket_duCsidCsi;
  Framework::DataSocketSink<CFreal> socket_duEtadCsi;
  Framework::DataSocketSink<CFreal> socket_duCsidEta;
  Framework::DataSocketSink<CFreal> socket_duEtadEta;
  
  Framework::DataSocketSink<CFreal> socket_dpdCsi;
  Framework::DataSocketSink<CFreal> socket_dpdEta;
  
  /// array of values for the speed components
  std::vector<RealVector> _speed;
  
  ///  velocity gradient
  RealVector _pdata;
   
  ///  cell average speed
  RealVector _cell_avg_speed;
   
  /// flag telling if a supersonic/subsonic transition happens in the cell
  bool _compression_through_sonicCell;
  
  /// mesh size
  CFreal _hsize;
  
  /// For method store_FixActiveCells() to know where to insert it.
  CFuint _cellID;
  
  /// dissipation coefficient
  CFreal _dissipEps;
    
  /// flag telling if the fix has to activated
  CFuint _activateFix;

  /// flag ....
  CFuint _variantChosen;
  
  /// Isotropic treatment of the dissipation?
//  CFuint _isotropic;
  
  /// Include energy dissipative terms?
  CFuint _includeEnergyDissipation;
  
  /// flag telling not to recompute the artificial viscous coefficient mu_s
  CFuint _freeze_mu_s;

  /// flag telling if the temperature gradients are included in the energy consistent fix
  CFuint _includeTemperatureGradient;
   
  /// Inform of the cells where the fix is active
  bool _store_fix_active_cells;
      
  /// Current global iteration
  static CFuint _present_iter;

  /// last iteration in which store_CarbuncleFixActive() and store_ArtViscCoeff() were accessed
  static CFuint _last_accessed_at_iter;

  /// last iteration in which store_CarbuncleFixActive() and store_ArtViscCoeff() were accessed
  static CFuint _last_accessed_at_iter__OUTER;
  
}; // end of class DiffTermCarbuncleFixEuler
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "DiffTermCarbuncleFixEuler.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_DiffTermCarbuncleFixEuler_hh
