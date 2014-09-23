#ifndef COOLFluiD_Physics_ICP_ICPReactionTerm_hh
#define COOLFluiD_Physics_ICP_ICPReactionTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an ICP reaction physical term
 *
 * @author Radek Honzatko
 * @author Emanuele Sartori
 * @author Andrea Lani
 */
template <typename BASE>
class ICPReactionTerm : public BASE {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   *
   * The data are:
   * SIGMA - electric conductivity
   * MU0 - permeability of free space
   * F - torch operating frequency
   */
  enum {SIGMA=0, MU0=1, F=2};

  /**
   * Constructor without arguments
   */
  ICPReactionTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ICPReactionTerm();

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Resize the physical data
   */
  virtual void resizePhysicalData(RealVector& physicalData);

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const 
  {
    return BASE::getDataSize() + 3;
  }
  
  /**
   * Get the permeability of free space
   */
  CFreal getPermeability() const
  {
    return m_permeability;
  }
  
  /**
   * Get the torch operating frequency in MHz
   */
  CFreal getFreqMHz() const
  {
    return m_freqMHz;
  }

  /**
   * Get the electric conductivity
   */
  CFreal getElConductivity(const CFreal& Tdim, const CFreal& pdim) const;

  /**
   * Use the 2D model considering the fluxes in the z direction
   */
  bool use2DModel() const
  {
    return m_use2DModel;
  }

  /*
   * Using 1D model, it is still possible to simulate
   * torch and chamber (in chamber, 1D approx is no longer
   * valid), setting a max x (z) value for Ep source...
   */ 
  CFreal xMaxUsing1DModel() const
  {
    return m_xMaxUsing1DModel;
  }

  /*
   * Output electric field on file
   */ 
  std::string outputFileEM() const
  {
    return m_tecplotFileEM;
  }
  
  /*
   * Save rate for the file where the electric field is written
   */
  int outputFileEM_SaveRate() const
  {
    return m_tecplotFileEMSaveRate;
  }

  /**
   * extra verbose rate
   */
  CFint extraVerboseRate() const
  {
    return m_extraVerboseRate;
  }

  /**
   * Use approximation used in PEGASE (verification purposes)
   */
  bool pegaseApproximation()
  {
    return m_pegaseValues;
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the name
   */
  static std::string getName() 
  {
    return "ICPReactionTerm";
  }

protected:
  
  /// permeability of the free space
  CFreal m_permeability;

  /// EM field - Lorentz force output file
  std::string m_tecplotFileEM;
  int m_tecplotFileEMSaveRate;

  /// operating frequency of the torch [MHz]
  CFreal m_freqMHz;

  bool m_pegaseValues;

  /// Flag telling if to use the 2D model
  bool m_use2DModel;

  /// max x coord considering 1D model
  CFreal m_xMaxUsing1DModel;

  /// nb of iteration between extra messages
  CFint m_extraVerboseRate;

}; // end of class ICPReactionTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ICP/ICPReactionTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_ICP_ICPReactionTerm_hh
