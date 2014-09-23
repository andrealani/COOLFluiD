#ifndef COOLFluiD_Numerics_FluctSplit_SUPG_ArtDiffStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_SUPG_ArtDiffStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include <boost/concept_check.hpp>


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Artificial Diffusion strategy
 *
 * @author Nadege Villedieu
 * @author Gabriel Maher (Hughes diffusion factor)
 */
class SUPG_ArtDiffStrategy : public ArtificialDiffusionStrategy {

public: // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  SUPG_ArtDiffStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~SUPG_ArtDiffStrategy();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    ArtificialDiffusionStrategy::configure(args);
  }
void addArtificialDiff(std::vector<RealVector>& residual) ;

  /**
   * Set up private data and data
   */
  virtual void setup();

//virtual void addArtificialDiff(std::vector<RealVector>& residual) ;

   RealVector& operator() (const std::vector<Framework::State*>& vars,
                            const RealVector& shapeF,
                            const RealMatrix& grad,
                            Framework::GeometricEntity* const geo);


private :

/// temporary storage of the average advection vector
  RealVector m_adv;
  /// storage for the temporary result of the artificual diffusion computation
  RealVector m_tmpAD;

 /// vector to store the accumulated computation of the residual
    RealVector m_v;

    /// matrix to store accumulated computation of the upwind parameter
    RealMatrix m_k;
 
 RealVector theta;

  /// array of cell states
  std::vector<Framework::State*> * m_cellStates;
  /// array of minimum states
  RealVector m_min_states;
  /// array of maximum states
  RealVector m_max_states;
 
  /// option for adding shock detection to the artificial diffusion
  bool m_with_shock_detect;
  
  
   /// artifical viscosity
  CFreal m_theta;
 
  CFreal m_userViscosity;
  
  /// Option for choosing original diffusion factor as defined by Brooks + Hughes (1982 + 1983)
  bool m_with_Hughes; 
  
  RealVector tau;
  RealVector fluxdot;
  RealVector unitVector;

};//End class SUPG_ArtDiffStrategy

    } // End namespace FluctSplit

}// End namespace COOLFluiD
//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SUPG_ArtDiffStrategy_hh
