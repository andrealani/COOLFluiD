#ifndef COOLFluiD_Numerics_FluctSplit_Acoustic_ArtDiffStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_Acoustic_ArtDiffStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"
#include "Framework/DataSocketSink.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a Artificial Diffusion strategy for acoustic simulations
 *
 * @author Lilla Koloszar
 *
 */
class Acoustic_ArtDiffStrategy : public ArtificialDiffusionStrategy {

public: // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Acoustic_ArtDiffStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~Acoustic_ArtDiffStrategy();

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

private :

   /// dimension
   CFuint  m_dim;

  /// artificial numerical viscosity
  CFreal m_nu;

   /// normal of eact state
   std::vector<RealVector> m_stateNormal;

 
};//End class Acoustic_ArtDiffStrategy

    } // End namespace FluctSplit

}// End namespace COOLFluiD
//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_Acoustic_ArtDiffStrategy_hh
