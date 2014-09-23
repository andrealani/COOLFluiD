#ifndef COOLFluiD_Numerics_FluctSplit_Null_ArtDiffStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_Null_ArtDiffStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/ArtificialDiffusionStrategy.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a Artificial Diffusion strategy
/// @author Nadege Villedieu
class FluctSplit_API Null_ArtDiffStrategy : public ArtificialDiffusionStrategy {

public: // methods

  /// Constructor.
  Null_ArtDiffStrategy(const std::string& name);

  /// Destructor.
  virtual ~Null_ArtDiffStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    ArtificialDiffusionStrategy::configure(args);
  }

  /// Set up private data and data
  virtual void setup();

//virtual void addArtificialDiff(std::vector<RealVector>& residual) ;

  RealVector& operator() (const std::vector<Framework::State*>& vars,
                            const RealVector& shapeF,
                            const RealMatrix& grad,
                            Framework::GeometricEntity* const geo);

  void addArtificialDiff(std::vector<RealVector>& residual) ;
private :

/// temporary storage of the average advection vector
  RealVector m_adv;
  /// storage for the temporary result of the artificual diffusion computation
  RealVector m_tmpAD;

/// vector to store the accumulated computation of the residual
    RealVector m_v;

    /// matrix to store accumulated computation of the upwind parameter
    RealMatrix m_k;

};//End class SUPG_ArtDiffStrategy

    } // End namespace FluctSplit

}// End namespace COOLFluiD
//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SUPG_ArtDiffStrategy_hh
