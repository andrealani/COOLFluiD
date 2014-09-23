#ifndef COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API CRD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  CRD_SplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~CRD_SplitStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data and data
  virtual void setup();

  /// Unsetup private data and data
  virtual void unsetup();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);
  
protected: // methods
  
  /// Compute the integral of the fluxes
  void computeFluxIntegral();

  /// Sets the current cell and calls the computation of the
  /// consistent state transformation.
  virtual void setCurrentCell();

protected: // data

  /// temporary storage of the fluctuation
  /// which is a contour integral of the fluxes
  RealVector m_phiT;

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;
  
  /// contour integrator
  Common::SafePtr<Framework::ContourIntegrator> m_contourIntegrator;
    
  /// interpolated states at quadrature points
  std::vector<Framework::State*> m_qdstates;

  /// interpolated extra variables at quadrature points
  std::vector<RealVector*> m_qExtraVars;

  /// array for the number of quadrature points in each cell
  /// @todo maybe this could be handled on th side of the integrator itself
  std::vector<CFuint> m_nbQPointsInCell;

  /// back up of update cell states
  std::vector<Framework::State*> m_statesBkp;

}; // class CRD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_CRD_SplitStrategy_hh
