#ifndef COOLFluiD_Numerics_FluctSplit_STM_CRD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_STM_CRD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {
      class SpaceTime_Splitter;

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Nadege Villedieu
class STM_CRD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  STM_CRD_SplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~STM_CRD_SplitStrategy();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data and data
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this strategy after the  processing phase
  virtual void unsetup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);
  
protected: // methods

  /// Compute the present contribution of integral of the Space-Time fluxes

  void computeSTFluxIntegral();

  /// Compute the past contribution of integral of the Space-Time fluxes

  void computeSTFluxIntegral_past();

  /// Sets the current cell and calls the computation of the
  /// consistent state transformation.
  void setCurrentCell();

private: // data

  /// The socket to use in this strategy for the update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// The socket to use in this strategy for the past normals
  Framework::DataSocketSink<InwardNormalsData*> socket_pastNormals;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// The socket for past volumes of the cells
  Framework::DataSocketSink<CFreal> socket_pastCellVolume;

  /// The socket for the speed of the cells
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

  /// the single splitter
  Common::SafePtr<SpaceTime_Splitter> m_splitter;

  /// temporary storage of the past states
  std::vector<Framework::State*> _pastStates;

  /// temporary storage of the consistent past states
  std::vector<Framework::State*> * tPastStates;

  ///Temporary vector
  RealVector _null;

  /// Contribution from the past to the residual of the cell
  RealVector m_flux_past_time;

  /// Contribution from the present to the residual of the cell
  RealVector m_flux_time;

  /// Contribution from the past to the residual of the cell
  RealVector m_flux_past_space;

  /// Contribution from the present to the residual of the cell
  RealVector m_flux_space;

  /// Volume of the current cell
  CFreal CellVolume;

/// interpolated states at quadrature points
  std::vector<Framework::State*> m_qdstates;

/// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// contour integrator
  Common::SafePtr<Framework::ContourIntegrator> m_contourIntegrator;

  /// array for the number of quadrature points in each cell
  /// @todo maybe this could be handled on th side of the integrator itself
  std::vector<CFuint> m_nbQPointsInCell;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;

  /// back up of update cell states
  std::vector<Framework::State*> m_statesBkp;

  ///Temporary vector of the contour integration
  RealVector m_phi;

  /// Storage for the residual (past states contribution)
  std::vector<RealVector> _pastResiduals;

  /// Storage for the residual (past states contribution of the N scheme used in B scheme)
  std::vector<RealVector> _pastResiduals_order1;

  /// temporary vector with past source term residual
  std::vector<RealVector> temp_residual;


}; // class STU_CRD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STM_CRD_SplitStrategy_hh
