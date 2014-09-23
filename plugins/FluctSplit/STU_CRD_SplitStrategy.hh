#ifndef COOLFluiD_Numerics_FluctSplit_STU_CRD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_STU_CRD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

      class SpaceTime_Splitter;

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Nadege Villedieu
/// This class is a special case of STUCRD strategy
/// It is done in order to optimize the computation time
class STU_CRD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  STU_CRD_SplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~STU_CRD_SplitStrategy();

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

  /// Compute the past or present contribution of integral of the Space-Time fluxes
  /// depending on the flag
  void computeSTFluxIntegral(bool isPast);
  
  /// Sets the current cell and calls the computation of the
  /// consistent state transformation.
  void setCurrentCell();

  /// Gqusslobatto quadrature rule
  void gausslobatto( std::vector<Framework::State*> *state, RealVector &m_phi);
  
  /// Past flux computation
  void pastFlux(CFreal Area) 
  {
    // compute the sum_i \int_T N_i u_i
    m_flux = -Area*(*(*m_linearStates)[0]);
    for (CFuint iState = 1 ; iState < _nbStatesInCell; ++ iState) {
      m_flux -= Area*(*(*m_linearStates)[iState]);
    }
  }
  
  /// Present flux computation
  void presentFlux(CFreal Area) 
  {
    // compute flux in time as the sum_i \int_T N_i u_i
    m_flux += Area*(*(*m_linearStates)[0]);
    for (CFuint iState = 1; iState < _nbStatesInCell; ++ iState) {
      m_flux += Area*(*(*m_linearStates)[iState]);
    }
  }
  
private: // data

  /// The socket to use in this strategy for the update coefficients
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
  RealVector m_flux;

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

  ///Temporary vector of the contour integration
  RealVector m_tmp;

  /// Number of states per cell
  CFuint _nbStatesInCell;
  
  /// Number of equations
  CFuint nbEqs;
  
  /// Dimension
  CFuint _dim;

 /// Storage for the residual (inter states contribution for the space residual)
  /// its size is nbcell and nbeqs*nbstaInsubcell
  RealVector m_pastResiduals; 

  /// weight of gauss-lobatto 
  RealVector m_w;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

 /// temporary face normal
  RealVector facenormal;
}; // class STU_CRD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STU_CRD_SplitStrategy_hh
