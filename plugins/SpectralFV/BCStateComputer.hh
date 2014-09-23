#ifndef COOLFluiD_Numerics_SpectralFV_BCStateComputer_hh
#define COOLFluiD_Numerics_SpectralFV_BCStateComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the ghost states corresponding
 * to a boundary condition
 *
 * @author Kris Van den Abeele
 */
class BCStateComputer : public SpectralFVMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,BCStateComputer > PROVIDER;

public:  // methods

  /// Constructor
  BCStateComputer(const std::string& name);

  /// Destructor
  ~BCStateComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCStateComputer";
  }

  /// Set up private data and data
  virtual void setup();

  /// function returning a boolean that is true if the boundary condition requires the spatial coordinates
  bool needsSpatialCoordinates()
  {
    return m_needsSpatCoord;
  }

  /// set the global ID of the boundary face
  void setFaceID(const CFuint faceID)
  {
    m_faceID = faceID;
    setBndFaceData();
  }

  /// adds a trs name
  void addTRSName(const std::string trsName)
  {
    m_trsNames.push_back(trsName);
  }

  /// get trs names
  Common::SafePtr< std::vector< std::string > > getTRSNames()
  {
    return &m_trsNames;
  }

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  virtual void computeGhostStates(const std::vector< Framework::State* >& intStates,
                                  std::vector< Framework::State* >& ghostStates,
                                  const std::vector< RealVector >& normals,
                                  const std::vector< RealVector >& coords) = 0;

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
   virtual void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                      std::vector< std::vector< RealVector* > >& ghostGrads,
                                      const std::vector< RealVector >& normals,
                                      const std::vector< RealVector >& coords) = 0;

protected: // methods

  /**
   * set data related to the current boundary face
   */
  virtual void setBndFaceData() {}

protected: // data

  /// boolean telling whether the boundary condition needs the coordinates of the flux points
  bool m_needsSpatCoord;

  /// global ID of the boundary face
  CFuint m_faceID;

  /// list of names of TRSs the BC applies to
  std::vector< std::string > m_trsNames;

}; // class BCStateComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BCStateComputer_hh

