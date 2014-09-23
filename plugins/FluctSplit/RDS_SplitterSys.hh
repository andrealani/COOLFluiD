#ifndef COOLFluiD_Numerics_FluctSplit_RDS_SplitterSys_hh
#define COOLFluiD_Numerics_FluctSplit_RDS_SplitterSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "Splitter.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Common/CFMultiMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class MatrixInverter; }

  namespace Framework { class GeometricEntity; }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a generic interface for RDS splitters for
/// system schemes
/// @author Andrea Lani
/// @author Tiago Quintino
class RDS_SplitterSys : public Splitter {

public:

  /// Constructor
  /// @see Splitter()
  RDS_SplitterSys(const std::string& name);

  /// Default destructor
  virtual ~RDS_SplitterSys();

  /// Set up
  virtual void setup();

  /// Compute the inflow parameters
  /// @post K+ and K- will be computed
  virtual void computeK(const std::vector<Framework::State*>& states,
			const InwardNormalsData* const normalsData);

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
private: // method
  
  /// Helper function just to compute K and be able to reuse the algorithm
  virtual void doComputeK(CFuint iState);
  
  /// Sets the correct block limits for a System Splitter
  /// Called by the RDS_Splitter constructor.
  /// @see _nbEquations
  /// @see _firstVarID
  /// @see _lastVarID
  virtual void setBlockData();
  
  /// Apply the carbuncle fix
  void applyCarbuncleFix();

  /// Apply carbuncle fix
  void applyCarbuncleFix(const std::vector<Framework::State*>& states,
  		 const InwardNormalsData* const normalsData);

  /// Update the correction for the carbuncle fix
  void updateDelta(Framework::GeometricEntity* const cell,
    bool reuseNShock);

  /// Compute the correction for the carbuncle fix within the current cell
  void computeCellDelta(Framework::GeometricEntity* const cell,
  CFint iState);

  /// Get the neighbor cell ID and update the delta considering the
  /// neighbor contribution
  void updateDelta(CFuint iState);

protected: // data
  
  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  ///  eignvalues
  std::vector<RealVector*>         m_eValues;
  
  /// temporary data for holding positive upwind parameter
  RealMatrix                       _tempKp;
  
  /// temporary data for holding negative upwind parameter
  RealMatrix                       _tempKm;

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       _inverter;

  /// positive upwind parameters
  std::vector<RealMatrix*>        _kPlus;

  /// negative upwind parameters
  std::vector<RealMatrix*>        _kMin;

  /// coefficient depending on the dimension
  CFreal                          m_kCoeff;

  /// temporary storage of nodal area
  CFreal                          m_nodeArea;
  
  /// maximum number of states in cell
  CFuint m_maxNbStatesInCell;

  /// array storing the local IDs in the 2D cell
  /// corresponding to the two nodes defining the shock face
  std::vector<CFuint> m_shockFaceNodes;

  /// corrections for carbuncle fix
  RealVector m_delta;

  /// cell builder
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> m_cellBuilder;

  /// node to cell connecivity
  Common::CFMultiMap<CFuint,CFuint> m_mapNode2CellID;

  CFreal m_nxShock;
  CFreal m_nyShock;
  std::vector<RealMatrix> m_lambda;
  std::vector<RealVector> m_nShock;
  /// temporary data for holding eignvaluesMinus
  RealVector                       m_eValuesP;

  /// temporary data for holding eignvaluesPlus
  RealVector                       m_eValuesM;

  /// temporary data for holding positive upwind parameter
  RealMatrix                       m_rightEv;

  /// temporary data for holding negative upwind parameter
  RealMatrix                       m_leftEv;
  
}; // end of class RDS_SplitterSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_RDS_SplitterSys_hh
