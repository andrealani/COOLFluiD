#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyIsoP3_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyIsoP3_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/P3Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy of higher order
/// elements which residual is computed with a contour integral.
/// @author Martin Vymazal
/// @author Nadege Villedieu
/// @author Tiago Quintino

class HOCRD_SplitStrategyIsoP3 : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD_SplitStrategyIsoP3(const std::string& name);

  /// Destructor.
  virtual ~HOCRD_SplitStrategyIsoP3();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /// Set up private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Compute the fluctuation
  /// @param residual the residual for each variable to distribute in each state
  virtual void computeFluctuation(std::vector<RealVector>& residual);

protected: // methods

  /// Compute the integral of the fluxes
  void computeHOCurvedFluctuation();

protected: // data

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// solution variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_solutionVar;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;

  /// matrix storing the scaled face normals generated on the fly
  std::vector<RealVector> m_scaledFaceNormals;

  /// fluctuation on each sub element
  std::vector<RealVector*> m_phisubT;

  /// direction of the faces
  MathTools::CFMat<CFreal> subelemfacedir;

  /// faces that compose each sub element
  MathTools::CFMat<CFuint> subelemtable;

  /// states that compose each sub face
  MathTools::CFMat<CFuint> subfacetable;

  ///to reconstruct subtriangle from subfaces, we need to 
  ///pick the left or right vertex of each subface, depending
  ///on orientation of subfaces
  CFuint m_subfaceVertex;

  /// states that compose each sub element
  std::vector<Framework::State*> substates;

  /// residuals for each sub element
  std::vector<RealVector> subresidual;

  /// fluxes on each sub face
  std::vector<RealVector> faceflux;

  /// quadrature points per face
  RealVector qd0;
  RealVector qd1;
  RealVector wqd;

  /// Linear shape functions
//   CFreal L0, L1, L2;

  /// Coordinates of reference triangle
  RealVector xi_ref;
  RealVector eta_ref;

  ///Number of quadrature points
  enum { nbQdPts = 5 };

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector facenormal;

  /// temporary subcell normals
  InwardNormalsData * m_subcell_normals;

  RealMatrix matrix_face_norms;
  RealMatrix matrix_node_norms;

  RealVector vector_face_areas;
  RealVector vector_node_areas;

  ///Object to compute normals of P3P3 triangle on the fly
  P3Normal m_CP3N;

private: // data

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

}; // class HOCRD_SplitStrategyIsoP3

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyIsoP3_hh
