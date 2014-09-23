#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy of higher order
/// elements which residual is computed with a contour integral.
/// @author Nadege Villedieu
/// @author Tiago Quintino
class HOCRD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD_SplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~HOCRD_SplitStrategy();

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
  void computeHOFluctuation();

protected: // data

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// solution variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_solutionVar;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;

  /// fluctuation on each sub element
  std::vector<RealVector*> m_phisubT;

  /// direction of the faces
  MathTools::CFMat<CFreal> subelemfacedir;

  /// faces that compose each sub element
  MathTools::CFMat<CFuint> subelemtable;

  /// states that compose each sub face
  MathTools::CFMat<CFuint> subfacetable;
  
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

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector facenormal;

private: // data

  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;

  /// temporary storage of the fluctuation
  /// which is a contour integral of the fluxes
  RealVector m_phiT;

}; // class HOCRD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategy_hh
