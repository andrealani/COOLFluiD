#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD_BT_ScalarSplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD_BT_ScalarSplitStrategy_hh

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
class HOCRD_BT_ScalarSplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD_BT_ScalarSplitStrategy(const std::string& name);

  /// Destructor.
  virtual ~HOCRD_BT_ScalarSplitStrategy();

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

  /// Compute the upwind parameter
    ///  kplus and kmin are some output because we want to store them
    ///  For each sub triangle instead of re-computing
  void computeK(const std::vector<Framework::State*>& states,
                std::vector<RealVector>& m_kPlus,
                std::vector<RealVector>& m_kMin);

  /// Distribute the residual using N scheme
void distributeN(std::vector<RealVector> & m_kPlus,std::vector<RealVector> & phiN);

  /// Distribute the residual using LDA scheme
void distributeLDA(std::vector<RealVector> & m_kPlus,std::vector<RealVector> & phiLDA);


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

  /// fluctuation on each sub element
  std::vector<RealVector> m_phiN1;

  /// fluctuation on each sub element
  std::vector<RealVector> m_phiN2;

  /// fluctuation on each sub element
  std::vector<RealVector> m_phiN3;

  /// fluctuation on each sub element
  std::vector<RealVector> m_phiN4;

/// fluctuation on each sub element
  std::vector<RealVector> m_phi;

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

  /// temporary face normal
  RealVector facenormal;

  RealVector m_sumKplusU;

  RealVector m_sumKplus;

  RealVector m_uTemp;

  RealVector m_uMin;

  RealVector m_temp;


  /// temporary data for computation of upwind parameters for 1st sub-triangle
  std::vector<RealVector> m_k1Plus;

  /// temporary data for computation of upwind parameters for 1st sub-triangle
  std::vector<RealVector> m_k1Min;

  /// temporary data for computation of upwind parameters for 2nd sub-triangle
  std::vector<RealVector> m_k2Plus;

  /// temporary data for computation of upwind parameters for 2nd sub-triangle
  std::vector<RealVector> m_k2Min;

  /// temporary data for computation of upwind parameters for 3rd sub-triangle
  std::vector<RealVector> m_k3Plus;

  /// temporary data for computation of upwind parameters for 3rd sub-triangle
  std::vector<RealVector> m_k3Min;

  /// temporary data for computation of upwind parameters for 4th sub-triangle
  std::vector<RealVector> m_k4Plus;

  /// temporary data for computation of upwind parameters for 4th sub-triangle
  std::vector<RealVector> m_k4Min;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> m_k;


  CFuint m_maxNbStatesInCell;

  /// adimensionalized normal vector
  RealVector               m_adimNormal;

  /// Tmporary vector containing the sum_{i \in Tk} |\Phi_i^N|
  CFreal   m_phitot;

  /// Vector of theta of the first sub-triangle
  RealVector               m_theta1;

/// Vector of theta of the second sub-triangle
  RealVector               m_theta2;

/// Vector of theta of the third sub-triangle
  RealVector               m_theta3;

/// Vector of theta of the fourth sub-triangle
  RealVector               m_theta4;

/// Vector of maximum of the four sub-triangle
  RealVector               m_theta;

private: // data

}; // class HOCRD_BT_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_ScalarBT_SplitStrategy_hh
