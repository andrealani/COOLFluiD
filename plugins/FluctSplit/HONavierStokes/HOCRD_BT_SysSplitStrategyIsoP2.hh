#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD_BT_SysSplitStrategyIsoP2_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD_BT_SysSplitStrategyIsoP2_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
#include "FluctSplit/P2Normal.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics { namespace NavierStokes { class EulerTerm; } }
  namespace MathTools { class MatrixInverter; }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy of higher order
/// elements which residual is computed with a contour integral.
/// @author Nadege Villedieu
/// @author Tiago Quintino
class HOCRD_BT_SysSplitStrategyIsoP2 : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD_BT_SysSplitStrategyIsoP2(const std::string& name);

  /// Destructor.
  virtual ~HOCRD_BT_SysSplitStrategyIsoP2();

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

  static void defineConfigOptions(Config::OptionList& options);
protected: // methods

  /// Compute the integral of the fluxes
  void computeHOCurvedFluctuation();
/// Compute the upwind parameter
    ///  kplus and kmin are some output because we want to store them
    ///  For each sub triangle instead of re-computing
  void computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus);
/// Distribute the residual using N scheme
void distributeN(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiN);

  /// Distribute the residual using LDA scheme
void distributeLDA(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiLDA);


void computeBlendingCoeff(CFreal & result);

protected: // data

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket with theta clending coefficient
  Framework::DataSocketSink<CFreal> socket_thetas;

  /// temporary gradient
  RealVector _grad;

  /// ID of the variable chosen for blending coefficient
  CFuint _varID;
/// min theta
  CFreal m_min_theta;

  /// store the thetas
  bool m_store_thetas;

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

  /// Linear shape functions
  CFreal L1, L2, L3;

  /// Coordinates of reference triangle
  RealVector xi_ref;
  RealVector eta_ref;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector facenormal;

  RealVector m_sumKplusU;

  RealVector m_sumKU;

  RealMatrix m_sumKplus;

  RealVector m_uTemp;

  RealVector m_uMin;

  RealMatrix m_tmp;


  /// temporary data for computation of upwind parameters for 1st sub-triangle
  std::vector<RealMatrix*> m_k1Plus;


  /// temporary data for computation of upwind parameters for 2nd sub-triangle
  std::vector<RealMatrix*> m_k2Plus;


  /// temporary data for computation of upwind parameters for 3rd sub-triangle
  std::vector<RealMatrix*> m_k3Plus;


  /// temporary data for computation of upwind parameters for 4th sub-triangle
  std::vector<RealMatrix*> m_k4Plus;


  /// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_k;

/// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_kMin;

  ///  eignvalues
  std::vector<RealVector*>         m_eValues;

  CFuint m_maxNbStatesInCell;

  /// adimensionalized normal vector
  RealVector               m_adimNormal;

  CFreal m_theta1;
  CFreal m_theta2;
  CFreal m_theta3;
  CFreal m_theta4;
  CFreal m_theta;


  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       m_inverter;

  RealMatrix  m_invK;

  RealVector  m_uInflow;
 /// temporary subcell normals
  InwardNormalsData * m_subcell_normals;

  RealMatrix matrix_face_norms;
  RealMatrix matrix_node_norms;

  RealVector vector_face_areas;
  RealVector vector_node_areas;

  ///Functor to compute normals of P2P2 triangle on the fly
  P2Normal m_CP2N;

 /// pmax - pmin
  CFreal _deltaP;

  /// reference length
  CFreal _length;

  /// reference speed
  CFreal _speed;

  /// name of the variable chosen for blending coefficient
  std::string _varName;

  /// temporary physical data array
  RealVector _pData;

  /// use shock-detector on all sub elements
  bool m_use_max_theta;


/// convective term
 Common::SafePtr<Physics::NavierStokes::EulerTerm> _cterm;

/// flag to see if we are dealing with an high order scheme
  bool _isHO;
  /// max number of subcells
  CFuint m_max_nbsubcells;

  /// Shock detector that we want to use
  std::string _sh_detector;

private: // data


}; // class HOCRD_BT_SysSplitStrategyIsoP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_BT_SysSplitStrategyIsoP2_hh
