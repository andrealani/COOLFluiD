#ifndef COOLFluiD_Numerics_FluctSplit_HOCRD3D_hh
#define COOLFluiD_Numerics_FluctSplit_HOCRD3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"
//#include "FluctSplit/P2Normal_Tetra.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace MathTools { class MatrixInverter; }
namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// RDS strategy for 3D Tetras P2
/// @author Antonino Bonanni
/// @author Tiago Quintino
class HOCRD3D : public FluctuationSplitStrategy {

public: // methods

  /// Constructor.
  HOCRD3D(const std::string& name);

  /// Destructor.
  virtual ~HOCRD3D();

// //-----------OLD------------------
//   /// Configure the object
//   virtual void configure(const Config::ConfigArgs& args)
//   {
//     FluctuationSplitStrategy::configure(args);
//   }
// //--------------------------------

//-----------NEW------------------
  /// Configure the object
  virtual void configure(Config::ConfigArgs& args)
  {
    FluctuationSplitStrategy::configure(args);
  }
//--------------------------------



  static void defineConfigOptions(Config::OptionList& options);

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
  void computeHOCurvedFluctuation(CFuint & i1);
  
  /// Compute the upwind parameter
    ///  kplus and kmin are some output because we want to store them
    ///  For each sub triangle instead of re-computing
  void computeK(const std::vector<Framework::State*>& states,
                std::vector<RealMatrix*>& m_kPlus);


  /// Distribute the residual using LDA scheme
void distributeLDA(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiLDA);


    /// Distribute the residual using LDA scheme
void distributeN(std::vector<RealMatrix*> & m_kPlus,std::vector<RealVector> & phiN);


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
  
/// states that compose each sub face
  MathTools::CFMat<CFuint> subelem_nodetable;
  /// states that compose each sub element
  std::vector<Framework::State*> substates;

  /// residuals for each sub element
  std::vector<RealVector> subresidual;

  /// fluxes on each sub face
  std::vector<RealVector> faceflux;

  /// quadrature points per face
  RealVector qd0;
  RealVector qd1;
  RealVector qd2;
  RealVector wqd;

  /// Coordinates of reference triangle
  RealVector xi_ref;
  RealVector eta_ref;
  RealVector zeta_ref;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector Facenormal0;
  RealVector Facenormal1;
  RealVector Facenormal2;
  RealVector Facenormal3;
  RealVector Facenormal;


  /// temporary subcell normals
  InwardNormalsData * m_subcell_normals;

  RealMatrix matrix_face_norms;
  RealMatrix matrix_node_norms;

  RealVector vector_face_areas;
  RealVector vector_node_areas;
 
  /// subfaces normals

  RealVector nf0;
  RealVector nf1;
  RealVector nf2;
  RealVector nf3;
  RealVector nf4;
  RealVector nf5;
  RealVector nf6;
  RealVector nf7;
  RealVector nf8;
  RealVector nf9;
  RealVector nf10;
  RealVector nf11;
  RealVector nf12;
  RealVector nf13;
  RealVector nf14;
  RealVector nf15;
  RealVector nf16;
  RealVector nf17;
  RealVector nf18;
  RealVector nf19;
  RealVector nf20;
  RealVector nf21;

//........only for "independent" computation..
RealVector advection;
RealVector Flux;
RealVector subFluc;
//..........................................

  /// subfaces Areas
  RealVector FacesAreas;

  /// surfaces normals stored in a matrix
  RealMatrix FacesMat;
  
/// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_kMin;
 /// temporary data for computation of upwind parameters for 1st sub-triangle
  std::vector<RealMatrix*> m_kPlus;

  /// temporary data for computation of upwind parameters
  std::vector<RealMatrix*> m_k;


  RealMatrix m_sumKplus;
  ///  eignvalues
  std::vector<RealVector*>         m_eValues;
 
 /// temporary data for holding positive upwind parameter
  RealMatrix                       m_rightEv;
 
 /// temporary data for holding negative upwind parameter
  RealMatrix                       m_leftEv;

 /// temporary data for holding eignvaluesMinus
 RealVector                       m_eValuesP;
 
  /// temporary data for holding eignvaluesPlus
       RealVector                       m_eValuesM;
       

/// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       m_inverter;

  RealMatrix  m_invK;

  RealVector m_uTemp;
  RealVector m_sumKplusU;
  /// adimensionalized normal vector
  RealVector               m_adimNormal;
//  ///Object to compute normals of P2P2 triangle on the fly
//  P2Normal_Tetra  m_CP2N;
  RealVector Fn0u0;
  RealVector Fn1u1;
  RealVector Fn2u2;
  RealVector Fn3u3;
  RealVector Fn0u4;
  RealVector Fn1u4;
  RealVector Fn1u5;
  RealVector Fn2u5;
  RealVector Fn0u6;
  RealVector Fn2u6;
  RealVector Fn1u7;
  RealVector Fn3u7;
  RealVector Fn3u8;
  RealVector Fn2u8;
  RealVector Fn3u9;
  RealVector Fn0u9;


  bool m_order;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _kPlus;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _kMin;

  /// temporary data for computation of upwind parameters
  std::vector<RealVector> _k;

private: // data


  /// the single splitter
  Common::SafePtr<Splitter> m_splitter;
  CFuint CellID;
  CFuint m_scheme;

  enum { nbQdPts = 4 };

}; // class HOCRD_SplitStrategyIsoP2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_HOCRD_SplitStrategyIsoP2_hh
