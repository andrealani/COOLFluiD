#ifndef COOLFluiD_Numerics_FiniteElement_ImplicitNewmarkComputeTimeResidual_hh
#define COOLFluiD_Numerics_FiniteElement_ImplicitNewmarkComputeTimeResidual_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteElementMethodData.hh"
#include "Common/CFMap.hh"
#include "Framework/State.hh"
#include "Framework/DataSocketSink.hh"
#include "FElemTypeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// sent to Domain to be executed in order to ComputeTimeResidual the MeshData.
class ImplicitNewmarkComputeTimeResidual : public FiniteElementMethodCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit ImplicitNewmarkComputeTimeResidual(const std::string& name);

  /// Destructor.
  virtual ~ImplicitNewmarkComputeTimeResidual();

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /// Execute Processing actions
  virtual void executeOnTrs();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  ///  Gets Alpha
  CFreal getAlpha() {  return _alpha; }

  ///  Gets Gamma
  CFreal getGamma() {  return _gamma; }

private: //functions

  /// Compute the residual
  void computeElementTimeResidual(Framework::GeometricEntity& cell,
                                  RealMatrix& elemMat,
                                  RealVector& elemVec,
                                  std::vector<RealVector>& residual);

  /// Add compute the term to add in the jacobian
  void computeJacobianTerm();

  /// Clean the perturbed residual vector
  void cleanOtherResidual()
  {
    fill(_otherResidual.begin(),_otherResidual.end(),0.0);
  }

private: // data

  /// Storage of the Alpha Coefficient
  CFreal _alpha;

  /// Storage of the Gamma Coefficient
  CFreal _gamma;

  ///flag for using the row-sum mass matrix lumping
  bool _lumpMassMatrix;

  /// Storage of the newmark coeficients
  CFreal _a3;
  CFreal _a4;
  CFreal _a5;

/// map of LSSMatrix accumulators, one for each cell type
  Common::CFMap<CFuint,FElemTypeData> m_map_femdata;

  /// socket for Rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for Rhs
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// socket for Rhs
  Framework::DataSocketSink<Framework::State*> socket_pastStatesD;

  /// socket for Rhs
  Framework::DataSocketSink<Framework::State*> socket_pastStatesD2;

  //direct datahandles (to avoid getting them at each residual computation
  Framework::DataHandle<Framework::State*> _pastStates;
  Framework::DataHandle<Framework::State*> _pastStatesD;
  Framework::DataHandle<Framework::State*> _pastStatesD2;

  Common::SafePtr<ComputeInertiaTerm> _inertiaTerm;

  /// Temporary storage of the integration result
  RealMatrix _integResultMat;

  /// storage for the temporary node residuals
  RealVector _tempRes;

  /// storage for the temporary perturbed states
  std::vector<RealVector> _otherResidual;

}; // class ImplicitNewmarkComputeTimeResidual

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ImplicitNewmarkComputeTimeResidual_hh

