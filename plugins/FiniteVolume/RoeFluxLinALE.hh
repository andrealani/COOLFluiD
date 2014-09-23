#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFluxLinALE_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFluxLinALE_hh

//////////////////////////////////////////////////////////////////////////////

#include "BaseRoeFlux.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/JacobianLinearizer.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux for Arbitrary Langragian Eulerian formulation
 * using the Cranck-Nicholson time integration and geometric parameter
 * averaging (see I.Lepot thesis, section 5.4)
 *
 * @author Thomas Wuilbaut
 *
 */
class RoeFluxLinALE : public BaseRoeFlux {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  RoeFluxLinALE(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeFluxLinALE();

  /**
   * Set up private data
   */
  virtual void setup()
  {
    CFAUTOTRACE;

    _sumFlux.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _result.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _rightEv.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
      Framework::PhysicalModelStack::getActive()->getNbEq());
    _leftEv.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
      Framework::PhysicalModelStack::getActive()->getNbEq());
    _eValues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _rightEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _leftEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _absEvalues.resize(Framework::PhysicalModelStack::getActive()->getNbEq());
    _absJacobStar.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
              Framework::PhysicalModelStack::getActive()->getNbEq());
    _jacobStar.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
              Framework::PhysicalModelStack::getActive()->getNbEq());
    _tempUnitNormal.resize(Framework::PhysicalModelStack::getActive()->getDim());
    _matResult.resize(Framework::PhysicalModelStack::getActive()->getNbEq(),
              Framework::PhysicalModelStack::getActive()->getNbEq());
    _meshSpeed.resize(Framework::PhysicalModelStack::getActive()->getDim());

    ///@todo this should disappear!!! Making BaseRoeFlux a MethodStrategy
    std::string dataHandleName = _namespace + "_pastNodes";
    _pastNodes = Framework::MeshDataStack::getActive()->getDataStorage()->getData<Framework::Node*>(dataHandleName);

    dataHandleName = _namespace + "_futureNodes";
    _futureNodes = Framework::MeshDataStack::getActive()->getDataStorage()->getData<Framework::Node*>(dataHandleName);
  }

  /**
   * Set the jacobian linearizer
   */
  void setJacobLinearizer
  (Common::SafePtr<Framework::JacobianLinearizer> linearizer)
  {
    _linearizer = linearizer;
  }

  /**
   * Set the variable set transformer (from solution to linearization
   * variables)
   */
  void setSolutionToLinearVarSetTransformer
  (Common::SafePtr<Framework::VarSetTransformer> transformer)
  {
    _solutionToLinearVarTrans = transformer;
  }

  /**
   * Set the variable set transformer (from update to solution
   * variables)
   */
  void setUpdateToSolutionVarSetTransformer
  (Common::SafePtr<Framework::VarSetTransformer> transformer)
  {
    _updateToSolutionVarTrans = transformer;
  }

  /**
   * Set the solution variable set
   */
  void setSolutionVarSet(Common::SafePtr<Framework::ConvectiveVarSet> solutionVarSet)
  {
    _solutionVarSet = solutionVarSet;
  }

  /**
   * Set the update variable set
   */
  void setUpdateVarSet(Common::SafePtr<Framework::ConvectiveVarSet> updateVarSet)
  {
    _updateVarSet = updateVarSet;
  }

  /**
   * Set the  diffusive variable set
   */
  virtual void setDiffusiveVarSet(Common::SafePtr<Framework::DiffusiveVarSet> diffusiveVarSet)
  {
    _diffusiveVarSet = diffusiveVarSet;
  }

  /**
   * Overloading of operator()
   */
  RealVector& compute();

  /**
   * Overloading of operator()
   */
  RealMatrix& computeJacobian();

protected:

  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues(Framework::FluxSplitterData& data);

  /**
   * Compute the entropy-corrected lambda
   */
  void computeLambdaCorr(CFreal& lambdaCorr) const
  {
    const CFuint nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    for (CFuint i = 0; i < nbEqs; ++i) {
      lambdaCorr = std::max(lambdaCorr,
       std::abs(_rightEvalues[i] - _leftEvalues[i]));
    }
  }

protected: //data

  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _solutionVarSet;

  /// acquaintance of the concrete variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVarSet;

  /// acquaintance of the concrete diffusive variable set
  Common::SafePtr<Framework::DiffusiveVarSet> _diffusiveVarSet;

  /// jacobian linearizer
  Common::SafePtr<Framework::JacobianLinearizer> _linearizer;

  /// variable set transformer from solution to linearization variables
  Common::SafePtr<Framework::VarSetTransformer> _solutionToLinearVarTrans;

  /// variable set transformer from update to solution variables
  Common::SafePtr<Framework::VarSetTransformer> _updateToSolutionVarTrans;

  /// array storing the sum of the right and left flux
  RealVector   _sumFlux;

  /// matrix of right eigenvectors
  RealMatrix   _rightEv;

  /// matrix of left eigenvectors
  RealMatrix   _leftEv;

  /// vector of eigenvalues
  RealVector   _eValues;

  /// vector of right state eigenvalues
  RealVector   _rightEvalues;

  /// vector of left state eigenvalues
  RealVector   _leftEvalues;

  /// vector of eigenvalues
  RealVector   _absEvalues;

  /// ID of the entropy correction type
  CFuint _entropyFixID;

  /// the jacobian matrix
  RealMatrix   _jacobStar;

  /// abs of the jacobian matrix
  RealMatrix   _absJacobStar;

  /// temporary unit normal
  RealVector   _tempUnitNormal;

  /// Mesh Speed
  RealVector   _meshSpeed;

  /// Mesh Speed projected on the unit normal
  CFreal   _vgn;

  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStates;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLR;

  /// pointer to temporary states in solution variables
  std::vector<Framework::State*>* _solutionStatesStar;

  /// vector storing the left and right states of a face
  std::vector<Framework::State*> _statesLRStar;

  /// array storing the temporary solution
  RealMatrix    _matResult;

  /// storage of the Past Nodes (at time n)
  Framework::DataHandle<Framework::Node*> _pastNodes;

  /// storage of the Future Nodes (at time n+1)
  Framework::DataHandle<Framework::Node*> _futureNodes;

}; // end of class RoeFluxLinALE

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFluxLinALE_hh
