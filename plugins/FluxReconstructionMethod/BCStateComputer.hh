#ifndef COOLFluiD_FluxReconstructionMethod_BCStateComputer_hh
#define COOLFluiD_FluxReconstructionMethod_BCStateComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include <map>

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the ghost states corresponding
 * to a boundary condition
 *
 * @author Kris Van den Abeele
 * @author Ray Vandenhoeck
 * @author Rayan Dhib
 * 
 */
class BCStateComputer : public FluxReconstructionSolverStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider< FluxReconstructionSolverData,BCStateComputer > PROVIDER;

public:  // methods

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  BCStateComputer(const std::string& name);

  /// Destructor
  ~BCStateComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCStateComputer";
  }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Set up private data and data
  virtual void setup();
  
   //Added from FV for parallelization
  /**
   * Set the preProcesses connectivity between faces belonging to different process
   *
   */
  virtual void preProcess();

  
  /// Unset up private data and data
  virtual void unsetup();

  /// function returning a boolean that is true if the boundary condition requires the spatial coordinates
  bool needsSpatialCoordinates()
  {
    return m_needsSpatCoord;
  }

  /// function returning a boolean that is true if the boundary condition requires the extra variables
  bool needsExtraVariables()
  {
    return m_needsExtraVars;
  }
  
  /// set the current face
  void setFace(Framework::GeometricEntity *const face)
  {
      m_face = face;
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

  /// set extra variables in the flux points
  void setExtraVars(std::vector< RealVector* >* extraVars)
  {
    m_extraVars = extraVars;
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

  /**
   * Sets the boundary values g_b of the gradient variables at the flux points of
   * the current boundary face: the interface value g^I_f of a boundary face,
   * which the correction (g_b - g^D_f) grad h_f of this face uses, with g^D_f
   * the gradient variables of the interior cell extrapolated to the flux
   * points. It enters the gradient of the cell corrected with all its faces,
   * the compact gradient the diffusive boundary flux uses, and the derivatives
   * of both in the Jacobian. g_b is a value of the gradient variables, not a
   * state.
   *
   * Default, component by component, with U the interior state at the flux
   * point and g(.) the gradient variables of a state:
   *
   *   g_b = g^D_f + 0.5*(g(U_ghost) - g(U))
   *
   * @param gradVarsFlxPnt  gradient variables extrapolated to the flux points
   * @param intStates       interior states at the flux points
   * @param ghostStates     ghost states at the flux points
   * @param unitNormals     unit normals at the flux points, pointing out of the domain
   * @param flxPntCoords    coordinates of the flux points
   * @param bndGradVars     the boundary values of the gradient variables
   */
  virtual void computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                  const std::vector< Framework::State* >& intStates,
                                  const std::vector< Framework::State* >& ghostStates,
                                  const std::vector< RealVector >& unitNormals,
                                  const std::vector< RealVector >& flxPntCoords,
                                  std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary states at the flux points of the current boundary face:
   * the states the diffusive boundary flux and the transport properties are
   * evaluated at. Default: U_b = 0.5*(U + U_ghost).
   */
  virtual void computeBndStates(const std::vector< Framework::State* >& intStates,
                                const std::vector< Framework::State* >& ghostStates,
                                const std::vector< RealVector >& unitNormals,
                                const std::vector< RealVector >& flxPntCoords,
                                std::vector< RealVector* >& bndStates);

  /**
   * Sets the gradients bndGrads[iFlx][iEq] of the diffusive boundary flux from
   * the compact face gradients intGrads[iFlx][iEq] of the interior cell; the
   * rules of the boundary condition on the normal derivatives are applied here.
   * bndStates are the states of computeBndStates.
   * Default: q_b = 0.5*(q + q_ghost), with q_ghost from computeGhostGradients.
   */
  virtual void computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                               std::vector< std::vector< RealVector* > >& bndGrads,
                               const std::vector< RealVector* >& bndStates,
                               const std::vector< RealVector >& unitNormals,
                               const std::vector< RealVector >& flxPntCoords);

  /**
   * Applies the rules on the boundary gradients of one flux point that need the
   * transport properties. Called per flux point after prepareFluxComputation,
   * on the gradients of computeBndGrads. Default: nothing.
   */
  virtual void constrainBndGrads(const RealVector& bndState,
                                 std::vector< RealVector* >& bndGrads,
                                 const RealVector& unitNormal,
                                 const RealVector& flxPntCoord)
  {
  }

  /**
   * Sets the transition flag of flux point iFlux of the current face (set with
   * setFace). Only the Gamma-Alpha transition model uses it: its wall boundary
   * condition reads it in computeGhostStates through transitionCriterion(). The
   * flags are kept per face, so every later evaluation of the same face reads
   * the flags set for it.
   */
  void setTransitionCriterion(const CFuint iFlux, const bool transition);

protected: // methods

  /**
   * Transition flag of flux point iFlux of the current face, false when no flag
   * was set for the face yet.
   */
  bool transitionCriterion(const CFuint iFlux) const;

  /**
 * Configures the command.
   */
  void configure ( Config::ConfigArgs& args );

  /**
   * Copies the gradients of every flux point: q_b = q.
   * @param intGrads  gradients [iFlx][iEq]
   * @param bndGrads  the copies [iFlx][iEq]
   */
  void copyGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                     std::vector< std::vector< RealVector* > >& bndGrads);

  /**
   * Removes the normal component of a gradient: grad -= (grad.n) n.
   * @param grad    the gradient
   * @param normal  unit normal n
   */
  void removeNormalComponent(RealVector& grad, const RealVector& normal);

  /**
   * Sets the boundary values of the gradient variables of a slip wall: the
   * gradient variables extrapolated to the flux points with the normal velocity
   * removed, g_b = g^D_f - (u.n) n on the velocity components.
   * @param gradVarsFlxPnt  gradient variables extrapolated to the flux points
   * @param unitNormals     unit normals at the flux points
   * @param velocityIDs     indices of the velocity components in the gradient variables
   * @param bndGradVars     the boundary values of the gradient variables
   */
  void setSlipWallBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                              const std::vector< RealVector >& unitNormals,
                              const std::vector< CFuint >& velocityIDs,
                              std::vector< RealVector* >& bndGradVars);

  /**
   * Sets the boundary gradients of a slip wall. On the velocity block J = grad(u)
   *
   *   J_b = (J + R J R)/2,   R = I - 2 n n^T,
   *
   * which removes both mixed normal-tangential blocks; every other gradient
   * loses its normal component.
   * @param intGrads     compact face gradients [iFlx][iEq]
   * @param bndGrads     the boundary gradients [iFlx][iEq]
   * @param unitNormals  unit normals at the flux points
   * @param velocityIDs  indices of the velocity components in the gradient variables
   */
  void setSlipWallBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                           std::vector< std::vector< RealVector* > >& bndGrads,
                           const std::vector< RealVector >& unitNormals,
                           const std::vector< CFuint >& velocityIDs);

private: // methods

  /**
   * Adds curvature to the boundary faces using the domain model
   */
    void addCurvatureToBndFaces();

protected: // data

  /// boolean telling whether the boundary condition needs the coordinates of the flux points
  bool m_needsSpatCoord;

  /// boolean telling whether the boundary condition needs extra variables in the flux points
  bool m_needsExtraVars;

  /// list of names of TRSs the BC applies to
  std::vector< std::string > m_trsNames;

  /// pointer to extra variables in the boundary points
  std::vector< RealVector* >* m_extraVars;

  /// boolean telling whether to use the domain model to add curvature to the boundary faces
  bool m_useDomainModel;
  
  /// variable for current face
  Framework::GeometricEntity* m_face;
  
  /// transition flags of the flux points of every face seen so far, by face ID (only used in Gamma-Alpha)
  std::map< CFuint, std::vector< bool > > m_transitionCriterion;

  /// number of flux points of a face, size of one entry of m_transitionCriterion
  CFuint m_nbrTransitionFlags;

  /// gradient variables of the interior states at the flux points, scratch of the default computeBndGradVars
  RealMatrix m_gradVarsFace;

  /// gradient variables of the ghost states, scratch of the default computeBndGradVars
  RealMatrix m_gradVarsGhost;

  /// state data pointers passed to the diffusive variable set
  std::vector< RealVector* > m_gradVarStatePtrs;

  /// tangential gradient of the normal velocity, scratch of setSlipWallBndGrads
  RealVector m_tangentialGradUn;

}; // class BCStateComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_FluxReconstructionMethod_BCStateComputer_hh
