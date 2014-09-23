#ifndef COOLFluiD_Numerics_SubSystemCoupler_FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC_hh

//////////////////////////////////////////////////////////////////////

#include "PreVariableTransformer.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace Framework {
    class PhysicalChemicalLibrary;
  }

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a FluidSolidStructMechHeatPre transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */

class FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC : public PreVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC();

  /**
   * Configuration
   */
  void configure ( Config::ConfigArgs& args )
  {
    PreVariableTransformer::configure(args);
  }

  /**
   * Returns the DataSocket's that this strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Un Setup the object
   */
  virtual void unsetup();

  /**
   * Set nb of states to transform
   */
  virtual void setNbStates(const CFuint nbOtherStates)
  {
    _nbOtherStates = nbOtherStates;
    if(_pastFluxes.size() == 0) _pastFluxes.resize(_nbOtherStates);
    if(_pastTemperatures.size() == 0) _pastTemperatures.resize(_nbOtherStates);
  }

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original);

  /**
   * Transform a vector
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
                                   const RealVector& coord,
                                   const RealVector& original,
                                   const RealVector& shapeFunctions)
  {
    return preTransform(faces, coord, original);
  }

  /**
   * Return the size of the transformed vector
   */
  CFuint getTransformedSize(const CFuint size)
  {
    return 3;
  }

private:

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;

  /// socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;

  /// socket to the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// socket to the face normals
  Framework::DataSocketSink<CFreal> socket_normals;

  /// socket to the isOutward flag
  Framework::DataSocketSink<CFint> socket_isOutward;

  /// corresponding diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVarSet;

  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;

  ///Vector of the fluxes of the previous iteration
  RealVector _pastFluxes;

  ///Vector of the temperatures of the previous iteration
  RealVector _pastTemperatures;

  ///Value for h
  CFreal _hConst;

  /// flag to know if the value of h is optimized or taken as a constant
  bool _isOptimizeH;

  /// flag to know if radiative flux is taken into account
  bool _isRadiativeFlux;

  /// Characteristic value of the body (black-body=0)
  CFreal _radiativeFluxEpsilon;

  ///Ref pressure
  CFreal _referencePressure;

  /// average State
  RealVector _avState;

  /// array of temporary values
  RealMatrix _values;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;
  
  ///Vector for the gradients
  std::vector<RealVector*> _gradients;

  ///normal to the face
  RealVector _normal;

}; // end of class FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FluidSolidStructMechHeatPreVariableTransformerRhoivtLTEFVMCC_hh
