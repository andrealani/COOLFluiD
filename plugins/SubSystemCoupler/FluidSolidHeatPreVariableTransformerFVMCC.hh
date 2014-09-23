#ifndef COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformerFVMCC_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformerFVMCC_hh

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

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////

/**
 * This class represents a FluidSolidHeatPre transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */

class FluidSolidHeatPreVariableTransformerFVMCC : public PreVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  FluidSolidHeatPreVariableTransformerFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~FluidSolidHeatPreVariableTransformerFVMCC();

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
    // we will communicate the temperature and the heat flux
    // so it returns 2
    return 2; 
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

  /// average State
  RealVector _avState;

  /// array of temporary values
  RealMatrix _values;
  
  /// array of temporary nodal states
  std::vector<RealVector*> _states;

  ///Vector for the gradients
  std::vector<RealVector*> _gradients;

  ///normal to the face
  RealVector m_normal;

  ///Value for h
  CFreal _hConst;

}; // end of class FluidSolidHeatPreVariableTransformerFVMCC

//////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FluidSolidHeatPreVariableTransformerFVMCC_hh
