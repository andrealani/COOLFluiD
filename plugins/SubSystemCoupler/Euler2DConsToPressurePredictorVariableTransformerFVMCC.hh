#ifndef COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressurePredictorVariableTransformerFVMCC_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressurePredictorVariableTransformerFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Euler2DConsToPressureVariableTransformer.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler2DConsToPressurePredictor transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DConsToPressurePredictorVariableTransformerFVMCC : public Euler2DConsToPressureVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);


  /**
   * Default constructor without arguments
   */
  Euler2DConsToPressurePredictorVariableTransformerFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DConsToPressurePredictorVariableTransformerFVMCC();

  /**
   * Sets Up the object
   */
  virtual void setup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = Euler2DConsToPressureVariableTransformer::needsSockets();

    result.push_back(&socket_nodes);
    result.push_back(&socket_states);
    result.push_back(&socket_gstates);

    return result;
  }

  /**
   * Transform a state into another one
   */
  virtual RealVector* preTransform(const std::vector<GeoEntityIdx>& faces,
		                  const RealVector& coord,
		                  const RealVector& original);

  /**
   * Transform a vector into another one
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
    return 2;
  }

private:

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;

  //Reference pressure (to substract from pressure obtained)
  CFreal _referencePressure;

  std::vector<CFreal> _pastPressures;

}; // end of class Euler2DConsToPressurePredictorVariableTransformerFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressurePredictorVariableTransformerFVMCC_hh
