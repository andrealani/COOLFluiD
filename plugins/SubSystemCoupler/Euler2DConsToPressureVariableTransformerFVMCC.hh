#ifndef COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformerFVMCC_hh
#define COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformerFVMCC_hh

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
 * This class represents a Euler2DConsToPressure transformer of variables
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DConsToPressureVariableTransformerFVMCC : public Euler2DConsToPressureVariableTransformer {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);


  /**
   * Default constructor without arguments
   */
  Euler2DConsToPressureVariableTransformerFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DConsToPressureVariableTransformerFVMCC();

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

private:

  //Reference pressure (to substract from pressure obtained)
  CFreal _referencePressure;

}; // end of class Euler2DConsToPressureVariableTransformerFVMCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_Euler2DConsToPressureVariableTransformerFVMCC_hh
