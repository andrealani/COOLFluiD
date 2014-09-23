#ifndef COOLFluiD_Numerics_FiniteVolume_RoeFluxALEBDF2_hh
#define COOLFluiD_Numerics_FiniteVolume_RoeFluxALEBDF2_hh

//////////////////////////////////////////////////////////////////////////////

#include "RoeFluxALE.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the Roe flux for Arbitrary Langragian Eulerian formulation
 * using the Three-Point Backward Difference time integration and geometric parameter
 * averaging (see I.Lepot thesis, section 5.4)
 *
 * @author Thomas Wuilbaut
 *
 */

class RoeFluxALEBDF2 : public RoeFluxALE {
public:

  /**
   * Constructor
   */
  RoeFluxALEBDF2(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~RoeFluxALEBDF2();
 
  /** 
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data
   */
  virtual void setup()
  {
    RoeFluxALE::setup();

    _pastUnitNormal.resize(Framework::PhysicalModelStack::getActive()->getDim());
  }
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:

  /**
   * Set the abs of the eigen values
   */
  virtual void setAbsEigenValues();
  
private:
  
  /// storage of the past nodes
  Framework::DataSocketSink <Framework::Node*> socket_pastPastNodes;
  
  /// storage of the averaged Normals 0.5*(n + n+1)
  Framework::DataSocketSink<CFreal> socket_avNormals;
  
  /// temporary past unit normal 0.5*(n + n-1)
  RealVector   _pastUnitNormal;
  
  /// use the shape functions
  bool _useShapeFunctions;
  
}; // end of class RoeFluxALEBDF2

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_RoeFluxALEBDF2_hh
