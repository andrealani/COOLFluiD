#ifndef COOLFluiD_Numerics_FiniteVolume_AUSMFluxALE_hh
#define COOLFluiD_Numerics_FiniteVolume_AUSMFluxALE_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an extension of the AUSM scheme suitable for unsteady
 * cases with moving meshes
 *
 * @author Thomas Wuilbaut
 * @author Andrea Lani
 *
 */
template <class BASE>
class AUSMFluxALE : public BASE {
public: // classes
  
  /**
   * Constructor
   */
  AUSMFluxALE(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~AUSMFluxALE();
  
  /**
   * Set up private data to prepare the simulation
   */
  virtual void setup();
  
  /**
   * Compute the mesh speed
   */
  virtual void computeMeshSpeed();
  
  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// Correct the given Mach number
  CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach - m_vgn/this->m_a12;
  }
  
  /**
   * Compute the flux
   */
  virtual void computeMassFluxImpl(RealVector& result);
  
private:
  
  /// storage of the Past Nodes (at time n)
  Framework::DataSocketSink<Framework::Node*> socket_pastNodes;
  
  /// storage of the Future Nodes (at time n+1)
  Framework::DataSocketSink<Framework::Node*> socket_futureNodes;
  
  /// normal mesh speed
  CFreal m_vgn;
  
  /// mesh speed
  RealVector m_meshSpeed;
  
}; // end of class AUSMFluxALE

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "AUSMFluxALE.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMFluxALE_hh
