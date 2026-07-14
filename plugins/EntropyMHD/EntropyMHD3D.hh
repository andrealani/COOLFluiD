#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3D_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectionPM.hh"
#include "EntropyMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * EntropyMHD 3D MHD physical model with GLM divergence cleaning.
 *
 * @author Rayan Dhib
 */
template <int DIM>
class EntropyMHD3D : public Framework::ConvectionPM<EntropyMHDTerm> {

public:

  /**
   * Constructor without arguments
   */
  EntropyMHD3D(const std::string& name);

  /**
   * Default destructor
   */
  ~EntropyMHD3D();

  /**
   * @return the space dimension of the SubSystem
   */
  CFuint getDimension() const;

  /**
   * @return the number of equations of the SubSystem
   */
  CFuint getNbEquations() const;

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("EntropyMHD" + Common::StringOps::to_str(DIM) + "D");
  }

  /**
   * Get the convective name
   */
  std::string getConvectiveName() const
  {
    return getTypeName();
  }

  /**
   * Get the diffusive name
   */
  std::string getDiffusiveName() const
  {
    return "Null";
  }

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

private:

  virtual void setReferenceValues() {}

  virtual void setReferenceTime()
  {
    _refTime = getRefLength();
  }

}; // end of class EntropyMHD3D

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

#include "EntropyMHD3D.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3D_hh
