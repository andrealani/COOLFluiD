#ifndef COOLFluiD_Physics_StructMechHeat_StructMechHeatPhysicalModel_hh
#define COOLFluiD_Physics_StructMechHeat_StructMechHeatPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Framework/VectorialFunction.hh"
#include "StructMech/StructMechPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a StructMechHeatPhysicalModel.
 *
 * @author Thomas Wuilbaut
 *
 */
class StructMechHeatPhysicalModel : public COOLFluiD::Physics::StructMech::StructMechPhysicalModel
{
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor without arguments
   */
  StructMechHeatPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechHeatPhysicalModel();

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("StructMechHeat" + Common::StringOps::to_str(getDimension()) + "D");
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
    std::string result;
    result = getTypeName() + "Diffusive";
    return result;
  }

  /**
   * Get the Linear Source name
   */
  std::string getSourceName() const
  {
    std::string result;
    result = getTypeName() + "Source";
    return result;
  }

  /**
   * Set up
   */
  void setup();

  /**
   * Set the maximum nb of states data
   */
  void setMaxNbStatesData(const CFuint maxNbStatesData)
  {
    throw Common::NotImplementedException (FromHere(),"StructMechHeatPhysicalModelStack::getActive()->setMaxNbStatesData()");
  }

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the conductivity material property
   */
  CFreal getConductivity() const
  {
    cf_assert(m_constantConductivity == true);

    return m_conductivity;
  }

  /**
   * Get the initial Temperature (uniform) for which there is no deformations
   */
  CFreal getInitialTemp() const
  {
    return m_initialTemp;
  }

  /**
   * Get the conductivity material property
   */
  CFreal getConductivity(const RealVector& coord, const RealVector& state);

  /**
   * Set current zone name
   */
  virtual void setCurrentZone(const std::string zoneName);

protected: // data

  /// Conductivity
  CFreal m_conductivity;

  ///initial Temperature (uniform) for which there is no deformations
  CFreal m_initialTemp;

  ///Nb Dim
  CFuint m_nbDim;

  ///Nb Eqs
  CFuint m_nbEqs;

  /// flag to know if conductivity is constant
  bool m_constantConductivity;

  /// output of the function computation
  RealVector m_Computedconductivity;

  /// a vector of string to hold the functions
  std::vector<std::string> m_functions;

  /// a vector of string to hold the functions
  std::vector<std::string> m_vars;

  /// the VectorialFunction to use
  Framework::VectorialFunction m_vFunction;

  /// storage for the temporary boundary point coordinates
  RealVector m_variables;

}; // end of class StructMechHeatPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMechHeat_StructMechHeatPhysicalModel_hh
