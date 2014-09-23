#ifndef COOLFluiD_Physics_Heat_HeatPhysicalModel_hh
#define COOLFluiD_Physics_Heat_HeatPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Common/CFMap.hh"

#include "Framework/PhysicalModelImpl.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

    class ConvHeatTerm;
    class DiffHeatTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a HeatPhysicalModel.
 *
 * @author Tiago Quintino
 *
 */
class HeatPhysicalModel : public Framework::PhysicalModelImpl
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
  HeatPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~HeatPhysicalModel();

  /**
   * Get the convective term
   */
  virtual Common::SafePtr<Framework::BaseTerm> getConvectiveTerm() const;

  /**
   * Get the diffusive term
   */
  virtual Common::SafePtr<Framework::BaseTerm> getDiffusiveTerm() const;

  /**
   * Get the source term
   *
   */
  virtual Common::SafePtr<Framework::BaseTerm> getSourceTerm() const
  {
    return CFNULL;
  }

  /**
   * Get the name of the Physical Model type
   * @return name in a std::string
   */
  std::string getTypeName() const
  {
    return std::string("Heat" + Common::StringOps::to_str(getDimension()) + "D");
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
    return "Null";
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
    throw Common::NotImplementedException (FromHere(),"HeatPhysicalModelStack::getActive()->setMaxNbStatesData()");
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
   * Get the conductivity material property
   */
  CFreal getConductivity(const RealVector& coord, const RealVector& state);

  /**
   * Set current zone name
   */
  virtual void setCurrentZone(const std::string zoneName);

  /**
   * This struct groups all the data that is typical of this PhysicalModel
   */
  struct PhysicalData {

    /// average speed
    RealVector avVel;

  }; // end PhysicalData

  /**
   * Get the Physical data
   */
  PhysicalData* getPhysicalData()
  {
    return &m_physicalData;
  }

protected: // data

  /// Conductivity
  CFreal m_conductivity;

  ///Nb Dim
  CFuint m_nbDim;

  ///Nb Eqs
  CFuint m_nbEqs;

  /// physical quantities to be computed in a linearized state
  PhysicalData m_physicalData;

  /// convective term
  std::auto_ptr<ConvHeatTerm> m_convectiveTerm;

  /// diffusive term
  std::auto_ptr<DiffHeatTerm> m_diffusiveTerm;

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

  /// names of the various zones of the mesh
  std::vector<std::string> m_zoneNames;

  /// map between zone name and a zoneID
  Common::CFMap<std::string, CFuint> m_zonesMap;

  /// zoneID of the current zone
  CFuint m_zoneID;

}; // end of class HeatPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_HeatPhysicalModel_hh
