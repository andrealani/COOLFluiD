#ifndef COOLFluiD_Physics_StructMech_StructMechPhysicalModel_hh
#define COOLFluiD_Physics_StructMech_StructMechPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NotImplementedException.hh"
#include "Framework/PhysicalModelImpl.hh"

#include "StructMech/Materials.hh"
#include "StructMech/MaterialPropertyLib.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

    class StructMechTerm;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for a StructMechPhysicalModel.
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 *
 */
class StructMechPhysicalModel : public Framework::PhysicalModelImpl
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
  StructMechPhysicalModel(const std::string& name);

  /**
   * Default destructor
   */
  ~StructMechPhysicalModel();

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
    return std::string("StructMech" + Common::StringOps::to_str(getDimension()) + "D");
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
    throw Common::NotImplementedException (FromHere(),"StructMechPhysicalModelStack::getActive()->setMaxNbStatesData()");
  }

  /**
   * Configures this configurable object.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Get the Stiffness Matrix
   */
  RealMatrix* getStiffnessMatrix();

  /**
   * Get the Density
   */
  CFreal getDensity() const
  {
    return m_material[m_zoneID]->computeDensity();
  }

  /**
   * Get the Density
   */
  bool isAnisotropic() const
  {
    return m_material[m_zoneID]->isAnisotropic();
  }

  /**
   * Get the Young Modulus
   */
  CFreal getYoung() const
  {
    return m_material[m_zoneID]->computeYoungModulus();
  }

  /**
   * Get the Thermal Expansion coeficient
   */
  CFreal getThermalExpansionCoef()
  {
    return m_material[m_zoneID]->computeThermalExpansionCoef();
  }

  /**
   * Get the Poisson coef
   */
  CFreal getPoisson()
  {
    return m_material[m_zoneID]->computePoissonCoef();
  }

  void setCurrentZone(const std::string zoneName)
  {
    PhysicalModelImpl::setCurrentZone(zoneName);

    m_zoneID = m_zonesMap.find(zoneName);
  }

private: // functions

  /**
   * Transform the Stiffness Matrix from local(material) to global coordinates
   */
  void stiffnessMatrixTransform();

protected: // data

  /// Materialsdynamic viscosity law
  std::vector< Common::SelfRegistPtr<MaterialPropertyLib> > m_material;

  ///Configuration String
  std::vector<std::string>  m_materialStr;

  /// names of the various zones of the mesh
  std::vector<std::string> m_zoneNames;

  /// map between zone name and a zoneID
  Common::CFMap<std::string, CFuint> m_zonesMap;

  /// zoneID of the current zone
  CFuint m_zoneID;

  /// the term used in this physical moder
  std::auto_ptr<StructMechTerm> m_structTerm;


}; // end of class StructMechPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_StructMechPhysicalModel_hh
