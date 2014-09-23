#ifndef COOLFluiD_Physics_Heat_DiffHeatTerm_hh
#define COOLFluiD_Physics_Heat_DiffHeatTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the interface for an Heat convective
 * physical term
 *
 * @author Tiago Quintino
 *
 */
class DiffHeatTerm : public Framework::BaseTerm
{
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Enumerator defining the mapping between
   * the variable name and its position in the
   * physical data
   */
  enum { T=0 };

  /**
   * Constructor without arguments
   */
  DiffHeatTerm(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~DiffHeatTerm();

  /**
   * Configures this object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set physical data
   */
  virtual void setupPhysicalData();

  /**
   * Physical data size
   */
  virtual CFuint getDataSize() const
  {
    return 1;
  }

protected: // data

}; // end of class DiffHeatTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_DiffHeatTerm_hh
