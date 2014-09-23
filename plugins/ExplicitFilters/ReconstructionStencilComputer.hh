#ifndef COOLFluiD_Numerics_ExplicitFilters_ReconstructionStencilComputer_hh
#define COOLFluiD_Numerics_ExplicitFilters_ReconstructionStencilComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "ExplicitFilters/StencilComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {

      class ReconstructionFilter;
            
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class provides an abstract interface for functors
   * computing the stencil for a specific kind of polynomial
   * reconstruction in ExplicitFilters
   *
   * @author Willem Deconinck
   */
class ReconstructionStencilComputer : public StencilComputer {

public:

  /**
   * Returns the DataSockets that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();

  /**
   * Defines the config options of this class
   * @param   options   config options of this class
   */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Configure the object
   * @param    args   configuration arguments
   */
  virtual void configure(Config::ConfigArgs& args);

  /**
   * Constructor
   */
   ReconstructionStencilComputer(const std::string& name);

  /**
   * Destructor
   */
  virtual ~ReconstructionStencilComputer();

  /**
   * Setup
   */
  virtual void setup();

  /**
   * Post-process a specific stencil
   */
  virtual void postProcessStencil(const CFuint& centreStateID);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ReconstructionStencilComputer";
  }

  CFuint CFround(const CFreal& x) {
    return floor(x+0.5);
  }
  
  CFdouble CFmathLog(const CFreal& base, const CFreal& x) {
    return (log(x)/log(base));
  }

protected: //data
  
private:
  
  Framework::DataSocketSource<std::vector<std::vector<Framework::State*> > > socket_stencilRings;
  
  CFuint m_nbDistanceZones;

}; // end of class ReconstructionStencilComputer

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace ExplicitFilters

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "ReconstructionFilter.hh"

#endif // COOLFluiD_Numerics_ExplicitFilters_ReconstructionStencilComputer_hh
