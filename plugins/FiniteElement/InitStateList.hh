#ifndef COOLFluiD_Numerics_FiniteElement_InitStateList_hh
#define COOLFluiD_Numerics_FiniteElement_InitStateList_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/VectorialFunction.hh"
#include "FiniteElementMethodData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a initalizing solution command
 *
 * @author Tiago Quintino
 * @author Andrea Lani
 */
class InitStateList : public FiniteElementMethodCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit InitStateList(const std::string& name);

  /**
   * Destructor.
   */
  ~InitStateList();


  /**
   * Execute Processing actions
   */
  virtual void execute();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this command after the processing phase
   */
  void unsetup() {}

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

protected:


  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// a vector of string to hold the functions
  std::vector<std::string> _functions;

  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  // the VectorialFunction to use
  Framework::VectorialFunction _vFunction;

  /// socket for State's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  ///file name
  std::string m_filename;

  ///List of states to apply the DirichletBC
  std::vector<CFuint> m_statesList;

}; // class InitStateList

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_InitStateList_hh

