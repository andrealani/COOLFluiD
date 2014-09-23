#ifndef COOLFluiD_Numerics_FluctSplit_StrongSubFunctionInletLinEuler3DCons_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSubFunctionInletLinEuler3DCons_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler3DCons.hh"
#include "MathTools/CFMat.hh"

#include "FluctSplit/SpaceTime_Splitter.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace MathTools { class MatrixInverter; }

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * Multidimensional characteristic bc for characteristic variables
 *
 * @author Lilla Koloszar
 *
 *
 *
 */

class StrongSubFunctionInletLinEuler3DCons : public FluctuationSplitCom {
public:

  /**
   * Constructor
   */
  StrongSubFunctionInletLinEuler3DCons(const std::string& name);

  /**
   * Default destructor
   */
  ~StrongSubFunctionInletLinEuler3DCons();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  void unsetup();

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure (Config::ConfigArgs& args);

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private:


/////////////////////// SANDBOX //////////////////////////////////////////////////////////////////////////

  /// the socket to the data handle of the normals
  Framework::DataSocketSink<InwardNormalsData*> socket_normals;

  /// handle to the neighbor cell
  Framework::DataSocketSink<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >    socket_faceNeighCell;

  /// the socket to the data handle of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// BC nodal normals
  std::vector< std::vector<RealVector> > _bcNormals;

  /// physical model (in conservative variables)
  Common::SelfRegistPtr<Physics::LinearizedEuler::LinEuler3DCons> _varSet;

  /// temporary data for holding the nb of variables
  CFuint m_nbEqs;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n}\f$
  RealVector _r1;

  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} \f$
  RealVector _r2;
  
  /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} \f$
  RealVector _r3;

   /// eigen vector corresponding to\f$\vec{u} \cdot \vec{n} + a \f$
  RealVector _r4;
  
    /// temporry for the variables
  RealVector m_var_values;

  /// a vector of string to hold the functions
  std::vector<std::string> m_function_inflow;

  /// a vector of string to hold the variables
  std::vector<std::string> m_vars_inflow;

  /// function for the mean flow
  Framework::VectorialFunction m_function_parser_inflow;
  
  std::vector<RealVector> inflow;
  
  

}; // end of class StrongSubFunctionInletLinEuler3DCons

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSlipWallEuler2DCons_hh
