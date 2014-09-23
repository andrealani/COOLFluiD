#ifndef COOLFluiD_Numerics_AeroCoef_Extract2DSectionCC_hh
#define COOLFluiD_Numerics_AeroCoef_Extract2DSectionCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FiniteVolume {
      class CellCenterFVMData;
      class DerivativeComputer;
    }

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class extracts a 2D section starting from a point on a TRS until
 * the end of the mesh along a sectionDirection for
 * @see CellCenterFVM
 *
 * @author Thomas Wuilbaut
 *
 */

class Extract2DSectionCC : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Extract2DSectionCC(const std::string& name);

  /**
   * Default destructor
   */
  ~Extract2DSectionCC();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the required values
   */
  void computeValues();

  /**
   * Extract Boundary Layer Profile
   */
  void extractSection(Framework::GeometricEntity* currFace);

private:

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the ghost state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  // the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;

   /// storage of the face normals
  Framework::DataSocketSink<CFreal> socket_normals;

  // pointer to the data of the cell centered FVM method
  Common::SafePtr<Numerics::FiniteVolume::CellCenterFVMData> _fvmccData;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> _updateVar;

  // temporary unit normal
  std::vector<CFreal> _sectionDirection;
  RealVector _sectionNormal;

  // Storage for choosing when to save the wall values file
  CFuint _saveRate;

  // append the iteration number in the output file
  bool _appendIter;

  // append the time in the output file
  bool _appendTime;

  // name of the output file
  std::string _outputFile;

  //Flag to extract along normal to TRS
  bool _extractAlongNormal;

  //coordinates for the extraction of section
  std::vector<CFreal> _extractCoord;
  RealVector _initExtract;

  CFreal _tolerance;

}; // end of class Extract2DSectionCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_Extract2DSectionCC_hh
