#ifndef COOLFluiD_Numerics_MatrixStabilityMethodData_hh
#define COOLFluiD_Numerics_MatrixStabilityMethodData_hh

//////////////////////////////////////////////////////////////////////////////

#include <fstream>

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Data Object that is accessed by the different
 * MatrixStabilityMethodWriterCom 's that compose the MatrixStabilityMethodWriter.
 *
 * @see RungeKuttaLSCom
 *
 * @author Kris Van den Abeele
 */
class MatrixStabilityMethodData : public Framework::ConvergenceMethodData {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments
   */
  MatrixStabilityMethodData(Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor
   */
  ~MatrixStabilityMethodData();

  /**
   * Configure the data from the supplied arguments.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "MatrixStabilityMethodData";
  }

  /**
   * @return m_nbrStates
   */
  CFuint getNbrStates()
  {
    return m_nbrStates;
  }

  /**
   * set m_nbrStates
   */
  void setNbrStates(const CFuint nbrStates)
  {
    m_nbrStates = nbrStates;
  }

  /**
   * @return m_setAllStatesToZero
   */
  bool getSetAllStatesToZero()
  {
    return m_setAllStatesToZero;
  }

  /**
   * set m_setAllStatesToZero
   */
  void setAllStatesToZero(const bool setStatesToZero)
  {
    m_setAllStatesToZero = setStatesToZero;
  }

  /**
   * @return m_setStateToZero
   */
  bool getSetStateToZero()
  {
    return m_setStateToZero;
  }

  /**
   * set m_setStateToZero
   */
  void setStateToZero(const bool setStateToZero)
  {
    m_setStateToZero = setStateToZero;
  }

  /**
   * @return m_stateIdx
   */
  CFuint getStateIdx()
  {
    return m_stateIdx;
  }

  /**
   * set m_stateIdx
   */
  void setStateIdx(const CFuint stateIdx)
  {
    m_stateIdx = stateIdx;
  }

  /**
   * @return m_outputFileName
   */
  std::string getOutputFileName()
  {
    return m_outputFileName;
  }

  /**
   * @return m_outputFile
   */
  Common::SafePtr< std::ofstream > getOutputFile()
  {
    return &m_outputFile;
  }

  /**
   * @return m_writeBinary
   */
  bool writeBinary()
  {
    return m_writeBinary;
  }

private: // data

  /// number of states in the mesh
  CFuint m_nbrStates;

  /// boolean telling whether to set all states to zero or to manipulate one state
  bool m_setAllStatesToZero;

  /// boolean telling whether to set one state to zero or to one
  bool m_setStateToZero;

  /// index of the state to manipulate
  CFuint m_stateIdx;

  /// name for output file
  std::string m_outputFileName;

  /// output file
  std::ofstream m_outputFile;

  /// boolean telling whether to write a binary file
  bool m_writeBinary;

}; // end of class MatrixStabilityMethodData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for MatrixStabilityMethod
typedef Framework::MethodCommand<MatrixStabilityMethodData> MatrixStabilityMethodCom;

/// Definition of a command provider for MatrixStabilityMethod
typedef Framework::MethodCommand<MatrixStabilityMethodData>::PROVIDER MatrixStabilityMethodComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MatrixStabilityMethodData_hh
