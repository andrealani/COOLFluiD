#ifndef COOLFluiD_Numerics_MatrixStabilityMethodWriter_AddMatrixColumnToFile_hh
#define COOLFluiD_Numerics_MatrixStabilityMethodWriter_AddMatrixColumnToFile_hh

//////////////////////////////////////////////////////////////////////////////

#include "MatrixStabilityMethodWriter/MatrixStabilityMethodData.hh"
#include "MathTools/RealVector.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MatrixStabilityMethodWriter {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This commands adds a matrix column to the file
   * @author Kris Van den Abeele
   */
class AddMatrixColumnToFile : public MatrixStabilityMethodCom {
public:

  /**
   * Constructor.
   */
  explicit AddMatrixColumnToFile(const std::string& name);

  /**
   * Destructor.
   */
  ~AddMatrixColumnToFile()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  /// handle to states
  Framework::DataSocketSink<CFreal> socket_rhs;

}; // class AddMatrixColumnToFile

//////////////////////////////////////////////////////////////////////////////

    } // namespace MatrixStabilityMethodWriter

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MatrixStabilityMethodWriter_AddMatrixColumnToFile_hh
