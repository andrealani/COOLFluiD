#ifndef COOLFluiD_IO_TecplotWriter_WriteInstantAndAvgSolutionHighOrder_hh
#define COOLFluiD_IO_TecplotWriter_WriteInstantAndAvgSolutionHighOrder_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FileWriter.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/ProxyDofIterator.hh"
#include "Framework/VarSetTransformer.hh"

#include "TecplotWriter/WriteSolution.hh"
#include "TecplotWriter/TecWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class Node; }

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action
/// to write the averaged and the instantaneous MeshData solution to a Tecplot format file
/// for visualization. Each high-order cell is written to a different zone,
/// (is in fact a miniature unstructured mesh)
/// This writer is suited for methods like discontinuous Galerkin, spectral volume,
/// spectral difference and related methods (the shape functions should be implemented though!!!)
///
/// @author Matteo Parsani

class WriteInstantAndAvgSolutionHighOrder : public WriteSolution {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit WriteInstantAndAvgSolutionHighOrder(const std::string& name);

  /**
   * Destructor.
   */
  ~WriteInstantAndAvgSolutionHighOrder()
  {
  }

  /**
    * Set up private data
   */
  void setup();

  /**
   * UnSet up private data and data of the aggregated classes
   * in this command after processing phase
   */
  void unsetup();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

protected:

  /// socket for PastStates's
  Framework::DataSocketSource< Framework::State*> socket_pastAvgStates;



}; // class WriteInstantAndAvgSolutionHighOrder

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_IO_TecplotWriter_WriteSolutionHighOrder_hh

