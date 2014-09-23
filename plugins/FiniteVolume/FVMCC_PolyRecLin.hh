#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRecLin_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRecLin_hh

//////////////////////////////////////////////////////////////////////////////

#include "PolyReconstructorLin.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a constant polynomial reconstructor for FVM
 *
 * @author Andrea Lani
 *
 */
class FVMCC_PolyRecLin : public PolyReconstructorLin {
public:

  /**
   * Constructor
   */
  FVMCC_PolyRecLin(const std::string& name);

  /**
   * Default destructor
   */
  ~FVMCC_PolyRecLin();

  /**
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /**
   * Update the weights when nodes are moving
   */
  virtual void updateWeights();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::PolyReconstructor<CellCenterFVMData>::configure(args);
  }

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = Framework::PolyReconstructor<CellCenterFVMData>::needsSockets();

    result.push_back(&socket_nodes);
    result.push_back(&socket_states);
    result.push_back(&socket_gstates);

    return result;
  }

  /**
   * Compute the gradients
   */
  virtual void computeGradients() = 0;

  /**
   * Compute the limiters for all the cells
   */
  virtual void computeLimiters();

private:

  /// pointer to the TRS of cells
  Common::SafePtr<Framework::TopologicalRegionSet> _cells;

  /// cell builder for cell centered FVM schemes
  Framework::GeometricEntityPool<Framework::CellTrsGeoBuilder> _cellBuilder;

protected:

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;


}; // end of class FVMCC_PolyRecLin

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRecLin_hh
