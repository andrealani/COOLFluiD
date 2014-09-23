#ifndef COOLFluiD_Numerics_FluctSplit_PeriodicBCImpl_hh
#define COOLFluiD_Numerics_FluctSplit_PeriodicBCImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/PeriodicBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a periodic boundary condition for implicit convergence
/// @author Tiago Quintino
class FluctSplit_API PeriodicBCImpl : public PeriodicBC {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  PeriodicBCImpl(const std::string& name);

  /// Default destructor
  virtual ~PeriodicBCImpl();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Configures this object with supplied arguments.
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /// Execute on a set of dofs
  virtual void executeOnTrs();

protected: // data

  /// the socket to the data handle of the boundary state neighbors
  Framework::DataSocketSink<std::valarray<Framework::State*> >  socket_bStatesNeighbors;

  /// storage for a block of values got from the jacobian matrix
  RealMatrix _jacobElem;

  /// storage for a block of values got from the jacobian matrix
  RealMatrix _jacobElema;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _irc;

  /// indexes of the row of each value in the block
  std::valarray<CFint> _ira;

  /// indexes of the columns of each value in the block
  std::valarray<CFint> _in;

  /// storage for a block of values to insert in the jacobian matrix
  RealMatrix _block;

}; // end of class PeriodicBCImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_PeriodicBCImpl_hh
