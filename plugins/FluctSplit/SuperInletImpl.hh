#ifndef COOLFluiD_Numerics_FluctSplit_SuperInletImpl_hh
#define COOLFluiD_Numerics_FluctSplit_SuperInletImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/SuperInlet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a supersonic inlet command
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API SuperInletImpl : public SuperInlet {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  SuperInletImpl(const std::string& name);

  /// Default destructor
  ~SuperInletImpl();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Set up the member data
  virtual void setup();

protected:

  /// Execute on the current TRS
  void executeOnTrs();

protected: // data

  /// the socket to the data handle of the boundary state neighbors
  Framework::DataSocketSink<std::valarray<Framework::State*> >
  socket_bStatesNeighbors;

  /// the socket to the data handle of the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the socket to the data handle of arrays of flags specifying if the
  /// time jacobian contribution of certain variables in boundary states
  /// have to be discarded
  Framework::DataSocketSink<std::vector<bool> > socket_discardTimeJacob;

  /// factor that multiplies the coefficient for the diagonal entry
  /// in the linear system matrix
  CFreal _diagCoeffFactor;

}; // end of class SuperInletImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SuperInletImpl_hh
