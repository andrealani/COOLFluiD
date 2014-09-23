#ifndef COOLFluiD_Numerics_FluctSplit_ArtViscCoeffUnSetup_hh
#define COOLFluiD_Numerics_FluctSplit_ArtViscCoeffUnSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/Storage.hh"
#include "Framework/ProxyDofIterator.hh"
#include "FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// executed in order to deallocate the socket holding
/// the values of the artificial viscosity coefficient
/// used for the carbuncle fix.
/// To be used in conjunction with ArtViscCoeffSetup
/// @author Nadege Villedieu
/// @author Tiago Quintino
/// @author JGM




class FluctSplit_API ArtViscCoeffUnSetup : public FluctuationSplitCom {
public: // member functions

  /// Constructor.
  explicit ArtViscCoeffUnSetup(const std::string& name);

  /// Destructor.
  virtual ~ArtViscCoeffUnSetup();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /// Executes this command
  virtual void execute();

protected: // data

  /// socket for theta blending coefficient
  Framework::DataSocketSink<CFreal> socket_ArtViscCoeff;

}; // class UnSetup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StdUnSetup_hh

