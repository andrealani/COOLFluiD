#ifndef COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega3D_hh
#define COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler3D
/// @author Andrea Lani
class FluctSplit_API StrongNoSlipWallAdiabaticNSKOmega3D : public FluctuationSplitCom {

public:

  /// Constructor.
  StrongNoSlipWallAdiabaticNSKOmega3D(const std::string& name);

  /// Default destructor
  ~StrongNoSlipWallAdiabaticNSKOmega3D();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

/// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
     static void defineConfigOptions(Config::OptionList& options);


protected:

  /// Execute on a set of dofs
  void executeOnTrs();

private:


  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // the socket to the data handle of the upadteCoeff
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  // the socket to the data handle of the isUpdated flags
  Framework::DataSocketSink<bool> socket_isUpdated;


  /// flag telling if initialization is performed
  bool _useForInitialization;

 /// the socket to the data handle of the rhs
   Framework::DataSocketSink<CFreal> socket_rhs;
 
    /// x-coponent of the wall velocity
    CFreal m_u0;
    /// y-coponent of the wall velocity
   CFreal m_v0;
 
   ///z-coponent of the wall velocity
   CFreal m_w0;
                 
}; // end of class StrongNoSlipWallAdiabaticNSKOmega3DImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongNoSlipWallAdiabaticNSKOmega3DImpl_hh
