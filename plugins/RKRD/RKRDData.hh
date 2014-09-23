#ifndef COOLFluiD_Numerics_RKRD_RKRDData_hh
#define COOLFluiD_Numerics_RKRD_RKRDData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Config/ConfigObject.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RKRD {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a Data Object that is accessed by the different
/// RKRDCom 's that compose the RKRD.
/// @see RKRDCom
/// @author Thomas Wuilbaut
class RKRDData : public Framework::ConvergenceMethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  RKRDData(Common::SafePtr<Framework::Method> owner);

  /// Destructor
  ~RKRDData();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName() { return "RKRDData"; }

  /// Access the order of the scheme
  CFuint getOrder() const { return m_order; }
  /// Access the order of the scheme
  CFuint& K() { return m_curr_k; }

  RealMatrix getAlpha() const { return m_alpha; }

  RealMatrix getBeta() const { return m_beta; }

private:

  /// order of the Runge-Kutta scheme
  CFuint m_order;

  /// vector of the coeficients of the RKRD method
  RealMatrix m_alpha;

  /// vector of the coeficients of the RKRD method
  RealMatrix m_beta;

  /// current kstep
  CFuint m_curr_k;
}; // end of class RKRDData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for RKRD
typedef Framework::MethodCommand<RKRDData> RKRDCom;

/// Definition of a command provider for RKRD
typedef Framework::MethodCommand<RKRDData>::PROVIDER RKRDComProvider;

//////////////////////////////////////////////////////////////////////////////

    } // namespace RKRD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RKRD_RKRDData_hh
