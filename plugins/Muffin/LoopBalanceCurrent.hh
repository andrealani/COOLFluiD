#ifndef COOLFluiD_Muffin_LoopBalanceCurrent_hh
#define COOLFluiD_Muffin_LoopBalanceCurrent_hh

#include "Muffin/SystemMITReM.hh"
#include "Muffin/Loop.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Looping over system commands by balancing current
class LoopBalanceCurrent : public Loop {

 public:  // functions

  /// Loop condition constructor
  LoopBalanceCurrent(const std::string& name);

  /// Loop condition destructor
  ~LoopBalanceCurrent() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  void setup();


 protected:  // functions

  /// LoopBalanceCurrent termination check (parent virtual implementation)
  bool finish(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs);

  // Get current status as string (parent virtual implementation)
  std::string getStatus(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs);


 private:  // data

  /// BalanceCurrent: total current over all TRS's [A]
  CFreal m_current_bal;

  /// BalanceCurrent: maximum absolute current over all TRS's [A]
  CFreal m_current_max;

  /// BalanceCurrent: total current logarithm L2 maximum level
  CFreal m_current_l2_level;

  /// MITReM system to access boundary currents from, name
  std::string m_sysmitrem_str;

  /// MITReM system to access boundary currents from, pointer
  Common::SafePtr< SystemMITReM > m_sysmitrem;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_LoopBalanceCurrent_hh

