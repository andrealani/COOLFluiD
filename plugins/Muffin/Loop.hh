#ifndef COOLFluiD_Muffin_Loop_hh
#define COOLFluiD_Muffin_Loop_hh

#include "Common/COOLFluiD.hh"
#include "Common/TaggedObject.hh"
#include "Muffin/Muffin.hh"
#include "Muffin/MuffinData.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Base class for looping over system commands
class Loop : public MuffinCom,
             public Common::TaggedObject {

 public:  // functions

  /// Loop condition constructor
  Loop(const std::string& name);

  /// Loop condition destructor
  virtual ~Loop() {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  virtual void setup();

  /// Execute the loop
  void execute();

  /// Set this loop as master
  void setMaster() { m_master = true; }

  /// If this loop is mastering simulation
  bool isMaster() const { return m_master; }


 protected:  // functions

  /// Log information, debug and error messages
  void log(const std::string& msg) { getMethodData().log("Loop " + getName() + ": " + msg); }
  void ver(const std::string& msg) { getMethodData().ver("Loop " + getName() + ": " + msg); }
  void err(const std::string& msg) { getMethodData().err("Loop " + getName() + ": " + msg); }

  /// Loop termination check
  virtual bool finish(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs);

  // Get current status as string
  virtual std::string getStatus(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs);

  /// Write convergence file
  void writeConvergence(CFuint i, CFreal t, const std::vector< CFreal >& logl2_states, const std::vector< CFreal >& logl2_rhs);


 protected:  // data

  /// Command pointers to loop on
  std::vector< Common::SafePtr< MuffinCom > > m_commands;


 private:  // data

  /// Convergence file to write to
  std::string m_cvg_file;

  /// Iterate: number of iterations
  CFuint m_niterations;

  /// ConvergeResidual/Solution: logarithm L2 maximum level
  CFreal m_clevel;

  /// ConvergeSolution: variable to be tested (name(s))
  std::vector< std::string > m_cv;

  /// ConvergeSolution: variable to be tested (index/indices)
  std::vector< CFuint > m_cvi;

  /// ConvergeResidual: equation to be tested (index/indices)
  std::vector< CFuint > m_cei;

  /// TimeStep: loop time duration
  CFreal m_tduration;

  /// TimeStep: loop time step
  CFreal m_tstep;

  /// TimeStep: loop starting time, set by first loop
  CFreal m_tstart;

  /// Command names to loop on, either Systems or other Loops
  std::vector< std::string > m_command_names;

  /// If this loop is mastering the simulation
  bool m_master;

};


/// Loop iteratively
class LoopIterate : public Loop {
 public:
  LoopIterate(const std::string& name) : Loop(name) { attachTag("LoopIterate"); }
};


/// Loop by time-stepping
class LoopTimeStep : public Loop {
 public:
  LoopTimeStep(const std::string& name) : Loop(name) { attachTag("LoopTimeStep"); }
};


/// Loop by converging linear system right-hand side vector
class LoopConvergeResidual : public Loop {
 public:
  LoopConvergeResidual(const std::string& name) : Loop(name) { attachTag("LoopConvergeResidual"); }
};


/// Loop by converging solution
class LoopConvergeSolution : public Loop {
 public:
  LoopConvergeSolution(const std::string& name) : Loop(name) { attachTag("LoopConvergeSolution"); }
};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_Loop_hh

