// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_SimulationStatus_hh
#define COOLFluiD_Framework_SimulationStatus_hh

//////////////////////////////////////////////////////////////////////////////
#include <boost/filesystem/path.hpp>

#include "Common/NonCopyable.hh"
#include "Common/CFMap.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a singleton object where global infos (number
/// of iterations and residual) about the status of the Simulation are held
/// @author Thomas Wuilbaut
/// @author Tiago Quintino
class Framework_API SimulationStatus : public Common::NonCopyable<SimulationStatus> {
public:

  /// @return the instance of this singleton
  static SimulationStatus& getInstance();

  /// Resets the SimulationStatus
  void resetAll()
  {
    resetResidual();
    resetTime();
    setNbIter(0);
  }

  /// Get a vector with the residual in each SubSystem
  std::vector<CFreal>& getResiduals()
  {
    return m_residual;
  }

  /// Get residual for regression testing
  CFreal& getLastResidual() {return m_lastResidual;}
  
//////////////////////////////////////////////////////////////////////////////
//
//   Coupling Residual related functions
//
//////////////////////////////////////////////////////////////////////////////

  /// Sets the residual names
  void addCouplingResidualNames(const std::string residualName);

  /// Sets the couplingResidual corresponding to residualName
  void setCouplingResidual(const CFreal residual, const std::string residualName);

  /// Get the value of the coupling residual corresponding to the name
  /// @post if m_couplingResidual vector has no size then
  ///       return 0
  /// @return the given residual
  CFreal getCouplingResidual(const std::string residualName);

  /// Get the value of all the coupling residuals
  /// @return the value of all the computed coupling residuals
  std::vector<CFreal>& getAllCouplingResiduals()
  {
    return m_couplingResiduals;
  }

  /// Reset to zero the coupling residuals
  void resetCouplingResidual()
  {
    for(CFuint i=0; i<m_couplingResiduals.size(); ++i)
    {
      m_couplingResiduals[i] = 0.0;
    }
  }



  /// Get the number of iterations
  CFuint getNbIter() const
  {
    return m_iter;
  }

  /// Reset to zero the residual
  void resetResidual();

  /// Reset to zero the time
  void resetTime();

  /// Reset the number of iterations
  void setNbIter(const CFuint p)
  {
    m_iter = p;
  }

  /// Increment the number of iterations
  void incrementNbIter()
  {
    ++m_iter;
  }

  /// Get the number of SubSystem's
  CFreal getNbSubSystems() const
  {
    return m_subSystems.size();
  }

  /// Get the SubSystem names
  std::vector<std::string>& getSubSystems()
  {
    return m_subSystems;
  }

  /// Get the number of SubSystem's
  void setSubSystems(const std::vector<std::string>& subSys);

  /// Set the info to append to outputfiles
  void setAppendIter(const bool iter)
  {
    m_appendIter = iter;
  }

  /// Needs to append the iteration number??
  bool isAppendIter()
  {
    return m_appendIter;
  }

  /// Set the flag: should the subsystem restart from the previous iteration
  void setRestart(const bool restart)
  {
    m_restart = restart;
  }

  /// Set the flag: should the subsystem restart from the previous iteration
  bool isRestart()
  {
    return m_restart;
  }

  /// Get the last output file for a given subsystem
  /// @param iSub index of the subsystem
  boost::filesystem::path getLastOutputFile(const std::string& subSystemName,const std::string& nspName)
  {
    std::string keyName = subSystemName + "_" + nspName;
    return m_lastOutputFiles.find(keyName)->second;
  }

  /// Set the last output file for a given subsystem
  /// @param outputFile name of the file
  /// @param nsp name of the namespace in the subSystem
  /// @param subSystemName name of the subSystem
  void setLastOutputFile(const boost::filesystem::path& outfile,const std::string& nspName, const std::string& subSystemName)
  {
    std::string keyName = subSystemName + "_" + nspName;
//CFout << "Setting m_lastOutputFiles[" <<keyName<< "] = " <<outfile.string() <<"\n";
    CFuint nbKeyFound = m_lastOutputFileConst.count(keyName);
    bool fixedValue = false;
    if(nbKeyFound > 0) fixedValue = m_lastOutputFileConst.find(keyName)->second;
    if(!fixedValue)
    {
      m_lastOutputFiles[keyName] = outfile;
    }
  }

  /// Set the last output file for a given subsystem
  /// @param outputFile name of the file
  /// @param nsp name of the namespace in the subSystem
  /// @param subSystemName name of the subSystem
  void setLastOutputFile(const boost::filesystem::path& outfile,const std::string& nspName, const std::string& subSystemName, const bool fixedValue)
  {
    std::string keyName = subSystemName + "_" + nspName;
//CFout << "Setting m_lastOutputFiles[" <<keyName<< "] = " <<outfile.string() <<"\n";

    m_lastOutputFiles[keyName] = outfile;
    m_lastOutputFileConst[keyName] = fixedValue;

  }

  /// Get the time of a subsystem
  CFreal getSimulationTime(const std::string subSysName);

private:

  /// Constructor
  SimulationStatus();

  /// Default destructor
  ~SimulationStatus();

private:

  /// the number of iterations
  CFuint m_iter;

  /// the current time (unsteady computations)
  CFreal m_currentTime;

  ///flag: should the iteration number be append to the filename
  bool m_appendIter;

  ///flag: should the subsystem restart from the previous iteration
  bool m_restart;

  /// the residual for the regression testing
  CFreal m_lastResidual;
  
  /// the Simulation residual number stored
  std::vector<CFreal> m_residual;

  /// the residual for the coupling between subsystems stored
  std::vector<CFreal> m_couplingResiduals;

  ///Map with names and ID of the coupling residuals
  Common::CFMap<std::string,CFuint>* m_couplingResidualNames;

  /// the Simulation time number stored
  std::vector<CFreal> m_time;

  /// the vector with the names of the SubSystems
  std::vector<std::string> m_subSystems;

  /// the vector with the names of the last output file of each namespace of each SubSystem
  std::map<std::string, boost::filesystem::path, std::less<std::string> > m_lastOutputFiles;

  std::map<std::string, bool, std::less<std::string> > m_lastOutputFileConst;
}; // end of class SimulationStatus

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_SimulationStatus_hh
