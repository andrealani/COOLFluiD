// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_ParalutionLSS_hh
#define COOLFluiD_Numerics_Paralution_ParalutionLSS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/LinearSystemSolver.hh"

#include "Paralution/ParalutionLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    class NumericalCommand;
    class BlockAccumulator;
  }

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a solver of linear systems using the PETSC library
/// @author Andrea Lani
class ParalutionLSS : public Framework::LinearSystemSolver {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  /// @param name missing documentation
  explicit ParalutionLSS(const std::string& name);


  /// Default destructor.
  virtual ~ParalutionLSS();

  /// Configures the method, by allocating its dynamic members.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets a vector with all the NumericalStrategy's this method will use.
  /// @return vector with the strategy pointers.
  virtual std::vector<Common::SafePtr<Framework::NumericalStrategy> > getStrategyList () const;

  /// Prints the Linear System to a file.
  void printToFile(const std::string prefix, const std::string suffix);

  /// Create a block accumulator with chosen internal storage
  /// @return a newly created block accumulator
  /// @post the block has to be deleted outside
  Framework::BlockAccumulator*  createBlockAccumulator
  (const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize, CFreal* ptr) const;
  
  /// Get the LSS system matrix
  Common::SafePtr<Framework::LSSMatrix> getMatrix() const
  {
    return &m_data->getMatrix();
  }

  /// Get the LSS solution vector
  Common::SafePtr<Framework::LSSVector> getSolVector() const
  {
    return &m_data->getSolVector();
  }
  
  /// Gets the LSS right hand side vector
  Common::SafePtr<Framework::LSSVector> getRhsVector() const
  {
    return &m_data->getRhsVector();
  }
  
protected: // abstract interface implementations

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /// Solve the linear system
  /// @see LinearSystemSolver::solveSys()
  void solveSysImpl();

  /// Sets up the data for the method commands to be applied.
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl();

  /// UnSets the data of the method.
  /// @see Method::setMethod()
  virtual void setMethodImpl();

private: // membe data

  ///The Setup command to use
  Common::SelfRegistPtr<ParalutionLSSCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<ParalutionLSSCom> m_unSetup;

  ///The command to use for computing the reference area.
  Common::SelfRegistPtr<ParalutionLSSCom> m_solveSys;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///name of the command that solves the linear system
  std::string m_solveSysStr;

  ///The data to share between ParalutionLSSCom commands
  Common::SharedPtr<ParalutionLSSData> m_data;
  
}; // class ParalutionLSS

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Paralution_ParalutionLSS_hh
