// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Pardiso_Pardiso_hh
#define COOLFluiD_Pardiso_Pardiso_hh

#include "Pardiso/PardisoData.hh"
#include "Framework/LinearSystemSolver.hh"

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
    class BlockAccumulator;
  }

  namespace Pardiso {


/// This class interfaces to Pardiso
class Pardiso : public Framework::LinearSystemSolver {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the options
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  explicit Pardiso(const std::string& name);

  /// Destructor
  ~Pardiso();

  /// Sets up the data for the method commands to be applied
  virtual void setMethodImpl();

  /// UnSets the data of the method
  virtual void unsetMethodImpl();

  /// Configures the method, by allocating its dynamic members
  virtual void configure ( Config::ConfigArgs& args );

  /// Solve the linear system
  void solveSysImpl();

  /// Prints the Linear System to a file
  void printToFile(const std::string prefix, const std::string suffix);

  /**
   * Create a block accumulator with chosen internal storage
   * @return a newly created block accumulator
   * @post the block has to be deleted outside
   */
  Framework::BlockAccumulator* createBlockAccumulator(
    const CFuint nbRows, const CFuint nbCols, const CFuint subBlockSize )
    const;

  /// Get the LSS system matrix
  Common::SafePtr< Framework::LSSMatrix > getMatrix() const {
    return &m_data->getMatrix();
  }

  /// Get the LSS solution vector
  Common::SafePtr< Framework::LSSVector > getSolVector() const {
    return &m_data->getSolVector();
  }

  /// Get the LSS right hand side vector
  Common::SafePtr< Framework::LSSVector > getRhsVector() const {
    return &m_data->getRhsVector();
  }

  /// Get the PardisoMatrix
  PardisoMatrix& getMatrix() {
    return m_data->getMatrix();
  }

  /// Get the PardisoVector for the solution
  PardisoVector& getSolVector() {
    return m_data->getSolVector();
  }

  /// Get the PardisoVector for the RHS
  PardisoVector& getRhsVector() {
    return m_data->getRhsVector();
  }


protected:

  /**
   * Get the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;


private:

  /// The Setup command to use
  Common::SelfRegistPtr< PardisoCom > m_setup;

  /// The UnSetup command to use
  Common::SelfRegistPtr< PardisoCom > m_unSetup;

  /// The command to use for computing the reference area
  Common::SelfRegistPtr< PardisoCom > m_solveSys;

  /// The Setup string for configuration
  std::string m_setupStr;

  /// The UnSetup string for configuration
  std::string m_unSetupStr;

  /// Name of the command that solves the linear system
  std::string m_solveSysStr;

  /// Data to share between PardisoCom commands
  Common::SharedPtr< PardisoData > m_data;

}; // class Pardiso


  } // namespace Pardiso
} // namespace COOLFluiD

#endif // COOLFluiD_Pardiso_Pardiso_hh

