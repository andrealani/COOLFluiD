// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Trilinos_TrilinosLSS_hh
#define COOLFluiD_Numerics_Trilinos_TrilinosLSS_hh

//////////////////////////////////////////////////////////////////////////////



#include "Trilinos.hh"
#include "Framework/LinearSystemSolver.hh"
#include "Framework/BlockAccumulator.hh"
#include "TrilinosLSSData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class NumericalCommand;
    class LSSMatrix;
    class LSSVector;
  }

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a solver of linear systems using the Tilinos library
 *
 * @author Andrea Lani
 *
 */
class TrilinosLSS : public Framework::LinearSystemSolver {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   *
   * @param name of this method used for configuring it
   */
  explicit TrilinosLSS(const std::string& name);


  /**
   * Default destructor.
   */
  virtual ~TrilinosLSS();

  /**
   * Configures the method, by allocating its dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Prints the Linear System to a file.
   */
  void printToFile(const std::string prefix, const std::string suffix);

  /**
   * Create a block accumulator with chosen internal storage
   * @return a newly created block accumulator
   * @post the block has to be deleted outside
   */
  Framework::BlockAccumulator* createBlockAccumulator(const CFuint nbRows,
                                             const CFuint nbCols,
                                             const CFuint subBlockSize,
                                             CFreal* ptr) const;

  /**
   * Get the LSS system matrix
   */
  Common::SafePtr<Framework::LSSMatrix> getMatrix() const
  {
    return Common::SafePtr<Framework::LSSMatrix>((Framework::LSSMatrix*)m_data->getMatrix());
  }

  /**
   * Get the LSS solution vector
   */
  Common::SafePtr<Framework::LSSVector> getSolVector() const
  {
    return Common::SafePtr<Framework::LSSVector>((Framework::LSSVector*)m_data->getSolVector());
  }

  /**
   * Get the LSS right hand side vector
   */
  Common::SafePtr<Framework::LSSVector> getRhsVector() const
  {
    return Common::SafePtr<Framework::LSSVector>((Framework::LSSVector*)m_data->getRhsVector());
  }

  /**
   * Get the Trilinos system matrix
   */
  TrilinosMatrix* getTrilinosMatrix()
  {
    return m_data->getMatrix();
  }

  /**
   * Get the Trilinos solution vector
   */
  TrilinosVector* getTrilinosSolVector()
  {
    return m_data->getSolVector();
  }

  /**
   * Get the Trilinos right hand side vector
   */
  TrilinosVector* getTrilinosRhsVector()
  {
    return m_data->getRhsVector();
  }

protected: // abstract interface implementations

  /**
   * Gets the Data aggregator of this method
   * @return SafePtr to the MethodData
   */
  virtual Common::SafePtr< Framework::MethodData > getMethodData () const;

  /**
   * Solve the linear system
   * @see LinearSystemSolver::solveSys()
   */
  void solveSysImpl();

  /**
   * Sets up the data for the method commands to be applied.
   * @see Method::unsetMethod()
   */
  virtual void unsetMethodImpl();

  /**
   * UnSets the data of the method.
   * @see Method::setMethod()
   */
  virtual void setMethodImpl();

private: // member data

  ///The Setup command to use
  Common::SelfRegistPtr<TrilinosLSSCom> m_setup;

  ///The UnSetup command to use
  Common::SelfRegistPtr<TrilinosLSSCom> m_unSetup;

  ///The command to use for computing the reference area.
  Common::SelfRegistPtr<TrilinosLSSCom> m_solveSys;

  ///The Setup string for configuration
  std::string m_setupStr;

  ///The UnSetup string for configuration
  std::string m_unSetupStr;

  ///name of the command that solves the linear system
  std::string m_solveSysStr;

  ///The data to share between TrilinosLSSCom commands
  Common::SharedPtr<TrilinosLSSData> m_data;

}; // class TrilinosLSS

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Trilinos_TrilinosLSS_hh
