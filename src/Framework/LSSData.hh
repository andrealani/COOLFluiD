// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_LSSData_hh
#define COOLFluiD_Framework_LSSData_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/LinearSystemSolver.hh"
#include "Framework/MethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a data object that is accessed by the different
/// LSSCom 's that compose the LSS.
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API LSSData : public Framework::MethodData {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  LSSData(Common::SafePtr<std::valarray<bool> > maskArray,
	  CFuint& nbSysEquations,
	  Common::SafePtr<Method> owner);
  
  /// Destructor.
  virtual ~LSSData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// Get the LSS idx mapping (from local to global LSS ids)
  /// Accessor that does not allow mutation for clients.
  Framework::LSSIdxMapping& getLocalToGlobalMapping()
  {
    return m_localToGlobal;
  }
  
  /// Get the LSS idx mapping (from local to locally updatable LSS ids)
  /// -1 is the default value (meaning not updatable)
  /// Accessor that does not allow mutation for clients.
  Framework::LSSIdxMapping& getLocalToLocallyUpdatableMapping()
  {
    return m_localToLocallyUpdateble;
  }
  
  /// Gets the size of the system of equations to solve
  const CFuint& getNbSysEquations() const
  {
    return m_nbSysEquations;
  }

  /// Mask array that specifies which equations are solved by the current LSS
  Common::SafePtr<std::valarray<bool> > getMaskArray() const
  {
    return m_maskArray;
  }

  /// Gets the Maximum number of iterations
  CFuint getMaxIterations() const
  {
    return m_maxIter;
  }

  /// Gets the save rate to write matrix and vector to file
  CFuint getSaveRate() const
  {
    return m_saveRate;
  }

  /// Gets the rate at which preconditioner must be recomputed
  CFuint getPreconditionerRate() const
  {
    return m_preconditionerRate;  
  }

  /// Returns if the convergence history of the solver should be outputed
  bool isOutput() const
  {
    return m_isOutput;
  }

  /// Returns the linear system should be saved to file when solving
  bool isSaveSystemToFile() const
  {
    return m_saveSystemToFile;
  }

  /**
   * Flag telling to use block preconditioner matrix
   */
  bool useBlockPreconditionerMatrix() const
  {
    return _useBlockPreconditionerMatrix;
  }
  
  /// Flag telling to use GPU support (if available)
  bool useGPU() const {return m_useGPU;}
  
  /// Flag telling to use node-based sparsity an assembly (instead of state-based)
  bool useNodeBased() const {return m_useNodeBased;}
  
 private: // data
  
  /// mapping local to global indices numbering
  Framework::LSSIdxMapping m_localToGlobal;

  /// mapping local to localy updatable indices numbering
  Framework::LSSIdxMapping m_localToLocallyUpdateble;

  /// mask array to specify which equations are solved by the current LSS
  Common::SafePtr<std::valarray<bool> > m_maskArray;

  /// the size of the system of equations to solve
  /// this has to be a reference because at construction time
  /// the m_nbSysEquations == 0 and will be changed at configuration
  /// time by LinearSystemSolver
  CFuint& m_nbSysEquations;

  /// Flag telling to use block preconditioner matrix
  bool _useBlockPreconditionerMatrix;

  /// maximum number of iterations for iterative solver
  CFuint m_maxIter;
  
  /// save rate to write matrix and vector to file
  CFuint m_saveRate;

  /// rate at which preconditioner must be recomputed 
  CFuint m_preconditionerRate;
 
  /// write output
  bool m_isOutput;
  
  /// save system to file
  bool m_saveSystemToFile;
  
  /// use GPU support (if available)
  bool m_useGPU;
  
  /// use node-based sparsity and assembly (instead of state-based)
  bool m_useNodeBased;
  
}; // end of class LSSData

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_LSSData_hh
