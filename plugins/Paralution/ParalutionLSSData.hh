// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_Paralution_ParalutionLSSData_hh
#define COOLFluiD_Numerics_Paralution_ParalutionLSSData_hh

//////////////////////////////////////////////////////////////////////////////

#include <paralution.hpp>

#include "Framework/LSSData.hh"

#include "Paralution/ParalutionVector.hh"
#include "Paralution/ParalutionMatrix.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

using namespace paralution;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a data object that is accessed by the different
* ParalutionLSSCom 's that compose the ParalutionLSS.
 *
 * @todo there is missing documentation in this class.
 *
 * @author Andrea Lani
 * @author Isaac Alonso
 *
 */
class ParalutionLSSData : public Framework::LSSData {
public:
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Default constructor without arguments.
   */
  ParalutionLSSData(Common::SafePtr<std::valarray<bool> > maskArray,
		    CFuint& nbSysEquations,
		    Common::SafePtr<Framework::Method> owner);

  /**
   * Destructor.
   */
  virtual ~ParalutionLSSData();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();
  
  /**
   * Gets the solution vector
   */
  ParalutionVector& getSolVector()
  {
    return _xVec;
  }

  /**
   * Gets the rhs vector
   */
  ParalutionVector& getRhsVector()
  {
    return _bVec;
  }

  /**
   * Gets the matrix
   */
  ParalutionMatrix& getMatrix()
  {
    return _aMat;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ParalutionLSS";
  }

  /**
   * Gets the relative tolerance
   */
  CFreal getRelativeTol() const
  {
    return _rTol;
  }

  /**
   * Gets the absolute tolerance
   */
  CFreal getAbsoluteTol() const
  {
    return _aTol;
  }

  /**
   * Gets the maximun number of iterations
   */
  CFuint getMaxIter() const
  {
    return _maxIter;
  }

  /**
   * Gets the name of the Krylov Solver
   */
  std::string getKSPName() const
  {
    return _kspTypeStr;
  }

  /** 
   * Gets the Krylov Solver
   */
  IterativeLinearSolver<LocalMatrix<CFreal>, LocalVector<CFreal>, CFreal >& getKSP(){
    return _ls;
  }

  /**
   * Gets the preconditioner
   */
  Preconditioner<LocalMatrix<CFreal>, LocalVector<CFreal>, CFreal >& getPreconditioner(){
    return _p;
  }

  bool getFirstIter(){ return _firstIter; }

  void setFirstIter(bool Iter){ _firstIter = Iter; }

  bool getUseGPU(){return _useGPU;}

  CFuint getVerbose(){return _verboseLevel;}

  CFuint getreBuildRatio(){return _reBuildRatio;}

private:

  /// Linear Solver
  GMRES<LocalMatrix<CFreal>, LocalVector<CFreal>, CFreal > _ls;

  /// Preconditioner
  Jacobi<LocalMatrix<CFreal>,LocalVector<CFreal>,CFreal> _p;

  /// solution vector
  ParalutionVector _xVec;

  /// rhs vector
  ParalutionVector _bVec;
  
  /// matrix
  ParalutionMatrix _aMat;
  
  /// Number of Krylov spaces
  CFuint _nbKsp;
  
  /// relative tolerance for iterative solver
  CFreal _rTol;

  /// absolute tolerance for iterative solver
  CFreal _aTol;

  /// Maximun number of iterations
  CFuint _maxIter;

  /// VerboseLevel for the Liner System Solver
  CFuint _verboseLevel;

  /// Number of iterations before reBuilding the preconditioner
  CFuint _reBuildRatio;
  
  /// First iteration
  bool _firstIter;

  /// Flag to solve the system on the accelerator
  CFuint _useGPU;
 
  /// Krylov solver name
  std::string _kspTypeStr;
      
}; // end of class ParalutionLSSData

//////////////////////////////////////////////////////////////////////////////

/// Definition of a command for Paralution
typedef Framework::MethodCommand<ParalutionLSSData> ParalutionLSSCom;

/// Definition of a command provider for Paralution
typedef Framework::MethodCommand<ParalutionLSSData>::PROVIDER ParalutionLSSComProvider;

//////////////////////////////////////////////////////////////////////////////


    } // namespace Paralution

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_Paralution_ParalutionLSSData_hh
