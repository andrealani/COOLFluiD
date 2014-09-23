// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DVarSet_hh
#define COOLFluiD_Physics_LinearAdvSys_LinearAdvSys3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/EquationSetData.hh"
#include "LinearAdvSys/LinearAdvSysVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LinearAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for a 2D Lineariazed Euler physical model.
 *
 * @author Tiago Guintino
 * @author Tomas Kopacek
 */
class LinearAdvSys3DVarSet : public LinearAdvSysVarSet {
public: // classes

  /**
   * Constructor
   * @see LinearAdvSysPhysicalModel
   */
  LinearAdvSys3DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  /**
   * Default destructor
   */
  virtual ~LinearAdvSys3DVarSet();

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup();

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal) = 0;

  /**
   * Set the matrix of the right eigenvectors and the matrix of the eigenvalues
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal) = 0;

  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const = 0; 
  /**
   * Get the normal speed
   */
    CFreal getNormalSpeed(const RealVector& data, const RealVector& normal) const
  {
      throw Common::NotImplementedException (FromHere(),"getNormalSpeed");
      return 0.;
  }

  /**
   * Get some data corresponding to the subset of equations related with
   * this variable set
   * @pre The most concrete ConvectiveVarSet will have to set these data
   */
  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData;
    return eqSetData;
  }
  
  /// Set the IDs corresponding to the velocity components in a State
  virtual void setStateVelocityIDs (std::vector<CFuint>& velIDs)
  { 
    velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2;  velIDs[YY] = 3; 
  }
  
protected:

  /// Compute the convective flux
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals);

  /// Compute the convective flux
  virtual void computeStateFlux(const RealVector& pdata);
 
  /**
   * Get the number of equations of this VarSet
   */
  CFuint getNbEqs() const
  {
    return 4;
  }

  /// Set the vector of the eigenValues
  virtual void computeEigenValues (const RealVector& pdata, const RealVector& normal, RealVector& eValues);

  /// Get the maximum eigenvalue
  virtual CFreal getMaxEigenValue(const RealVector& pdata, const RealVector& normal);

  /// Get the maximum absolute eigenvalue
  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal);

}; // end of class LinearAdvSys2DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdvSys_LinearAdvSys2DVarSet_hh
