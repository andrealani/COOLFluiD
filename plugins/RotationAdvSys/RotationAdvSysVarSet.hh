// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_RotationAdvSys_RotationAdvSysVarSet_hh
#define COOLFluiD_Physics_RotationAdvSys_RotationAdvSysVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "RotationAdvSysTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace RotationAdvSys {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a variable set for an Linearized Euler physical model.
 *
 * @author Nadege Villedieu
 */
class RotationAdvSysVarSet : public Framework::ConvectiveVarSet {
public: // classes
  
  typedef RotationAdvSysTerm PTERM;

 
  /**
   * Constructor
   * @see LinearAdvSysPhysicalModel
   */
  RotationAdvSysVarSet(Common::SafePtr<Framework::BaseTerm> term):
    Framework::ConvectiveVarSet(term),
   _model(term.d_castTo<RotationAdvSysTerm>())
  {
  }
 
  

  /**
   * Default destructor
   */
  virtual ~RotationAdvSysVarSet(){
  
   }

  /**
   * Set up the private data and give the maximum size of states physical
   * data to store
   */
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
  }

  /**
   * Gets the block separator for this variable set
   */
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"RotationAdvSysVarSet::computeJacobians()");
  }

  /**
   * Split the jacobian
   */
  virtual void splitJacobian(RealMatrix& jacobPlus,
                          RealMatrix& jacobMin,
                          RealVector& eValues,
                          const RealVector& normal)
  {
    throw Common::NotImplementedException (FromHere(),"RotationAdvSysVarSet::splitJacobian()");
  }

  /**
   * Get the speed
   */
  virtual CFreal getSpeed(const Framework::State& state) const = 0;

  /**
   * Give dimensional values to the adimensional state variables
   */
  virtual void setDimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "RotationAdvSysVarSet::setDimensionalValues() not implemented");
  }

  /**
   * Give adimensional values to the dimensional state variables
   */
  virtual void setAdimensionalValues(const Framework::State& state,
                                    RealVector& result)
  {
    throw Common::NotImplementedException
      (FromHere(), "RotationAdvSysVarSet::setAdimensionalValues() not implemented");
  }

  /**
   * Get the model
   */
  Common::SafePtr<RotationAdvSysTerm> getModel() const
  {
    return _model;
  }
  
protected:


  /// Compute the convective flux
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals)=0;

  /// Compute the convective flux
  virtual void computeStateFlux(const RealVector& pdata)=0;

  /// acquaintance of the model
  Common::SafePtr<RotationAdvSysTerm> _model;

}; // end of class LinearAdvSysVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearAdvSys

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdvSys_LinearAdvSysVarSet_hh
