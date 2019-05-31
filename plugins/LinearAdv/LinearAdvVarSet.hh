// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.
#ifndef COOLFluiD_Physics_LinearAdv_LinearAdvVarSet_hh
#define COOLFluiD_Physics_LinearAdv_LinearAdvVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/EquationSetData.hh"
#include "LinearAdvTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {


  namespace Physics {
namespace LinearAdv {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a LinearAdv physical model 2D for conservative
/// variables
/// @author Andrea Lani
class LinearAdvVarSet : public Framework::ConvectiveVarSet {
  
public: // classes
  typedef LinearAdvVarSet PVARSET;
  typedef LinearAdvTerm PTERM;
    
  /**
   * Constructor
   * @see LinearAdvPhysicalModel
   */
  LinearAdvVarSet(Common::SafePtr<Framework::BaseTerm> term) :
    Framework::ConvectiveVarSet(term),
    _model(term.d_castTo<LinearAdvTerm>())
  {
  }

  /// Default destructor
  virtual ~LinearAdvVarSet()
  {
  }

  /// Set up the private data and give the maximum size of states physical
  /// data to store
  virtual void setup()
  {
    Framework::ConvectiveVarSet::setup();
  }

  /// Gets the block separator for this variable set
  virtual CFuint getBlockSeparator() const = 0;

  /**
   * Set the jacobians
   */
  virtual void computeJacobians()
  {
    throw Common::NotImplementedException (FromHere(),"LinearAdvVarSet::computeJacobians()");
  }

  /// Set the scalar part of the jacobian
  virtual void computeScalarJacobian(const RealVector& normal,
                              RealVector& jacob) = 0;

  /// Set the matrix of the right eigenvectors and the matrix of the eigenvalues
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                     RealMatrix& leftEv,
                                     RealVector& eValues,
                                     const RealVector& normal) = 0;

  /// Get the model
  Common::SafePtr<LinearAdvTerm> getModel() const
  {
    return _model;
  }

protected:


  /// Compute the convective flux
  virtual void computeFlux(const RealVector& pdata, const RealVector& normals)=0;

  /// Compute the convective flux
  virtual void computeStateFlux(const RealVector& pdata)=0;

  /// acquaintance of the model
  Common::SafePtr<LinearAdvTerm> _model;

}; // end of class LinearAdv2DVarSet

//////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
}
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_LinearAdvVarSet_hh
