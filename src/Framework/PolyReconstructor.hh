// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PolyReconstructor_hh
#define COOLFluiD_Framework_PolyReconstructor_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Storage.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/CFSide.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for all polynomial reconstructors
/// for FVM
/// @author Andrea Lani
template < typename METHODDATA >
class PolyReconstructor : public MethodStrategy<METHODDATA> {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider
  <METHODDATA,PolyReconstructor<METHODDATA> > PROVIDER;

  /// Constructor
  PolyReconstructor(const std::string& name);

  /// Default destructor
  virtual ~PolyReconstructor();

  /// Compute the gradients
  virtual void computeGradients() = 0;

  /// Set private data that will be used during the computation
  virtual void setup();
  
  /// Unsetup private data that will be used during the computation
  virtual void unsetup();
  
  /**
   * Prepare the reconstruction step
   */
  virtual void prepareReconstruction() {}
  
  /// Set the flag imposing a null gradient
  void setZeroGradient(std::vector<bool> *const zeroGradient)
  {
    cf_assert(zeroGradient != CFNULL);
    _zeroGradient = zeroGradient;
  }
  
  /// Extrapolate the solution in the quadrature points of
  /// the given face
  void extrapolate(GeometricEntity* const face)
  {
    _isBoundaryFace = face->getState(1)->isGhost();
    extrapolateImpl(face);
  }
  
  /// Extrapolate the solution in the quadrature points of
  /// the given face
  /// @param face         pointer to the face object
  /// @param iVar         ID of the variable to reconstruct
  /// @param leftOrRight  0(left) or 1(right)
  void extrapolate(GeometricEntity* const face, CFuint iVar, CFuint leftOrRight)
  {
    _isBoundaryFace = face->getState(1)->isGhost();
    extrapolateImpl(face, iVar, leftOrRight);
  }
  
  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );
  
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "PolyReconstructor";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// Restore the back up values corresponding to the given variable ID
  /// AL: don't make this function virtual: the profiler indicates that
  /// the high call frequency makes this function expensive if it is
  /// not inlined!!
  void restoreValues(CFuint iVar, CFuint leftOrRight)
  {
    getValues(leftOrRight)[iVar] = getBackupValues(leftOrRight)[iVar];
  }
  
  /// Update the weights when nodes are moving
  virtual void updateWeights()
  {
    throw Common::NotImplementedException (FromHere(),"PolyReconstructor::updateWeights()");
  }
  
  /// Get the extrapolated values
  std::vector<State*>& getExtrapolatedValues() {return _extrapValues;}
  
  /// Get the extrapolated physical data
  std::vector<RealVector>& getExtrapolatedPhysicaData() {return _extrapPhysData;}
  
  /// Get the backup physical data
  std::vector<RealVector>& getBackupPhysicaData() {return _backupPhysData;}
  
protected: // helper functions
  
  /// Allocate reconstruction data needed for the flux evaluation
  virtual void allocateReconstructionData() = 0;
  
  /// Copy the backup values LEFT and RIGHT
  void copyBackupValues()
  {
    getValues(LEFT) = getBackupValues(LEFT);
    getValues(RIGHT) = getBackupValues(RIGHT);
  }
  
  /// Tells if the face is on the boundary
  bool isBoundaryFace() const
  {
    return _isBoundaryFace;
  }
  
  /// Extrapolate the solution in the quadrature points of
  /// the given face
  virtual void extrapolateImpl(GeometricEntity* const face) = 0;
  
  /// Extrapolate the solution in the quadrature points of
  /// the given face
  virtual void extrapolateImpl(GeometricEntity* const face,
			       CFuint iVar,
			       CFuint leftOrRight) = 0;
  
  /// Get the left extrapolated solution values
  State& getValues(CFuint leftOrRight) {return *_extrapValues[leftOrRight];}
  
  /// Get the back up extrapolated solution values
  RealVector& getBackupValues(CFuint leftOrRight) {return _backupValues[leftOrRight];}
     
  /// Get the extrapolated coordinates
  const std::vector<Node*>& getCoord() const {return _extrapCoord;}
  
private: // data

  /// this flag tell if the processed face is on the boundary
  bool _isBoundaryFace;
  
protected:
  
  /// array of extrapolated coordinates
  std::vector<Node*> _extrapCoord;
  
  /// array of extrapolated solution values
  std::vector<State*> _extrapValues;
  
  /// back up left values
  std::vector<RealVector> _backupValues;
  
  /// array of extrapolated physical data values
  std::vector<RealVector> _extrapPhysData;
  
  /// back up of extrapolated physical data values
  std::vector<RealVector> _backupPhysData;
  
  /// socket for limiter values
  Framework::DataSocketSink< CFreal> socket_limiter;
  
  /// flag imposing a null gradient
  std::vector<bool>* _zeroGradient;

  // iteration threshold after which historical modification
  /// is applied to the limiter
  CFuint _limitIter;

  /// residual threshold after which historical modification
  /// is applied to the limiter
  CFreal _limitRes;

  /// interactive factor to multiply the gradient
  CFreal _gradientFactor;
  
  /// flag telling iof to freeze the limiter value
  bool _freezeLimiter;

  /// apply a fix for the BC treatment
  bool _bcFix;

}; // end of class PolyReconstructor

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

#include "PolyReconstructor.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PolyReconstructor_hh
