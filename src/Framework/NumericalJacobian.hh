// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NumericalJacobian_hh
#define COOLFluiD_Framework_NumericalJacobian_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "Config/ConfigObject.hh"
#include "Framework/PhysicalModel.hh"
#include "MathTools/MathChecks.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/MathFunctions.hh"

#ifdef CF_HAVE_CUDA
#include "Common/CUDA/CudaEnv.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class computes the variable perturbation and the finite difference
/// for evaluating numerical jacobian.
/// @author Andrea Lani
class Framework_API NumericalJacobian : public Config::ConfigObject, 
					public Common::NonCopyable<NumericalJacobian> {
public:

#ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename PHYS>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE ~DeviceConfigOptions() {}
    
    CFreal originalValue; /// backup copy of the original value
    CFreal refValues[PHYS::NBEQS]; /// array of reference Values 
    CFreal tol; /// tolerance
  };
  
  /// nested class defining a functor
  template <typename PHYS>
  class DeviceFunc {
  public:
    typedef NumericalJacobian BASE;
        
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<PHYS>* dco) : m_dco(dco) {}
    
    /// Compute eps for the numericsl evaluation of the jacobian
    /// and perturb the given component of the state vector
    HOST_DEVICE void perturb(const CFuint iVar, CFreal* value)
    {
      m_dco->originalValue = *value; setEps(iVar, *value); *value += m_eps;
    }    
    
    /// Restore the original value in the perturbed component
    HOST_DEVICE void restore(CFreal* value) const {*value = m_dco->originalValue;}
    
    /// Compute the derivative
    template <typename A1, typename A2, typename A3>
    HOST_DEVICE void computeDerivative(A1* res, A2* pertRes, A3* diff) const
    {
      const CFreal invEps = 1./m_eps; 
      (*diff) = (*pertRes - *res)*invEps;
    }
    
  private:
    
    /// @return 1.0 if value >= 0.0
    HOST_DEVICE CFreal sign(const CFreal value) {return (value < 0.0) ? -1.0 : 1.0;}
    
    /// compute epsilon
    HOST_DEVICE void setEps(const CFuint iVar, const CFreal value)
    {
      const CFreal absv = (value>=0.) ? value : -value; 
      const CFreal absr = (m_dco->refValues[iVar]>0.) ? m_dco->refValues[iVar] : -m_dco->refValues[iVar];
      const CFreal maxa = (absv > absr) ? absv : absr; 
      m_eps = m_dco->tol*sign(value)*maxa;
    }
    
    /// device options
    DeviceConfigOptions<PHYS>* m_dco;
    
    /// epsilon for numerical perturbation
    CFreal m_eps;
  };
  
  /// copy the local configuration options to the device
  template <typename PHYS>
  void copyConfigOptionsToDevice(DeviceConfigOptions<PHYS>* dco) 
  {
    CudaEnv::copyHost2Dev(&dco->originalValue, &_originalValue, 1);
    cf_always_assert(_refValues.size() == PHYS::NBEQS);
    CudaEnv::copyHost2Dev(dco->refValues, &_refValues[0], PHYS::NBEQS);
    CudaEnv::copyHost2Dev(&dco->tol, &_tol, 1);
  }   
  
  /// copy the local configuration options to the device
  template <typename PHYS>
  void copyConfigOptions(DeviceConfigOptions<PHYS>* dco) 
  {
    dco->originalValue = _originalValue;
    for (CFuint i = 0; i < PHYS::NBEQS; ++i) {dco->refValues[i] = _refValues[i];}
    dco->tol = _tol;
  }   
#endif
  
  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  NumericalJacobian(const std::string& name);

  /// Compute eps for the numericsl evaluation of the jacobian
  /// and perturb the given component of the state vector
  void perturb(const CFuint iVar, CFreal& value)
  {
    // store the original value of the to-be-perturbed
    // component of the state vector
    _originalValue = value;

    // compute eps for this component
    setEps(iVar, value);
    value += _eps;
  }

  /// Restore the original value in the perturbed component
  void restore(CFreal& value) const
  {
    value = _originalValue;
  }

  /// Set the reference values
  void setRefValues(RealVector& refValues)
  {
    const CFuint size = refValues.size();
    _refValues.resize(size);
    _refValues = refValues;
  }

  /// Compute the derivative
  template <class ARRAY>
  void computeDerivative(const ARRAY& res,
			 const ARRAY& pertRes,
			 ARRAY& diff) const
  {
    cf_assert(MathTools::MathChecks::isNotZero(_eps));
    const CFreal invEps = 1./_eps;
    diff = (pertRes - res)*invEps;
  }
  
  /// Compute the derivative
  template <class ARRAY>
  void computeDerivativePrint(const ARRAY& res,
			      const ARRAY& pertRes,
			      ARRAY& diff) const
  {
    cf_assert(MathTools::MathChecks::isNotZero(_eps));
    const CFreal invEps = 1./_eps;
    diff = (pertRes - res)*invEps;
    
    /* if (res.size() == 13) {
       CFLog(INFO, "pertRes = " << pertRes << "\n");
      CFLog(INFO, "res     = " << res << "\n");
      CFLog(INFO, "diff    = " << diff << "\n");
      CFLog(INFO, "_eps    = " << _eps << ", invEps = " << invEps << " \n");
      
      diff = pertRes - res;
      CFLog(INFO, "diff 2  = " << diff << "\n");
      diff *= invEps;
      CFLog(INFO, "diff 3  = " << diff << "\n");
      }*/
  }
  
  /// Compute the derivative for a slice
  void computeDerivative(RealSliceVector res,
			 RealSliceVector pertRes,
			 RealSliceVector diff) const
  {
    cf_assert(MathTools::MathChecks::isNotZero(_eps));
    const CFreal invEps = 1./_eps;
    diff = (pertRes - res)*invEps;

  }
  
  /// Get epsilon
  CFreal getEps() const
  {
    return _eps;
  }
  
  /// compute epsilon
  CFreal computeEps(const CFuint iVar, const CFreal value) const
  {
    return _tol*MathTools::MathFunctions::sign(value)*
      std::max(std::abs(value), std::abs(_refValues[iVar]));
  }
  
private: // helper method

  /// Copy constructor
  NumericalJacobian(const NumericalJacobian& other);

  /// Overloading of the assignment operator
  const NumericalJacobian& operator=
    (const NumericalJacobian& other);

  /// Set epsilon
  void setEps(const CFuint iVar, const CFreal value)
  {
    _eps = _tol*MathTools::MathFunctions::sign(value)*
      std::max(std::abs(value), std::abs(_refValues[iVar]));

    // std::max(std::abs(value), MathTools::MathConsts::CFrealEps());

    // _eps = MathTools::MathFunctions::sign(value)*
    //       std::max(MathTools::MathConsts::CFrealEps(), _tol*std::abs(value));
  }

private: // data

  /// perturbation value
  CFreal _eps;

  /// backup copy of the original value
  CFreal _originalValue;

  /// RealVector of reference Values
  RealVector _refValues;

  // tolerance for eps computation
  CFreal _tol;

}; // end of class NumericalJacobian

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NumericalJacobian_hh
