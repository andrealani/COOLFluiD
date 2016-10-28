// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeNorm_hh
#define COOLFluiD_Framework_ComputeNorm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/NumericalStrategy.hh"
#include "Framework/Storage.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a generic interface for computing the norm
/// @author Andrea Lani
/// @author Tiago Quintino
class Framework_API ComputeNorm : public NumericalStrategy {

public: // typedefs

  typedef Environment::ConcreteProvider<ComputeNorm,1> PROVIDER;
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  ComputeNorm(const std::string& name);

  /// Default destructor
  virtual ~ComputeNorm();

  /// Calculates the norms
  virtual RealVector compute () = 0;

  /// Setup the object
  virtual void setup();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Get the ID's of all the variables whose norm is computed
  std::vector<CFuint> getComputedNormVarIDs()
  {
    return m_compute_var_id;
  }

  /// Get the ID of the variable to monitor
  CFuint getMonitoredVarIndex()
  {
    return m_var_idx;
  }

  /// Get the number of residuals returned
  CFuint getNbComputedResiduals()
  {
    return m_compute_var_id.size();
  }

  /// Get the flag telling that you use the normalized global residual 
  // for the convergence method
  bool getGlobalRes()
  {
    return m_global_res;
  }
  
  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeNorm";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// index of the varID in the m_compute_var_id vector
  CFuint m_var_idx;
  /// ID of the variable of which to monitor the norm
  CFuint m_var;
  /// list of ID's of the variables of which to compute the norm
  std::vector<CFuint> m_compute_var_id;
  /// values of the computed residuals
  RealVector m_residuals;
  /// temporay index iterator
  CFuint m_var_itr;
  /// Flag to normalize the residuals
  bool m_normalizedRes;
  /// Values to normalized the residual
  std::vector<CFreal> m_refVals;
  /// Flag to use the global residual for the convergece method
  bool m_global_res;

}; // end of class ComputeNorm

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(ComputeNorm) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeNorm_hh
