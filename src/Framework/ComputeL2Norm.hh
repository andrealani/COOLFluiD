// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeL2Norm_hh
#define COOLFluiD_Framework_ComputeL2Norm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalReduce.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DynamicDataSocketSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is computes the L2 norm of a vector
/// @author Andrea Lani
/// @author Dries Kimpe
/// @author Tiago Quintino
class Framework_API ComputeL2Norm : public ComputeNorm {

public: // typedefs

  /// For global reduce
  typedef CFreal GR_RESULTTYPE;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  ComputeL2Norm(const std::string& name);

  /// Default destructor
  virtual ~ComputeL2Norm();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Calculates the norms
  virtual RealVector compute ();

  /// Retrieves the value for the global reduce of the result
  CFreal GR_GetLocalValue () const;

  /// Global reduce of the result
  static void GR_Combine (const CFreal & S1, const CFreal & S2, CFreal & S3);

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

private: // data

  /// object to perform the global reduce
  Framework::GlobalReduce<ComputeL2Norm> m_gr;
  /// The set of data sockets to be used by the strategy
  Framework::DynamicDataSocketSet<> sockets_norm;
  /// socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;
  /// name of the vector on which to apply the norm
  std::string m_vecnorm_name;
  
  /// tolerance on the residual
  CFreal m_tolerance;
  
}; // end of class ComputeL2Norm

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeL2Norm_hh
