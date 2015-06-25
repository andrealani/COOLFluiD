// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ComputeAllNorms_hh
#define COOLFluiD_Framework_ComputeAllNorms_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalReduce.hh"
#include "Framework/ComputeNorm.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/State.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is a functor computing the L2, L1 and L infinity norms
/// for the rhs of the solution
/// @author Tiago Quintino
class Framework_API ComputeAllNorms : public ComputeNorm {

public: // functions

  /// For global reduce
  typedef CFreal GR_RESULTTYPE;

  /// Default constructor without arguments
  ComputeAllNorms(const std::string& name);

  /// Default destructor
  virtual ~ComputeAllNorms();

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
  Framework::GlobalReduce<ComputeAllNorms> m_gr;
  /// socket for rhs
  Framework::DataSocketSink< CFreal> socket_rhs;
  /// socket for states
  Framework::DataSocketSink<Framework::State*, Framework::GLOBAL> socket_states;

}; // end of class ComputeAllNorms

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ComputeAllNorms_hh
