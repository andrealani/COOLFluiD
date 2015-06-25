// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullComputeNorm_hh
#define COOLFluiD_Framework_NullComputeNorm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeNorm.hh"
#include "Common/CFLog.hh"
#include "Framework/GlobalReduce.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class is a Null object of the functor class that computes
/// the Norm the SubSystem.
/// @author Tiago Quintino
class Framework_API NullComputeNorm : public ComputeNorm {
public:

  /// For global reduce
  typedef CFreal GR_RESULTTYPE;

  /// Default constructor without arguments
  /// @see ComputeNorm()
  NullComputeNorm(const std::string& name) : ComputeNorm(name), _gr(*this)
  {
    CFLog(NOTICE,"NullComputeNorm() name=" << name << "\n");
  }

  /// Default destructor
  ~NullComputeNorm() {}

  /// Calculate the norm
  RealVector compute ();

  /// For globalreduce
  CFreal GR_GetLocalValue () const;

  /// for globalreduce
  static void GR_Combine (const CFreal & S1, const CFreal & S2, CFreal & S3);

private:

  Framework::GlobalReduce<NullComputeNorm> _gr;
  
}; // end of class NullComputeNorm

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullComputeNorm_hh
