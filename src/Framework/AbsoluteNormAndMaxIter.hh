// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_AbsoluteNormAndMaxIter_hh
#define COOLFluiD_Framework_AbsoluteNormAndMaxIter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StopCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a StopConditon for the convergence process
/// that instructs it to stop when a certain norm of the residual is reached.
/// @author Tiago Quintino
class Framework_API AbsoluteNormAndMaxIter : public StopCondition {

public: // methods

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments
  /// @see StopCondition()
  AbsoluteNormAndMaxIter(const std::string& name);

  /// Default destructor
  ~AbsoluteNormAndMaxIter();

  /// returns true if we have to evaluate the stopcondition
  /// needs to combine values from different CPU domains
  virtual bool IsGlobal () const;

  /// Take the combined value from all CPU's and decide if the SubSystem
  /// should stop
  virtual bool isAchieved (const ConvergenceStatus& status) ;

private: // data

  /// target value of the relative norm
  CFreal m_norm;

  /// maximum number of iterations
  CFuint m_maxIter;

}; // end of class AbsoluteNormAndMaxIter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_AbsoluteNormAndMaxIter_hh
