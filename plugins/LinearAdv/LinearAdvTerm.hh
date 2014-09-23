// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_LinearAdv_LinearAdvTerm_hh
#define COOLFluiD_Physics_LinearAdv_LinearAdvTerm_hh

////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseTerm.hh"

////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
namespace LinearAdv {

////////////////////////////////////////////////////////////////////////////

// This class represents the interface for a LinearAdvTerm.
// @author Andrea Lani
// @author Tiago Quintino
class LinearAdvTerm : public Framework::BaseTerm {

public:

  // Enumerator defining the mapping between
  // the variable name and its position in the
  // physical data
  enum {u=0, VX=1, VY=2, VZ=3};

  // Constructor without arguments
  LinearAdvTerm(const std::string& name);

  // Default destructor
  virtual ~LinearAdvTerm();

  // Set physical data
  virtual void setupPhysicalData();

  // Resize the physical data
  virtual void resizePhysicalData(RealVector& physicalData);

  // Physical data size
  virtual CFuint getDataSize() const
  {
    return 4;
  }

  // Get the name
  static std::string getName()
  {
    return "LinearAdvTerm";
  }

}; // end of class LinearAdvTerm

////////////////////////////////////////////////////////////////////////////

} // namespace LinearAdv
  }
} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinearAdv_LinearAdvTerm_hh
