// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_GlobalTypeTrait_hh
#define COOLFluiD_Framework_GlobalTypeTrait_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Framework {

//////////////////////////////////////////////////////////////////////////////

template <typename T>
class GlobalTypeTrait
{
public:
  typedef T GTYPE;
};

template <>
class GlobalTypeTrait<State*>
{
public:
  typedef State::GTYPE GTYPE;
};

template <>
class GlobalTypeTrait<Node*>
{
public:
  typedef Node::GTYPE GTYPE;
};

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} //  namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_GlobalTypeTrait_hh
