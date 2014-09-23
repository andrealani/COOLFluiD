// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NonLinearAdv2DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NonLinearAdv {

//////////////////////////////////////////////////////////////////////////////

CFreal NonLinearAdv2DVarSet::getMaxEigenValueImpl(State& state,
						 const RealVector& normal)
{
 /// @todo broken after release 2009.3
/*  const RealVector& data = *_currData;
  return data[NonLinearAdvTerm::VX]*normal[XX] +
         data[NonLinearAdvTerm::VY]*normal[YY];*/
 return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NonLinearAdv2DVarSet::getMaxAbsEigenValueImpl(State& state,
             const RealVector& normal)
{
 /// @todo broken after release 2009.3
/*  const RealVector& data = *_currData;
  return fabs(data[NonLinearAdvTerm::VX]*normal[XX] +
         data[NonLinearAdvTerm::VY]*normal[YY]);*/
 return 0.0;
}

//////////////////////////////////////////////////////////////////////////////

void NonLinearAdv2DVarSet::setEigenValuesImpl(State& state,
					     const RealVector& normal,
					     RealVector& result)
{
  cf_assert(state.size() == 1);
  cf_assert(result.size() == 1);
  /// @todo broken after release 2009.3
/*  const RealVector& data = *_currData;

  result[0] = data[NonLinearAdvTerm::VX]*normal[XX] +
              data[NonLinearAdvTerm::VY]*normal[YY];*/
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NonLinearAdv

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
