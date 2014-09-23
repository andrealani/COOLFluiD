// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Burgers2DVarSet.hh"
#include "BurgersTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Burgers {

//////////////////////////////////////////////////////////////////////////////

 CFreal Burgers2DVarSet::getMaxEigenValue(const RealVector& pdata, const RealVector& normal)
{
   return pdata[BurgersTerm::VX]*normal[XX] + normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

 CFreal Burgers2DVarSet::getMaxAbsEigenValue(const RealVector& pdata, const RealVector& normal)
{
  return fabs(pdata[BurgersTerm::VX]*normal[XX] + normal[YY]);
}

//////////////////////////////////////////////////////////////////////////////

void Burgers2DVarSet::computeEigenValues(const RealVector& pdata, const RealVector& normal, RealVector& eValues)
{
   eValues[0] = pdata[BurgersTerm::VX]*normal[XX] + normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Burgers

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
