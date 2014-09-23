// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/NullDomainModel.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullDomainModel,
                            Framework::DomainModel,
                            FrameworkLib, 1>
aNullDomainModelProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullDomainModel::NullDomainModel(const std::string& name) :
  Framework::DomainModel(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullDomainModel::~NullDomainModel()
{
}

//////////////////////////////////////////////////////////////////////////////

Framework::DomainModel::TRidx NullDomainModel::getNbTopoDefs () const
{
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void NullDomainModel::computeParamCoord (const TRidx idx, const XVector& coord , PVector& pcoord) const
{
  CFLogWarn ("Calling computeParamCoord() on a Null DomainModel\n");
}

//////////////////////////////////////////////////////////////////////////////


void NullDomainModel::computeCoord (const TRidx idx, const PVector& pcoord, XVector& coord) const
{
  CFLogWarn ("Calling computeCoord() on a Null DomainModel\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullDomainModel::compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const
{
  CFLogWarn ("Calling compute1stDeriv() on a Null DomainModel\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullDomainModel::compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const
{
  CFLogWarn ("Calling compute2ndDeriv() on a Null DomainModel\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullDomainModel::computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const
{
  CFLogWarn ("Calling computeAll() on a Null DomainModel\n");
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
