// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullVarSet.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullVarSet,
			    ConvectiveVarSet,
			    FrameworkLib, 1> 
nullVarSetProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullVarSet::NullVarSet(Common::SafePtr<Framework::BaseTerm> term) : ConvectiveVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

NullVarSet::~NullVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computeJacobians()
{
  CFLog(VERBOSE, "NullVarSet::computeJacobians() => this is a Null Var"
	<< "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computeEigenValuesVectors(RealMatrix& rightEv,
				       RealMatrix& leftEv,
				       RealVector& eValues,
				       const RealVector& normal)
{
  CFLog(VERBOSE,"NullVarSet::computeEigenValuesVectors() => "
	<< "this is a Null VarSet" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computeEigenValues(const RealVector& pdata,
				    const RealVector& normal,
				    RealVector& eValues)
{
  CFLog(VERBOSE,"NullVarSet::computeEigenValues() => "
	<< "this is a Null VarSet" << "\n");
}

//////////////////////////////////////////////////////////////////////////////

CFreal NullVarSet::getMaxEigenValue(const RealVector& pdata,
				    const RealVector& normal)
{
  CFLog(VERBOSE,"NullVarSet::getMaxEigenValue() => "
	<< "this is a Null VarSet" << "\n");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

CFreal NullVarSet::getMaxAbsEigenValue(const RealVector& pdata,
				       const RealVector& normal)
{
  CFLog(VERBOSE,"NullVarSet::getMaxAbsEigenValue() => "
	<< "this is a Null VarSet" << "\n");
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computeFlux(const RealVector& pdata,
			     const RealVector& normals)
{
  CFLog(VERBOSE,"NullVarSet::computeFlux() \n");
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computeStateFlux(const RealVector& pdata)
{
  CFLog(VERBOSE,"NullVarSet::computeFlux() \n");
}

//////////////////////////////////////////////////////////////////////////////

void NullVarSet::computePerturbedPhysicalData(const Framework::State& state,
					      const RealVector& pdataBkp,
					      RealVector& pdata,
					      CFuint iVar)
{
  CFLog(VERBOSE,"NullVarSet::computePerturbedPhysicalData() \n");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Framework
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
