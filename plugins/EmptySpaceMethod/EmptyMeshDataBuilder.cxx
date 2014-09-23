// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/MeshData.hh"

#include "EmptySpaceMethod/EmptyMeshDataBuilder.hh"
#include "EmptySpaceMethod/Empty.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace EmptySpaceMethod {


//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< EmptyMeshDataBuilder,MeshDataBuilder,EmptyModule,1 >
  emptyMeshDataBuilderProvider("Empty");

//////////////////////////////////////////////////////////////////////////////

EmptyMeshDataBuilder::EmptyMeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
}

//////////////////////////////////////////////////////////////////////////////

EmptyMeshDataBuilder::~EmptyMeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void EmptyMeshDataBuilder::setMaxNbStatesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void EmptyMeshDataBuilder::setMaxNbNodesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void EmptyMeshDataBuilder::setMaxNbFacesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void EmptyMeshDataBuilder::setMapGeoToTrs()
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace EmptySpaceMethod
}  // namespace COOLFluiD

