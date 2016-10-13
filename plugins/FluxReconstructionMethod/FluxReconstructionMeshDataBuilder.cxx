// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstructionMeshDataBuilder.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {


//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider< FluxReconstructionMeshDataBuilder,MeshDataBuilder,FluxReconstructionModule,1 >
  fluxReconstructionMeshDataBuilderProvider("FluxReconstruction");

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionMeshDataBuilder::FluxReconstructionMeshDataBuilder(const std::string& name) :
  MeshDataBuilder(name)
{
}

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionMeshDataBuilder::~FluxReconstructionMeshDataBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionMeshDataBuilder::setMaxNbStatesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbStatesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionMeshDataBuilder::setMaxNbNodesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbNodesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionMeshDataBuilder::setMaxNbFacesInCell()
{
  MeshDataStack::getActive()->Statistics().setMaxNbFacesInCell(
    PhysicalModelStack::getActive()->getNbEq()+1 );
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionMeshDataBuilder::setMapGeoToTrs()
{
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

