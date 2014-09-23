// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.



#include "Framework/MeshStatistics.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

MeshStatistics::MeshStatistics() :
  m_nbFaces(0),
  m_maxNbStatesInCell(0),
  m_maxNbNodesInCell(0),
  m_maxNbFacesInCell(0)
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshStatistics::reset()
{
  m_nbFaces           = 0;
  m_maxNbStatesInCell = 0;
  m_maxNbNodesInCell  = 0;
  m_maxNbFacesInCell  = 0;
}

//////////////////////////////////////////////////////////////////////////////

void MeshStatistics::setMaxNbStatesInCell(const CFuint maxNbStatesInCell)
{
  m_maxNbStatesInCell =
    std::max(maxNbStatesInCell, m_maxNbStatesInCell) ;
}

//////////////////////////////////////////////////////////////////////////////

void MeshStatistics::setMaxNbNodesInCell(const CFuint maxNbNodesInCell)
{
  m_maxNbNodesInCell =
    std::max(maxNbNodesInCell, m_maxNbNodesInCell) ;
}

//////////////////////////////////////////////////////////////////////////////

void MeshStatistics::setMaxNbFacesInCell(const CFuint maxNbFacesInCell)
{
  m_maxNbFacesInCell =
    std::max(maxNbFacesInCell, m_maxNbFacesInCell);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

