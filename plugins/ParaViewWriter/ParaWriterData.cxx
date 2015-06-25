// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "ParaViewWriter/ParaViewWriter.hh"
#include "ParaViewWriter/ParaWriterData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace IO {

    namespace ParaViewWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<ParaWriterData>, ParaWriterData, ParaViewWriterModule> NullParaWriterComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void ParaWriterData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("printExtraValues","flag telling if to print extra values");
   options.addConfigOption< std::string >("updateVar","update variables set");
   options.addConfigOption< std::vector<std::string> >("SurfaceTRS","List of TRS's to be writen in the surface file.");
   options.addConfigOption< bool >("SurfaceOnly","Print only the surface data chosen in SurfaceTRS");
   options.addConfigOption< bool >("VectorAsComponents","Switch to write velocity by components or coupled.");
}

//////////////////////////////////////////////////////////////////////////////

ParaWriterData::ParaWriterData(Common::SafePtr<Framework::Method> owner)
  : OutputFormatterData(owner),
    m_filepath(),
    m_updateVarStr(),
    m_updateVarSet(),
    m_stdTrsGeoBuilder()
{
  addConfigOptionsTo(this);
  m_updateVarStr = "Null";
  setParameter("updateVar",&m_updateVarStr);

  m_printExtraValues = false;
  setParameter("printExtraValues",&m_printExtraValues);

  m_surface_only = false;
  setParameter("SurfaceOnly",&m_surface_only);

  m_surfTRS = std::vector<std::string>();
  setParameter("SurfaceTRS",&m_surfTRS);

  m_writeVectorAsComponents=false;
  setParameter("VectorAsComponents",&m_writeVectorAsComponents);
}

//////////////////////////////////////////////////////////////////////////////

ParaWriterData::~ParaWriterData()
{
}

//////////////////////////////////////////////////////////////////////////////

void ParaWriterData::configure ( Config::ConfigArgs& args )
{
  OutputFormatterData::configure(args);

  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string provider = "Null";
  if (m_updateVarStr != "Null") {
    provider = (physModel->getConvectiveName() != "Null") ?
      (physModel->getConvectiveName() + m_updateVarStr) :
      (physModel->getDiffusiveName() + m_updateVarStr) ;
  }

  CFLog(VERBOSE, "ParaWriterData::UpdateVarStr = " << provider << "\n");

  m_updateVarSet.reset(Environment::Factory<ConvectiveVarSet>::getInstance().
                      getProvider(provider)->create(physModel->getImplementor()->getConvectiveTerm()));

  cf_assert(m_updateVarSet.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

CFuint ParaWriterData::getVTKCellTypeID(CFGeoShape::Type shape,CFuint geoOrder)
{
  switch (geoOrder)
  {
    case 0:
    case 1:
    {
      switch (shape)
      {
        case CFGeoShape::LINE:
        {
          return 3;
        }
        case CFGeoShape::TRIAG:
        {
          return 5;
        }
        case CFGeoShape::QUAD:
        {
          return 9;
        }
        case CFGeoShape::TETRA:
        {
          return 10;
        }
        case CFGeoShape::PYRAM:
        {
          return 14;
        }
        case CFGeoShape::PRISM:
        {
          return 13;
        }
        case CFGeoShape::HEXA:
        {
          return 12;
        }
        default:
        {
          throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
        }
      }
    } break;
    case 2:
    {
      switch (shape)
      {
        case CFGeoShape::LINE:
        {
          return 21;
        }
        case CFGeoShape::TRIAG:
        {
          return 22;
        }
        case CFGeoShape::QUAD:
        {
          return 23;
        }
        case CFGeoShape::TETRA:
        {
          return 24;
        }
        case CFGeoShape::HEXA:
        {
          return 25;
        }
        default:
        {
          throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
        }
      }
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

vector< RealVector > ParaWriterData::getOutputPntsMappedCoords(CFGeoShape::Type shape,CFuint solOrder)
{
  // output variable
  vector< RealVector > nodeMappedCoords;

  // for zeroth order solution polynomial, use same points as for first order
  if (solOrder == 0)
  {
    solOrder = 1;
  }

  // number of points needed for representing a polynomial of order degree solOrder
  const CFuint nbrNodes1D = solOrder + 1;

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        RealVector coords(1);
        coords[KSI] = -1.0 + iKsi*2.0/solOrder;
        nodeMappedCoords.push_back(coords);
      }
    } break;
    case CFGeoShape::TRIAG:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = iKsi*1.0/solOrder;
        for (CFuint iEta = 0; iEta < nbrNodes1D-iKsi; ++iEta)
        {
          RealVector coords(2);
          coords[KSI] = ksi;
          coords[ETA] = iEta*1.0/solOrder;
          nodeMappedCoords.push_back(coords);
        }
      }
    } break;
    case CFGeoShape::QUAD:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = -1.0 + iKsi*2.0/solOrder;
        for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
        {
          RealVector coords(2);
          coords[KSI] = ksi;
          coords[ETA] = -1.0 + iEta*2.0/solOrder;
          nodeMappedCoords.push_back(coords);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      /// @warn: for tetra this is only implemented for P1
      RealVector coords(3);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 1.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 1.0;
      coords[ZTA] = 0.0;
      nodeMappedCoords.push_back(coords);
      coords[KSI] = 0.0;
      coords[ETA] = 0.0;
      coords[ZTA] = 1.0;
      nodeMappedCoords.push_back(coords);
    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputPntsMappedCoords() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      for (CFuint iKsi = 0; iKsi < nbrNodes1D; ++iKsi)
      {
        const CFreal ksi = -1.0 + iKsi*2.0/solOrder;
        for (CFuint iEta = 0; iEta < nbrNodes1D; ++iEta)
        {
          const CFreal eta = -1.0 + iEta*2.0/solOrder;
          for (CFuint iZta = 0; iZta < nbrNodes1D; ++iZta)
          {
            RealVector coords(3);
            coords[KSI] = ksi;
            coords[ETA] = eta;
            coords[ZTA] = -1.0 + iZta*2.0/solOrder;
            nodeMappedCoords.push_back(coords);
          }
        }
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }

  return nodeMappedCoords;
}

//////////////////////////////////////////////////////////////////////////////

vector< vector< CFuint > > ParaWriterData::getOutputCellNodeConn(CFGeoShape::Type shape,CFuint solOrder)
{
  // output variable
  vector< vector< CFuint > > cellsNodesConn;

  // for zeroth order solution polynomial, use same points as for first order
  if (solOrder == 0)
  {
    solOrder = 1;
  }

  // number of points needed for representing a polynomial of order degree solOrder
  const CFuint nbrNodes1D = solOrder + 1;

  switch (shape)
  {
    case CFGeoShape::LINE:
    {
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        vector< CFuint > cellNodesConn(2);
        cellNodesConn[0] = iKsi;
        cellNodesConn[1] = iKsi+1;
        cellsNodesConn.push_back(cellNodesConn);
      }
    } break;
    case CFGeoShape::TRIAG:
    {
      CFuint nodeIdx = 0;
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < solOrder-iKsi-1; ++iEta)
        {
          vector< CFuint > cellNodesConn(3);
          cellNodesConn[0] = nodeIdx + iEta;
          cellNodesConn[1] = nodeIdx + nbrNodes1D - iKsi + iEta;
          cellNodesConn[2] = nodeIdx + iEta + 1;
          cellsNodesConn.push_back(cellNodesConn);

          cellNodesConn[0] = nodeIdx + nbrNodes1D - iKsi + iEta;
          cellNodesConn[1] = nodeIdx + nbrNodes1D - iKsi + iEta + 1;
          cellNodesConn[2] = nodeIdx + iEta + 1;
          cellsNodesConn.push_back(cellNodesConn);
        }
        
        vector< CFuint > cellNodesConn(3);
        cellNodesConn[0] = nodeIdx + solOrder-iKsi-1;
        cellNodesConn[1] = nodeIdx + nbrNodes1D -iKsi + solOrder - 1 - iKsi;
        cellNodesConn[2] = nodeIdx + solOrder-iKsi;
        cellsNodesConn.push_back(cellNodesConn);
        
        nodeIdx += nbrNodes1D-iKsi;
      }
    } break;
    case CFGeoShape::QUAD:
    {
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < solOrder; ++iEta)
        {
          vector< CFuint > cellNodesConn(4);
          cellNodesConn[0] = (iKsi  )*nbrNodes1D + iEta  ;
          cellNodesConn[1] = (iKsi+1)*nbrNodes1D + iEta  ;
          cellNodesConn[2] = (iKsi+1)*nbrNodes1D + iEta+1;
          cellNodesConn[3] = (iKsi  )*nbrNodes1D + iEta+1;
          cellsNodesConn.push_back(cellNodesConn);
        }
      }
    } break;
    case CFGeoShape::TETRA:
    {
      /// @warn: for tetra this is only implemented for P1
      vector< CFuint > cellNodesConn(4);
      cellNodesConn[0] = 0;
      cellNodesConn[1] = 1;
      cellNodesConn[2] = 2;
      cellNodesConn[3] = 3;
      cellsNodesConn.push_back(cellNodesConn);
    } break;
    case CFGeoShape::PYRAM:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputCellNodeConn() for PYRAM");
    } break;
    case CFGeoShape::PRISM:
    {
      throw Common::NotImplementedException (FromHere(),"ParaWriterData::getOutputCellNodeConn() for PRISM");
    } break;
    case CFGeoShape::HEXA:
    {
      const CFuint nbrNodes1DSq = nbrNodes1D*nbrNodes1D;
      for (CFuint iKsi = 0; iKsi < solOrder; ++iKsi)
      {
        for (CFuint iEta = 0; iEta < solOrder; ++iEta)
        {
          for (CFuint iZta = 0; iZta < solOrder; ++iZta)
          {
            vector< CFuint > cellNodesConn(8);
            cellNodesConn[0] = (iKsi  )*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta  ;
            cellNodesConn[1] = (iKsi+1)*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta  ;
            cellNodesConn[2] = (iKsi+1)*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta  ;
            cellNodesConn[3] = (iKsi  )*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta  ;
            cellNodesConn[4] = (iKsi  )*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta+1;
            cellNodesConn[5] = (iKsi+1)*nbrNodes1DSq + (iEta  )*nbrNodes1D + iZta+1;
            cellNodesConn[6] = (iKsi+1)*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta+1;
            cellNodesConn[7] = (iKsi  )*nbrNodes1DSq + (iEta+1)*nbrNodes1D + iZta+1;
            cellsNodesConn.push_back(cellNodesConn);
          }
        }
      }
    } break;
    default:
    {
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown GeoShape");
    }
  }

  return cellsNodesConn;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace ParaViewWriter

  } // namespace IO

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

