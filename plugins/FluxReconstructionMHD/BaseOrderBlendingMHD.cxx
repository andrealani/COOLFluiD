// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "MHD/MHDTerm.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BaseOrderBlendingMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BaseOrderBlendingMHD, FluxReconstructionSolverData, FluxReconstructionMHDModule>
    OrderBlendingMHDFRProvider("OrderBlendingMHD");

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlendingMHD::BaseOrderBlendingMHD(const std::string& name) :
  BaseOrderBlending(name)
{
}

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlendingMHD::~BaseOrderBlendingMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlendingMHD::extractMonitoredField()
{
  // Handle MHD-specific "B2" (magnetic pressure) via symbolic physical data slots.
  // MHDTerm::BX/BY/BZ are the canonical indices for all MHD variable sets — no
  // hardcoded state indices and no dependency on conservative vs primitive layout.
  if (m_modalMonitoredExpression == "B2")
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_obUpdateVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_obPData);
      const CFreal bx = m_obPData[MHDTerm::BX];
      const CFreal by = m_obPData[MHDTerm::BY];
      const CFreal bz = m_obPData[MHDTerm::BZ];
      m_tempSolPntVec[iSol] = bx*bx + by*by + bz*bz;
    }
    return;
  }

  // Any other expression is handled by the physics-agnostic base class.
  BaseOrderBlending::extractMonitoredField();
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
