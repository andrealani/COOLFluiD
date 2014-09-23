// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CheckNodeNumberingTetra.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

const CFuint CheckNodeNumberingTetra::nbNodesPerElem = 4;

const CFuint CheckNodeNumberingTetra::dim = 3;

//////////////////////////////////////////////////////////////////////////////

CheckNodeNumberingTetra::CheckNodeNumberingTetra(const Table<CFreal>* const nodes) :
  CheckNodeNumbering(nodes)
{
  _nodalCoord.resize(nbNodesPerElem,dim);
}

//////////////////////////////////////////////////////////////////////////////

CheckNodeNumberingTetra::~CheckNodeNumberingTetra()
{
}

//////////////////////////////////////////////////////////////////////////////

void CheckNodeNumberingTetra::
checkElementNodalNumberingImpl()
{
  CFAUTOTRACE;

  const CFuint nbElem = _elementNode->nbRows();

  // loop over all elements in read from THOR file
  for (CFuint iElem = 0; iElem < nbElem; ++iElem) {

    CFLogDebugMax("Checking Element: " << iElem << "\n");

    cf_assert(_elementNode->nbCols(iElem) == nbNodesPerElem);

    copyCoord<nbNodesPerElem,dim>(iElem);

    /// The volume is calculated and if it is negative then nodes
    /// 1 and 2 are swapped in the element-node connectivity.
    /// This is a definite solution hence it is not checked afterwards.
    CFreal positive = _calculator.checkTetraNumbering(_nodalCoord);
    if(!positive) {
      CFLogDebugMax( "Negative Volume - Swapping nodes and states in tetrahedra." << "\n");

      swap((*_elementNode )(iElem,1),(*_elementNode )(iElem,2));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace THOR2CFmesh

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

