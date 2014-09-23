// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CheckNodeNumberingPyram.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/CFLog.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace THOR2CFmesh {

//////////////////////////////////////////////////////////////////////////////

const CFuint CheckNodeNumberingPyram::nbNodesPerElem = 5;

const CFuint CheckNodeNumberingPyram::dim = 3;

//////////////////////////////////////////////////////////////////////////////

CheckNodeNumberingPyram::CheckNodeNumberingPyram(const Table<CFreal>* const nodes) :
  CheckNodeNumbering(nodes)
{
  _nodalCoord.resize(nbNodesPerElem,dim);
}

//////////////////////////////////////////////////////////////////////////////

CheckNodeNumberingPyram::~CheckNodeNumberingPyram()
{
}

//////////////////////////////////////////////////////////////////////////////

void CheckNodeNumberingPyram::checkElementNodalNumberingImpl()
{
  CFLogDebugMin(
        "CheckNodeNumberingPyram::checkElementNode() => start " << "\n");

  const CFuint nbElem = _elementNode->nbRows();

  // loop over all elements in read from THOR file
  for (CFuint iElem = 0; iElem < nbElem; ++iElem) {

    CFLogDebugMax("Checking Element: " << iElem << "\n");

    cf_assert(_elementNode->nbCols(iElem) == nbNodesPerElem);

    copyCoord<nbNodesPerElem,dim>(iElem);

    /// The volume is calculated and if some partial part
    /// is negative then nodes
    ///  1 and 3 are swapped in the element-node connectivity.
    /// This might not be a definite solution but
    /// we only apply this one. We calculate again and if it didnt solve
    /// another exception is thrown and we bail out...
    bool positive = _calculator.checkPyramNumbering(_nodalCoord);

    if (!positive) {
      CFLogDebugMax("Negative Volume - Swapping nodes and states in pyramid." << "\n");

      swap((*_elementNode )(iElem,1),(*_elementNode )(iElem,3));

      // need to recopy the nodes into temporary matrix
      copyCoord<nbNodesPerElem,dim>(iElem);

      positive = _calculator.checkPyramNumbering(_nodalCoord);
      if(!positive){
        throw Framework::NegativeVolumeException (FromHere(),"Found negative partial volume in Pyramid.");
      }
    }
  }

  CFLogDebugMin("CheckNodeNumberingPyram::checkElementNode() done" <<  "\n");
}

//////////////////////////////////////////////////////////////////////////////

     } // namespace THOR2CFmesh

 

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
