// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethod.hh"


#include "ALE_FVMGeometricAverage.hh"
#include "Framework/Node.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ALE_FVMGeometricAverage, NewtonIteratorData, NewtonMethodModule> ale_FVMGeometricAverageProvider("ALE_FVMGeometricAverage");

//////////////////////////////////////////////////////////////////////////////

void ALE_FVMGeometricAverage::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();

  // Set the futureNodes = Nodes
  cf_assert (nodes.size() == futureNodes.size());
  for (CFuint i=0; i < nodes.size();++i){
    *nodes[i] = *futureNodes[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ALE_FVMGeometricAverage::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_futureNodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
