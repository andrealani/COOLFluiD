#include "Framework/State.hh"
#include "Environment/ObjectProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/SetQuadLagrangeP1FluxReconstructionP1StateCoord.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SetQuadLagrangeP1FluxReconstructionP1StateCoord,
               Framework::SetElementStateCoord,
               FluxReconstructionModule>
SetQuadLagrangeP1FluxReconstructionP1StateCoord("QuadLagrangeP1FluxReconstructionP1");

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP1StateCoord::operator() (const vector<Framework::Node*>& nodes,
                                         vector<Framework::State*>& states)
{
  cf_assert(states.size() == 4);
  cf_assert(nodes.size() == 4);
  
//   FluxReconstructionElementData* frElemData = new QuadFluxReconstructionElementData(1);
// 
//   Common::SafePtr< std::vector< RealVector > > solPnts = frElemData->getSolPntsLocalCoords();
// 
//   for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
//   {
//     m_solPnts1D[iSol] = (*solPnts1D)[iSol];
//   }
//   
//   _tempCoord = 0.5*(*nodes[0] + *nodes[3]);
//   states[1]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//   _tempCoord = 0.5*(*nodes[0] + *nodes[1]);
//   states[3]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//   _tempCoord = 0.25*(*nodes[0] + *nodes[1] + *nodes[2] + *nodes[3]);
//   states[4]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//   _tempCoord = 0.5*(*nodes[2] + *nodes[3]);
//   states[5]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
//   _tempCoord = 0.5*(*nodes[1] + *nodes[2]);
//   states[7]->setSpaceCoordinates(new Framework::Node(_tempCoord,false));
// 
//   delete frElemData;
// 
  // assign nodes to the states
//  states[0]->setSpaceCoordinates(nodes[0]);
//  states[1]->setSpaceCoordinates(nodes[3]);
//  states[2]->setSpaceCoordinates(nodes[1]);
//  states[3]->setSpaceCoordinates(nodes[2]);
//   
}

//////////////////////////////////////////////////////////////////////////////

void SetQuadLagrangeP1FluxReconstructionP1StateCoord::update(const vector<Framework::Node*>& nodes,
                                    vector<Framework::State*>& states)
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace ShapeFunctions

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
