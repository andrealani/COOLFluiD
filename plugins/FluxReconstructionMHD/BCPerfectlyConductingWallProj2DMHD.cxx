#include "Framework/MethodStrategyProvider.hh"

#include "MHD/MHD2DProjectionVarSet.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BCPerfectlyConductingWallProj2DMHD.hh"

#include "Common/NotImplementedException.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCPerfectlyConductingWallProj2DMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCPerfectlyConductingWallProj2DMHDProvider("PerfectlyConductingWallProj2DMHD");

//////////////////////////////////////////////////////////////////////////////

BCPerfectlyConductingWallProj2DMHD::BCPerfectlyConductingWallProj2DMHD(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

BCPerfectlyConductingWallProj2DMHD::~BCPerfectlyConductingWallProj2DMHD()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCPerfectlyConductingWallProj2DMHD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

void BCPerfectlyConductingWallProj2DMHD::computeGhostStates(const vector< State* >& intStates,
                                         vector< State* >& ghostStates,
                                         const std::vector< RealVector >& normals,
                                         const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 9);
    cf_assert(ghostState.size() == 9);
    CFLog(VERBOSE, "state: " << intState << "\n");
    // set the physical data starting from the inner state
    m_varSet->computePhysicalData(intState,m_intSolPhysData);
    CFLog(VERBOSE, "physData: " << m_intSolPhysData << "\n");

    const CFreal vn = m_intSolPhysData[MHDTerm::VX]*normal[XX] + 
      m_intSolPhysData[MHDTerm::VY]*normal[YY];// + m_intSolPhysData[MHDTerm::VZ]*normal[ZZ];
    
    const CFreal bn = m_intSolPhysData[MHDTerm::BX]*normal[XX] +
      m_intSolPhysData[MHDTerm::BY]*normal[YY];// + m_intSolPhysData[MHDTerm::BZ]*normal[ZZ]; 
    
    m_ghostSolPhysData[MHDTerm::RHO] = m_intSolPhysData[MHDTerm::RHO];
    m_ghostSolPhysData[MHDTerm::VX] = m_intSolPhysData[MHDTerm::VX] - 2.0*vn*normal[XX];
    m_ghostSolPhysData[MHDTerm::VY] = m_intSolPhysData[MHDTerm::VY] - 2.0*vn*normal[YY];
    m_ghostSolPhysData[MHDTerm::VZ] = m_intSolPhysData[MHDTerm::VZ];// - 2.0*vn*normal[ZZ];
    m_ghostSolPhysData[MHDTerm::BX] = m_intSolPhysData[MHDTerm::BX] - 2.0*bn*normal[XX];
    m_ghostSolPhysData[MHDTerm::BY] = m_intSolPhysData[MHDTerm::BY] - 2.0*bn*normal[YY];
    m_ghostSolPhysData[MHDTerm::BZ] = m_intSolPhysData[MHDTerm::BZ];// - 2.0*bn*normal[ZZ];
    m_ghostSolPhysData[MHDTerm::P] = m_intSolPhysData[MHDTerm::P];
    m_ghostSolPhysData[MHDTerm::V] = m_intSolPhysData[MHDTerm::V];
    m_ghostSolPhysData[MHDTerm::A] = m_intSolPhysData[MHDTerm::A];
    m_ghostSolPhysData[MHDTerm::B] = m_intSolPhysData[MHDTerm::B];
    m_ghostSolPhysData[MHDProjectionTerm::PHI] = m_intSolPhysData[MHDProjectionTerm::PHI];
  

    // set the ghost state from its physical data
    m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCPerfectlyConductingWallProj2DMHD::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                            std::vector< std::vector< RealVector* > >& ghostGrads,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  const CFuint nbrGradVars = intGrads[0].size();

  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    const RealVector& normal = normals[iState];

    for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
    {
      const RealVector& gradInt = *intGrads[iState][iGradVar];
      RealVector& gradGhost = *ghostGrads[iState][iGradVar];

      gradGhost = gradInt;

      // Apply sign flip to normal component for Dirichlet-enforced variables:
      if (iGradVar == MHDTerm::VY || iGradVar == MHDTerm::BY)
      {
        // Reflect the component normal to the wall
        CFreal gradDotN = gradInt[XX]*normal[XX] + gradInt[YY]*normal[YY];
        gradGhost[XX] -= 2.0 * gradDotN * normal[XX];
        gradGhost[YY] -= 2.0 * gradDotN * normal[YY];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCPerfectlyConductingWallProj2DMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get MHD 3D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo<MHD2DProjectionVarSet>();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in BCSuperOutletProj3DMHD!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // get face builder
  m_faceBuilder = getMethodData().getFaceBuilder();

  m_nbrFaceFlxPnts = frLocalData[0]->getFaceFlxPntsFaceLocalCoords()->size();
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  CFuint nbrFaceFlxPntsMax= (order+1)*(order+1);

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

