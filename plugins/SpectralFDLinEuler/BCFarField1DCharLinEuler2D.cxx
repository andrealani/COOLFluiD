#include "Framework/MethodStrategyProvider.hh"
#include "Common/ShouldNotBeHereException.hh"

#include "LinEuler/LinEuler2DVarSet.hh"
#include "LinEuler/LinEulerTerm.hh"

#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"
#include "SpectralFDLinEuler/BCFarField1DCharLinEuler2D.hh"

#include "Common/NotImplementedException.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCFarField1DCharLinEuler2D,SpectralFDMethodData,BCStateComputer,SpectralFDLinEulerModule >
    BCFarField1DCharLinEuler2DProvider("FarField1DCharLinEuler2D");

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
BCFarField1DCharLinEuler2D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_meanflow);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

BCFarField1DCharLinEuler2D::BCFarField1DCharLinEuler2D(const std::string& name) :
  BCStateComputer(name),
  socket_meanflow("meanflow"),
  m_linEulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()/*,
  m_tanVel()*/
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

BCFarField1DCharLinEuler2D::~BCFarField1DCharLinEuler2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCFarField1DCharLinEuler2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());
  cf_assert(nbrStates <= m_extraVars->size());
  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // normal
    const RealVector& normal = normals[iState];

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);
    cf_assert(intState.size() == 4);
    cf_assert(ghostState.size() == 4);
    cf_assert((*m_extraVars)[iState]->size() == 4);

     /// acquaintance of the PhysicalModel
     Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model =  m_linEulerVarSet->getModel();

     CFuint stateID = intState.getLocalID();

     RealVector meanflow_state = _model->getMeanFlowState(stateID);

    // set the extra variables
    m_linEulerVarSet->setExtraPhysicalVars(&meanflow_state);

    // set the physical data starting from the inner state
    m_linEulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // get mean flow data
    const CFreal rho0   = m_intSolPhysData[LinEulerTerm::rho0 ];
    const CFreal U0     = m_intSolPhysData[LinEulerTerm::U0   ];
    const CFreal V0     = m_intSolPhysData[LinEulerTerm::V0   ];
    const CFreal C0     = m_intSolPhysData[LinEulerTerm::c    ];
    const CFreal C0Sq   = C0*C0;
    const CFreal V0n    = U0*normal[XX] + V0*normal[YY];

    // check flow condition
    if (V0n < -C0) // supersonic inflow (normal points outward)
    {
      // for now, no incoming perturbations
      m_ghostSolPhysData[LinEulerTerm::rho] = -m_intSolPhysData[LinEulerTerm::rho];
      m_ghostSolPhysData[LinEulerTerm::u]   = -m_intSolPhysData[LinEulerTerm::u  ];
      m_ghostSolPhysData[LinEulerTerm::v]   = -m_intSolPhysData[LinEulerTerm::v  ];
      m_ghostSolPhysData[LinEulerTerm::p]   = -m_intSolPhysData[LinEulerTerm::p  ];
    }
    else if (V0n < 0.) // subsonic inflow
    {
      // pressure -> from internal solution
      const CFreal p = m_intSolPhysData[LinEulerTerm::p];

      // incoming characteristics are zero for now
      // density perturbation (from 2nd characteristic var (entropy) == 0)
      const CFreal rho = p/C0Sq; // zero entropy perturbation relation

      // tangential velocity component (from 3rd characteristic var == 0) (is forced to zero -> this causes reflections)
//       m_tanVel[XX] = 0.0;
//       m_tanVel[YY] = 0.0;

      // normal velocity component (from 4th characteristic var == 0)
      const CFreal vn = -p/(rho0*C0); // acoustic impedance relation

      // velocity components
      const CFreal u = vn*normal[XX]; // + m_tanVel[XX];
      const CFreal v = vn*normal[YY]; // + m_tanVel[YY];

      // set ghost solution values
      m_ghostSolPhysData[LinEulerTerm::rho] = 2.0*rho - m_intSolPhysData[LinEulerTerm::rho];
      m_ghostSolPhysData[LinEulerTerm::u  ] = 2.0*u   - m_intSolPhysData[LinEulerTerm::u  ];
      m_ghostSolPhysData[LinEulerTerm::v  ] = 2.0*v   - m_intSolPhysData[LinEulerTerm::v  ];
      m_ghostSolPhysData[LinEulerTerm::p  ] = p;
    }
    else if (V0n < C0) // subsonic outflow
    {
      // density and velocity --> from internal solution (2nd, 3rd and 4th characteristic from internal)
      const CFreal rho = m_intSolPhysData[LinEulerTerm::rho];
      const CFreal u   = m_intSolPhysData[LinEulerTerm::u  ];
      const CFreal v   = m_intSolPhysData[LinEulerTerm::v  ];

      // normal velocity component
      const CFreal vn = u*normal[XX] + v*normal[YY];

      // pressure (from 1st characteristic var == 0)
      const CFreal p = rho0*C0*vn;// acoustic impedance relation

      // set ghost solution values
      m_ghostSolPhysData[LinEulerTerm::rho] = rho;
      m_ghostSolPhysData[LinEulerTerm::u  ] = u  ;
      m_ghostSolPhysData[LinEulerTerm::v  ] = v  ;
      m_ghostSolPhysData[LinEulerTerm::p  ] = 2.0*p   - m_intSolPhysData[LinEulerTerm::p];
    }
    else if (V0n >= C0) // supersonic outflow
    {
      m_ghostSolPhysData[LinEulerTerm::rho] = m_intSolPhysData[LinEulerTerm::rho];
      m_ghostSolPhysData[LinEulerTerm::u  ] = m_intSolPhysData[LinEulerTerm::u  ];
      m_ghostSolPhysData[LinEulerTerm::v  ] = m_intSolPhysData[LinEulerTerm::v  ];
      m_ghostSolPhysData[LinEulerTerm::p  ] = m_intSolPhysData[LinEulerTerm::p  ];
    }
    else
    {
      throw Common::ShouldNotBeHereException(FromHere(),"BCFarField1DCharLinEuler2D::computeGhostStates(): could not determine flow condition");
    }

    // set the ghost state from its physical data
    m_linEulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCFarField1DCharLinEuler2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                               std::vector< std::vector< RealVector* > >& ghostGrads,
                                               const std::vector< RealVector >& normals,
                                               const std::vector< RealVector >& coords)
{
  throw Common::NotImplementedException(FromHere(),"BCFarField1DCharLinEuler2D::computeGhostGradients(): this should not be needed for linearized Euler!!!");
}

//////////////////////////////////////////////////////////////////////////////

void BCFarField1DCharLinEuler2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // extra variables required
  m_needsExtraVars = true;

  // get linearized Euler 2D varset
  m_linEulerVarSet = getMethodData().getUpdateVar().d_castTo<LinEuler2DVarSet>();
  if (m_linEulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException(FromHere(),"Update variable set is not LinEuler2DVarSet in BCFarField1DCharLinEuler2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_linEulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_linEulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // resize m_tanVel
//   m_tanVel.resize(2);

  //Get mean flow socket
  socket_meanflow.setParentNamespace( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
