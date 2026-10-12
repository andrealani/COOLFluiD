#include <cmath>

#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Common/PE.hh"

#include "Environment/Factory.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/ManufacturedForcingTNEQ.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ManufacturedForcingTNEQ, FluxReconstructionSolverData, FluxReconstructionNEQModule>
ManufacturedForcingTNEQProvider("ManufacturedForcingTNEQ");

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("SpeciesNames","Species names in the order of the chemistry library (N2 O2 NO N O in any order).");
  options.addConfigOption< std::vector<CFreal> >("RhoBar","Bar partial densities [kg/m3], same order as SpeciesNames.");
  options.addConfigOption< CFreal >("UBar","Velocity scale ub [m/s] (also used for v).");
  options.addConfigOption< CFreal >("TBar","Temperature scale Tb [K].");
  options.addConfigOption< CFreal >("TvBar","Vibrational temperature scale Tvb [K].");
  options.addConfigOption< CFreal >("Epsilon","Amplitude e of the perturbations (0: uniform state).");
  options.addConfigOption< CFreal >("FDStep","Step h [m] of the eighth-order central differences of the flux.");
}

//////////////////////////////////////////////////////////////////////////////

ManufacturedForcingTNEQ::ManufacturedForcingTNEQ(const std::string& name) :
  TNEQSourceTerm(name),
  m_speciesNames(),
  m_rhoBar(),
  m_uBar(),
  m_TBar(),
  m_TvBar(),
  m_eps(),
  m_fdStep(),
  m_mms(),
  m_inputToUpdateVar(),
  m_prim(),
  m_inputState(CFNULL),
  m_fdState(CFNULL),
  m_msStates(),
  m_msCellStates(),
  m_pData(),
  m_unitVectors(),
  m_fluxPlus(),
  m_fluxMinus(),
  m_divF(),
  m_divF2(),
  m_srcAtMS(),
  m_forcing(),
  m_isComputed(),
  m_checksDone(false),
  m_maxForcing(),
  m_maxSource(),
  m_maxDivF(),
  m_maxDivFDiff()
{
  addConfigOptionsTo(this);

  // defaults: Freno et al. (2021), Sec. 4.2.2, species in the PLATO air5 order
  m_speciesNames = vector<std::string>();
  m_speciesNames.push_back("N");
  m_speciesNames.push_back("O");
  m_speciesNames.push_back("N2");
  m_speciesNames.push_back("NO");
  m_speciesNames.push_back("O2");
  setParameter("SpeciesNames",&m_speciesNames);

  m_rhoBar = vector<CFreal>();
  m_rhoBar.push_back(1.0e-4);
  m_rhoBar.push_back(1.0e-4);
  m_rhoBar.push_back(7.7e-3);
  m_rhoBar.push_back(1.0e-4);
  m_rhoBar.push_back(2.0e-3);
  setParameter("RhoBar",&m_rhoBar);

  m_uBar = 11484.10;
  setParameter("UBar",&m_uBar);

  m_TBar = 5000.;
  setParameter("TBar",&m_TBar);

  m_TvBar = 1000.;
  setParameter("TvBar",&m_TvBar);

  m_eps = 0.05;
  setParameter("Epsilon",&m_eps);

  m_fdStep = 1.0e-3;
  setParameter("FDStep",&m_fdStep);
}

//////////////////////////////////////////////////////////////////////////////

ManufacturedForcingTNEQ::~ManufacturedForcingTNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  TNEQSourceTerm::configure(args);

  // the forcing is constant: no numerical Jacobian, and the analytical one of the parent would add the TNEQ
  // Jacobian a second time
  if (m_addJacob || m_useAnaJacob)
  {
    CFLog(NOTICE, "ManufacturedForcingTNEQ: AddJacob and AnalyticalJacob ignored, the forcing has no Jacobian\n");
  }
  m_addJacob = false;
  m_useAnaJacob = false;

  m_mms.setup(m_speciesNames, m_rhoBar, m_uBar, m_TBar, m_TvBar, m_eps);
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::setup()
{
  CFAUTOTRACE;

  TNEQSourceTerm::setup();

  // checks: 2D, dimensional run, one equation per species, u, v, T and Tv
  SafePtr<PhysicalModelImpl> model = PhysicalModelStack::getActive()->getImplementor();
  if (m_dim != 2)
  {
    throw Common::BadValueException (FromHere(),"ManufacturedForcingTNEQ: 2D only");
  }
  if (model->isAdimensional() || std::abs(model->getRefLength() - 1.) > 0.)
  {
    throw Common::BadValueException (FromHere(),"ManufacturedForcingTNEQ: needs a dimensional run with reference length 1");
  }
  const CFuint nbSpecies = m_eulerVarSet->getModel()->getNbScalarVars(0);
  if (nbSpecies != m_mms.nbSpecies() || m_nbrEqs != nbSpecies + 4)
  {
    throw Common::BadValueException (FromHere(),"ManufacturedForcingTNEQ: needs 5 species and the equations rho_s, u, v, T, Tv");
  }

  // the names must match the library molar masses
  RealVector mm(nbSpecies);
  m_library->getMolarMasses(mm);
  m_mms.checkMolarMasses(mm);

  // the transformer from the manufactured primitive variables (RhoivtTv) to the update variables
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  const std::string provider = VarSetTransformer::getProviderName(physModel->getNameImplementor(), "RhoivtTv", updateVarStr);
  m_inputToUpdateVar = Environment::Factory<VarSetTransformer>::getInstance().getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);

  m_prim.resize(m_nbrEqs);
  m_inputState = new State(RealVector(m_nbrEqs));
  m_fdState = new State(RealVector(m_nbrEqs));
  m_msStates.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_msStates[iSol] = new State(RealVector(m_nbrEqs));
  }

  SafePtr<BaseTerm> convTerm = model->getConvectiveTerm();
  convTerm->resizePhysicalData(m_pData);

  m_unitVectors.resize(2);
  for (CFuint iDim = 0; iDim < 2; ++iDim)
  {
    m_unitVectors[iDim].resize(2);
    m_unitVectors[iDim] = 0.;
    m_unitVectors[iDim][iDim] = 1.;
  }

  m_fluxPlus.resize(m_nbrEqs);
  m_fluxMinus.resize(m_nbrEqs);
  m_divF.resize(m_nbrEqs);
  m_divF2.resize(m_nbrEqs);
  m_srcAtMS.resize(m_nbrSolPnts*m_nbrEqs);

  // forcing storage, one entry per state and equation
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
  const CFuint nbStates = rhs.size()/m_nbrEqs;
  m_forcing.assign(nbStates*m_nbrEqs, 0.);
  m_isComputed.assign(nbStates, false);
  m_checksDone = false;

  m_maxForcing.resize(m_nbrEqs);
  m_maxSource.resize(m_nbrEqs);
  m_maxDivF.resize(m_nbrEqs);
  m_maxDivFDiff.resize(m_nbrEqs);
  m_maxForcing = 0.;
  m_maxSource = 0.;
  m_maxDivF = 0.;
  m_maxDivFDiff = 0.;

  CFLog(INFO, "ManufacturedForcingTNEQ: Freno et al. (2021) Eq. (19), e = " << m_eps << ", ub = " << m_uBar
        << " m/s, Tb = " << m_TBar << " K, Tvb = " << m_TvBar << " K, FD step " << m_fdStep << " m\n");
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::unsetup()
{
  CFAUTOTRACE;

  deletePtr(m_inputState);
  deletePtr(m_fdState);
  for (CFuint iSol = 0; iSol < m_msStates.size(); ++iSol)
  {
    deletePtr(m_msStates[iSol]);
  }
  m_msStates.clear();

  TNEQSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::execute()
{
  TNEQSourceTerm::execute();

  if (m_checksDone) return;
  m_checksDone = true;

  // global maxima of the checks of the first pass
  const std::string nsp = getMethodData().getNamespace();
#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel())
  {
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    RealVector loc(m_nbrEqs);
    loc = m_maxForcing;  MPI_Allreduce(&loc[0], &m_maxForcing[0],  m_nbrEqs, MPI_DOUBLE, MPI_MAX, comm);
    loc = m_maxSource;   MPI_Allreduce(&loc[0], &m_maxSource[0],   m_nbrEqs, MPI_DOUBLE, MPI_MAX, comm);
    loc = m_maxDivF;     MPI_Allreduce(&loc[0], &m_maxDivF[0],     m_nbrEqs, MPI_DOUBLE, MPI_MAX, comm);
    loc = m_maxDivFDiff; MPI_Allreduce(&loc[0], &m_maxDivFDiff[0], m_nbrEqs, MPI_DOUBLE, MPI_MAX, comm);
  }
#endif
  if (PE::GetPE().GetRank(nsp) == 0)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal relFD = (m_maxDivF[iEq] > 0.) ? m_maxDivFDiff[iEq]/m_maxDivF[iEq] : 0.;
      CFLog(INFO, "ManufacturedForcingTNEQ: eq " << iEq << " max|f| " << m_maxForcing[iEq]
            << " max|S(U_MS)| " << m_maxSource[iEq] << " max|divF| " << m_maxDivF[iEq]
            << " FD check max|divF_h - divF_2h|/max|divF_h| " << relFD << "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::addSourceTerm(RealVector& resUpdates)
{
  const CFuint nbrSol = m_cellStates->size();

  // the forcing of a cell is computed the first time the cell is visited
  if (!m_isComputed[(*m_cellStates)[0]->getLocalID()])
  {
    computeForcingInCell();
  }

  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      resUpdates[m_nbrEqs*iSol + iEq] = m_forcing[m_nbrEqs*stateID + iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::computeForcingInCell()
{
  const CFuint nbrSol = m_cellStates->size();
  cf_assert(nbrSol <= m_msStates.size());

  // manufactured states at the solution points
  m_msCellStates.resize(nbrSol);
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    const Node& node = (*m_cellStates)[iSol]->getCoordinates();
    setManufacturedState(node[XX], node[YY], *m_msStates[iSol]);
    m_msCellStates[iSol] = m_msStates[iSol];
  }

  // TNEQ source at the manufactured states: the parent works on the states of m_cellStates, so it is pointed
  // to the manufactured states for this call
  std::vector<State*>* cellStates = m_cellStates;
  m_cellStates = &m_msCellStates;
  m_srcAtMS = 0.;
  TNEQSourceTerm::addSourceTerm(m_srcAtMS);
  m_cellStates = cellStates;

  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    const Node& node = (*m_cellStates)[iSol]->getCoordinates();
    computeDivFlux(node[XX], node[YY], m_fdStep, m_divF);

    // FD check with the double step, first pass only
    if (!m_checksDone)
    {
      computeDivFlux(node[XX], node[YY], 2.*m_fdStep, m_divF2);
    }

    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      const CFreal src = m_srcAtMS[m_nbrEqs*iSol + iEq];
      const CFreal f = m_divF[iEq] - src;
      m_forcing[m_nbrEqs*stateID + iEq] = f;

      if (!m_checksDone)
      {
        m_maxForcing[iEq]  = std::max(m_maxForcing[iEq], std::abs(f));
        m_maxSource[iEq]   = std::max(m_maxSource[iEq], std::abs(src));
        m_maxDivF[iEq]     = std::max(m_maxDivF[iEq], std::abs(m_divF[iEq]));
        m_maxDivFDiff[iEq] = std::max(m_maxDivFDiff[iEq], std::abs(m_divF[iEq] - m_divF2[iEq]));
      }
    }
    m_isComputed[stateID] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::setManufacturedState(const CFreal x, const CFreal y, State& state)
{
  m_mms.evaluate(x, y, m_prim);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    (*m_inputState)[iEq] = m_prim[iEq];
  }
  const State& upd = *m_inputToUpdateVar->transform(m_inputState);
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    state[iEq] = upd[iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::computeDivFlux(const CFreal x, const CFreal y, const CFreal h, RealVector& divF)
{
  // first derivative, eighth order: f' = sum_k c_k (f(x + k h) - f(x - k h))/h, k = 1..4
  const CFreal c[4] = {4./5., -1./5., 4./105., -1./280.};
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();

  divF = 0.;
  for (CFuint iDim = 0; iDim < 2; ++iDim)
  {
    const RealVector& e = m_unitVectors[iDim];
    for (CFuint k = 1; k <= 4; ++k)
    {
      const CFreal d = k*h;

      setManufacturedState(x + d*e[XX], y + d*e[YY], *m_fdState);
      updateVarSet->computePhysicalData(*m_fdState, m_pData);
      m_fluxPlus = updateVarSet->getFlux()(m_pData, e);

      setManufacturedState(x - d*e[XX], y - d*e[YY], *m_fdState);
      updateVarSet->computePhysicalData(*m_fdState, m_pData);
      m_fluxMinus = updateVarSet->getFlux()(m_pData, e);

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        divF[iEq] += c[k-1]*(m_fluxPlus[iEq] - m_fluxMinus[iEq])/h;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ManufacturedForcingTNEQ::getSToStateJacobian(const CFuint iState)
{
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    m_stateJacobian[iEq] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
