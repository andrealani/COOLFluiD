#include <cmath>
#include <fstream>
#include <iomanip>

#include "Common/CFLog.hh"
#include "Common/BadValueException.hh"
#include "Common/PE.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/ComputeErrorMMSTNEQ.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeErrorMMSTNEQ, FluxReconstructionSolverData, FluxReconstructionNEQModule>
ComputeErrorMMSTNEQProvider("ComputeErrorMMSTNEQ");

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMMSTNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("SpeciesNames","Species names in the order of the chemistry library (N2 O2 NO N O in any order).");
  options.addConfigOption< std::vector<CFreal> >("RhoBar","Bar partial densities [kg/m3], same order as SpeciesNames.");
  options.addConfigOption< CFreal >("UBar","Velocity scale ub [m/s] (also used for v).");
  options.addConfigOption< CFreal >("TBar","Temperature scale Tb [K].");
  options.addConfigOption< CFreal >("TvBar","Vibrational temperature scale Tvb [K].");
  options.addConfigOption< CFreal >("Epsilon","Amplitude e of the perturbations (0: uniform state).");
  options.addConfigOption< std::string >("OutputFile","File the errors are appended to, one line per residual evaluation.");
  options.addConfigOption< CFuint >("QuadPointsPerDir","Gauss-Legendre points per direction for the error integrals (0: P+3).");
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorMMSTNEQ::ComputeErrorMMSTNEQ(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_rhs("rhs"),
  m_speciesNames(),
  m_rhoBar(),
  m_uBar(),
  m_TBar(),
  m_TvBar(),
  m_eps(),
  m_outputFile(),
  m_nbQuadPerDir(),
  m_mms(),
  m_nbrEqs(),
  m_nbSpecies(),
  m_logVariables(false),
  m_quadLocalCoords(),
  m_quadWeights(),
  m_solPolyValsAtQuad(),
  m_solWeights(),
  m_nbCalls(0)
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

  m_outputFile = "mms_errors.dat";
  setParameter("OutputFile",&m_outputFile);

  m_nbQuadPerDir = 0;
  setParameter("QuadPointsPerDir",&m_nbQuadPerDir);
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorMMSTNEQ::~ComputeErrorMMSTNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > > ComputeErrorMMSTNEQ::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_rhs);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMMSTNEQ::gaussLegendre(const CFuint n, std::vector<CFreal>& x, std::vector<CFreal>& w)
{
  // Newton iteration on the Legendre polynomial P_n, p1 = P_n(z), p2 = P_{n-1}(z), dp = P_n'(z)
  const CFreal pi = 3.141592653589793;
  x.resize(n);
  w.resize(n);
  for (CFuint i = 0; i < n; ++i)
  {
    CFreal z = std::cos(pi*(i + 0.75)/(n + 0.5));
    CFreal dp = 1.;
    for (CFuint it = 0; it < 100; ++it)
    {
      CFreal p1 = 1., p2 = 0.;
      for (CFuint j = 1; j <= n; ++j)
      {
        const CFreal p3 = p2;
        p2 = p1;
        p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/j;
      }
      dp = n*(z*p1 - p2)/(z*z - 1.);
      const CFreal dz = p1/dp;
      z -= dz;
      if (std::abs(dz) < 1.e-16) break;
    }
    // derivative at the converged root for the weight
    CFreal p1 = 1., p2 = 0.;
    for (CFuint j = 1; j <= n; ++j)
    {
      const CFreal p3 = p2;
      p2 = p1;
      p1 = ((2.*j - 1.)*z*p2 - (j - 1.)*p3)/j;
    }
    dp = n*(z*p1 - p2)/(z*z - 1.);
    // the roots come in descending order
    x[n-1-i] = z;
    w[n-1-i] = 2./((1. - z*z)*dp*dp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMMSTNEQ::setup()
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::setup();

  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // checks: 2D, dimensional run
  SafePtr<PhysicalModelImpl> model = PhysicalModelStack::getActive()->getImplementor();
  if (PhysicalModelStack::getActive()->getDim() != 2)
  {
    throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: 2D only");
  }
  if (model->isAdimensional() || std::abs(model->getRefLength() - 1.) > 0.)
  {
    throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: needs a dimensional run with reference length 1");
  }

  // variable layout: rho_s (or ln rho_s), u, v, T (or ln T), Tv (or ln Tv)
  SafePtr< MultiScalarVarSet<Euler2DVarSet> > varSet =
    getMethodData().getUpdateVar().d_castTo< MultiScalarVarSet<Euler2DVarSet> >();
  m_nbSpecies = varSet->getModel()->getNbScalarVars(0);
  const vector<std::string>& names = getMethodData().getUpdateVar()->getVarNames();
  m_logVariables = (names[0] == "lnrho0");
  const bool linVariables = (names[0] == "rho0");
  if ((!m_logVariables && !linVariables) || m_nbrEqs != m_nbSpecies + 4 ||
      names[m_nbSpecies] != "u" || names[m_nbSpecies + 1] != "v" ||
      names[m_nbSpecies + 2] != (m_logVariables ? "lnT" : "T"))
  {
    throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: update variables must be RhoivtTv or LogRhoivLogTTv with one Tv");
  }

  m_mms.setup(m_speciesNames, m_rhoBar, m_uBar, m_TBar, m_TvBar, m_eps);
  if (m_nbSpecies != m_mms.nbSpecies())
  {
    throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: the mixture must have the 5 species of the solution");
  }
  RealVector mm(m_nbSpecies);
  model->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>()->getMolarMasses(mm);
  m_mms.checkMolarMasses(mm);

  // quadrature points: tensor product of n Gauss-Legendre points, n = P+3 by default
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  if (frLocalData.size() != 1 || frLocalData[0]->getShape() != CFGeoShape::QUAD)
  {
    throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: quadrilaterals only");
  }
  const CFuint order = static_cast<CFuint>(frLocalData[0]->getPolyOrder());
  const CFuint nQuad = (m_nbQuadPerDir > 0) ? m_nbQuadPerDir : order + 3;
  vector<CFreal> xq, wq;
  gaussLegendre(nQuad, xq, wq);
  m_quadLocalCoords.clear();
  m_quadWeights.clear();
  for (CFuint i = 0; i < nQuad; ++i)
  {
    for (CFuint j = 0; j < nQuad; ++j)
    {
      RealVector c(2);
      c[KSI] = xq[i];
      c[ETA] = xq[j];
      m_quadLocalCoords.push_back(c);
      m_quadWeights.push_back(wq[i]*wq[j]);
    }
  }
  m_solPolyValsAtQuad = frLocalData[0]->getSolPolyValsAtNode(m_quadLocalCoords);

  // solution point weights: the solution points must be the Gauss-Legendre points of order P
  vector<CFreal> xs, ws;
  gaussLegendre(order + 1, xs, ws);
  SafePtr< vector< RealVector > > solCoords = frLocalData[0]->getSolPntsLocalCoords();
  m_solWeights.resize(solCoords->size());
  for (CFuint iSol = 0; iSol < solCoords->size(); ++iSol)
  {
    CFreal w = 1.;
    for (CFuint iDim = 0; iDim < 2; ++iDim)
    {
      CFint idx = -1;
      for (CFuint k = 0; k < xs.size(); ++k)
      {
        if (std::abs((*solCoords)[iSol][iDim] - xs[k]) < 1.e-12) idx = k;
      }
      if (idx < 0)
      {
        throw Common::BadValueException (FromHere(),"ComputeErrorMMSTNEQ: the solution points are not Gauss-Legendre points");
      }
      w *= ws[idx];
    }
    m_solWeights[iSol] = w;
  }

  m_nbCalls = 0;

  CFLog(INFO, "ComputeErrorMMSTNEQ: " << nQuad << "x" << nQuad << " Gauss-Legendre points per cell, "
        << (m_logVariables ? "LogRhoivLogTTv" : "RhoivtTv") << " update variables, output " << m_outputFile << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMMSTNEQ::execute()
{
  CFAUTOTRACE;

  const CFuint nbFields = m_nbrEqs;
  const CFuint TID = m_nbSpecies + 2;

  // accumulated per field: integral |e|, integral e^2, max |e|; per equation: max |R|, integral |R|
  vector<CFreal> errL1(nbFields, 0.), errL2(nbFields, 0.), errLinf(nbFields, 0.);
  // per field: x and y of the quadrature point where |e| is largest
  vector<CFreal> xLinf(nbFields, 0.), yLinf(nbFields, 0.);
  vector<CFreal> resLinf(m_nbrEqs, 0.), resL1(m_nbrEqs, 0.);
  CFreal area = 0.;
  CFuint nbCells = 0, nbStates = 0;

  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  SafePtr< vector< RealVector > > solCoords = frLocalData[0]->getSolPntsLocalCoords();

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  RealVector stateQ(m_nbrEqs);
  RealVector prim(m_nbrEqs);
  RealVector exact(m_nbrEqs);

  const CFuint startIdx = (*elemType)[0].getStartIdx();
  const CFuint endIdx   = (*elemType)[0].getEndIdx();
  for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
  {
    geoData.idx = elemIdx;
    GeometricEntity* cell = geoBuilder->buildGE();
    vector<State*>* cellStates = cell->getStates();
    const CFuint nbrSol = cellStates->size();

    if ((*cellStates)[0]->isParUpdatable())
    {
      ++nbCells;
      nbStates += nbrSol;

      // errors at the quadrature points
      const std::valarray<CFreal> jacobQuad = cell->computeGeometricShapeFunctionJacobianDeterminant(m_quadLocalCoords);
      for (CFuint iQ = 0; iQ < m_quadLocalCoords.size(); ++iQ)
      {
        stateQ = 0.;
        for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
        {
          stateQ += m_solPolyValsAtQuad[iQ][iSol]*(*(*cellStates)[iSol]);
        }
        for (CFuint i = 0; i < m_nbrEqs; ++i)
        {
          const bool isLog = m_logVariables && (i < m_nbSpecies || i >= TID);
          prim[i] = isLog ? std::exp(stateQ[i]) : stateQ[i];
        }

        const RealVector coord = cell->computeCoordFromMappedCoord(m_quadLocalCoords[iQ]);
        m_mms.evaluate(coord[XX], coord[YY], exact);

        const CFreal wJ = m_quadWeights[iQ]*jacobQuad[iQ];
        area += wJ;
        for (CFuint i = 0; i < nbFields; ++i)
        {
          const CFreal e = std::abs(prim[i] - exact[i])/m_mms.barValue(i);
          errL1[i] += wJ*e;
          errL2[i] += wJ*e*e;
          if (e > errLinf[i])
          {
            errLinf[i] = e;
            xLinf[i] = coord[XX];
            yLinf[i] = coord[YY];
          }
        }
      }

      // residual at the solution points
      const std::valarray<CFreal> jacobSol = cell->computeGeometricShapeFunctionJacobianDeterminant(*solCoords);
      for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
      {
        const CFuint stateID = (*cellStates)[iSol]->getLocalID();
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          const CFreal r = rhs[stateID*m_nbrEqs + iEq];
          resLinf[iEq] = std::max(resLinf[iEq], std::abs(r)/jacobSol[iSol]);
          resL1[iEq] += m_solWeights[iSol]*std::abs(r);
        }
      }
    }

    geoBuilder->releaseGE();
  }

  // reductions over the ranks
  const std::string nsp = getMethodData().getNamespace();
#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel())
  {
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    vector<CFreal> loc;
    loc = errL1;   MPI_Allreduce(&loc[0], &errL1[0],   nbFields, MPI_DOUBLE, MPI_SUM, comm);
    loc = errL2;   MPI_Allreduce(&loc[0], &errL2[0],   nbFields, MPI_DOUBLE, MPI_SUM, comm);
    // location of the maximum: gather (|e|, x, y) of every rank and keep the largest |e|
    const CFuint nbRanks = PE::GetPE().GetProcessorCount(nsp);
    vector<CFreal> locMax(3*nbFields), allMax(3*nbFields*nbRanks);
    for (CFuint i = 0; i < nbFields; ++i)
    {
      locMax[3*i] = errLinf[i]; locMax[3*i+1] = xLinf[i]; locMax[3*i+2] = yLinf[i];
    }
    MPI_Allgather(&locMax[0], 3*nbFields, MPI_DOUBLE, &allMax[0], 3*nbFields, MPI_DOUBLE, comm);
    for (CFuint r = 0; r < nbRanks; ++r)
    {
      for (CFuint i = 0; i < nbFields; ++i)
      {
        const CFuint k = 3*(r*nbFields + i);
        if (allMax[k] > errLinf[i])
        {
          errLinf[i] = allMax[k]; xLinf[i] = allMax[k+1]; yLinf[i] = allMax[k+2];
        }
      }
    }
    loc = resLinf; MPI_Allreduce(&loc[0], &resLinf[0], m_nbrEqs, MPI_DOUBLE, MPI_MAX, comm);
    loc = resL1;   MPI_Allreduce(&loc[0], &resL1[0],   m_nbrEqs, MPI_DOUBLE, MPI_SUM, comm);
    CFreal areaL = area;
    MPI_Allreduce(&areaL, &area, 1, MPI_DOUBLE, MPI_SUM, comm);
    CFuint nL = nbCells;
    MPI_Allreduce(&nL, &nbCells, 1, MPI_UNSIGNED, MPI_SUM, comm);
    nL = nbStates;
    MPI_Allreduce(&nL, &nbStates, 1, MPI_UNSIGNED, MPI_SUM, comm);
  }
#endif

  if (PE::GetPE().GetRank(nsp) == 0)
  {
    const vector<std::string> fields = m_mms.fieldNames();
    const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

    std::ofstream out;
    if (m_nbCalls == 0)
    {
      out.open(m_outputFile.c_str());
      out << "# Freno et al. (2021) MMS errors, e = (solution - exact)/bar value; residual R per equation\n";
      out << "# call iter nbCells nbStates area";
      for (CFuint i = 0; i < nbFields; ++i)
      {
        out << " L1_" << fields[i] << " L2_" << fields[i] << " Linf_" << fields[i];
      }
      for (CFuint i = 0; i < nbFields; ++i)
      {
        out << " xLinf_" << fields[i] << " yLinf_" << fields[i];
      }
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        out << " resLinf_" << iEq << " resL1_" << iEq;
      }
      out << "\n";
    }
    else
    {
      out.open(m_outputFile.c_str(), std::ios::app);
    }
    out << std::setprecision(10) << std::scientific;
    out << m_nbCalls << " " << iter << " " << nbCells << " " << nbStates << " " << area;
    for (CFuint i = 0; i < nbFields; ++i)
    {
      out << " " << errL1[i] << " " << std::sqrt(errL2[i]) << " " << errLinf[i];
    }
    for (CFuint i = 0; i < nbFields; ++i)
    {
      out << " " << xLinf[i] << " " << yLinf[i];
    }
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      out << " " << resLinf[iEq] << " " << resL1[iEq];
    }
    out << "\n";
    out.close();

    CFLog(INFO, "ComputeErrorMMSTNEQ: call " << m_nbCalls << " Linf(" << fields[0] << ") " << errLinf[0]
          << " Linf(T) " << errLinf[TID] << " area " << area << "\n");
  }

  ++m_nbCalls;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
