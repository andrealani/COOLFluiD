#include "Framework/MethodCommandProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"//3D
#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/ComputeErrorEuler.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider<
    ComputeErrorEuler,FluxReconstructionSolverData,FluxReconstructionNavierStokesModule >
  ComputeErrorEulerProvider("ComputeErrorEuler");

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorEuler::ComputeErrorEuler(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_states("states"),
  m_eulerVarSet(CFNULL),
  m_cellAvgSolCoefs(),
  m_cell(),
  m_cellStates(),
  m_tpIntegrator(),
  m_quadPntCoords(),
  m_quadCoefs(),
  m_cellBuilder(CFNULL)
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorEuler::~ComputeErrorEuler()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  ComputeErrorEuler::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::execute()
{
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal R = m_eulerVarSet->getModel()->getR();
  CFreal s;
  CFreal ht;
  
  // data from testcases
  if (true)
  {
    const CFreal TtBump = 0.00365795;
    const CFreal ptBump = 1.186212306;
    const CFreal MBump = 0.5;
    const CFreal Tt = TtBump/m_eulerVarSet->getModel()->getTempRef();
    const CFreal pt = ptBump/m_eulerVarSet->getModel()->getPressRef();
    const CFreal M = MBump;
  
    // get some data from the physical model
    
    const CFreal p = pt/pow(1.0+(gamma-1.0)/2.0*M*M, gamma/(gamma-1.0));
    const CFreal T = Tt/(1.0+(gamma-1.0)/2.0*M*M);
    const CFreal rho = p/(R*T);
    const CFreal a = pow(gamma*p/rho,0.5);
    s = p*pow(rho,-gamma);
    ht = gamma/(gamma-1.0)*R*T+a*a*M*M/2.0;
    //CFLog(NOTICE, "gamma:" << gamma << "\n");
    //CFLog(NOTICE, "R:" << R << "\n");
    //CFLog(NOTICE, "p:" << p << "\n");
    //CFLog(NOTICE, "rho:" << rho << "\n");
    //CFLog(NOTICE, "rhoU:" << M*a << "\n");
  }
  else
  {
    const CFreal MCyl = 0.38;
    const CFreal pCyl = (gamma-1.0)*(2.60108 - 0.5*1.0*0.449622064*0.449622064);
    const CFreal rho = 1.0;
    const CFreal T = pCyl/(rho*R);
    const CFreal a = pow(gamma*pCyl/rho,0.5);
    s = pCyl*pow(rho,-gamma);
    ht = gamma/(gamma-1.0)*R*T+a*a*MCyl*MCyl/2.0;
  }
  CFreal es = 0.0;
  CFreal eht = 0.0;
  CFreal esLi = 0.0;
  CFreal ehtLi = 0.0;
  CFreal esL1 = 0.0;
  CFreal ehtL1 = 0.0;
  CFreal totaleht = 0.0;
  CFreal totales = 0.0;
  CFuint nbStates = 0;
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;
  
  // loop over element types, for the moment there should only be one
  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoDataCell.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();
      if ((*m_cellStates)[0]->isParUpdatable())
      {
      
      nbStates += m_cellStates->size();

      vector< RealVector > nodeCoord;
      vector< Node* >* nodes = m_cell->getNodes();
      for (CFuint iNode = 0; iNode < 4; ++iNode)//3D
      {
	nodeCoord.push_back(*((*nodes)[iNode]->getData()));
      }

      vector< CFreal > quadPntWeights = m_tpIntegrator.getQuadPntsWheights(nodeCoord);
      const CFuint nbrQPnts = m_quadPntCoords.size();
      cf_assert(quadPntWeights.size() == nbrQPnts);
      
      vector< RealVector > statesInQuadPnts;

      for (CFuint iPnt = 0; iPnt < nbrQPnts; ++iPnt)
      {
	RealVector temp(4);//3D
	temp = 0.0;
	statesInQuadPnts.push_back(temp);
	if(m_cell->getID() == 1)
	{
	  //CFLog(VERBOSE,"Weights: " << quadPntWeights[iPnt] << "\n");
	}
      }

      // extrapolate the left and right states to the flx pnts
      for (CFuint iPnt = 0; iPnt < nbrQPnts; ++iPnt)
      {
        for (CFuint iSol = 0; iSol < m_cellStates->size(); ++iSol)
        {
          statesInQuadPnts[iPnt] += m_quadCoefs[iPnt][iSol]*(*((*m_cellStates)[iSol]));
        }
      }
      
      vector< State* > quadStates;

      for (CFuint iPnt = 0; iPnt < nbrQPnts; ++iPnt)
      {
        quadStates.push_back(new State(statesInQuadPnts[iPnt]));
      }

      for (CFuint iSol = 0; iSol < nbrQPnts; ++iSol)
      {
        // set the physical data starting from the inner state
        m_eulerVarSet->computePhysicalData((*((*m_cellStates)[iSol])),m_solPhysData);
	//m_eulerVarSet->computePhysicalData(*(quadStates[iSol]),m_solPhysData);

        const CFreal M2 = m_solPhysData[EulerTerm::V]/m_solPhysData[EulerTerm::A];
        const CFreal p2 = m_solPhysData[EulerTerm::P];
        const CFreal T2 = m_solPhysData[EulerTerm::T];
        const CFreal rho2 = m_solPhysData[EulerTerm::RHO];
        const CFreal V2 = m_solPhysData[EulerTerm::V];
        const CFreal s2 = p2*pow(rho2,-gamma);
        const CFreal ht2 = gamma/(gamma-1.0)*R*T2+V2*V2/2.0;
    
        // L2
        CFreal errorht = pow((ht - ht2)/ht,2);
        const CFreal errors = pow((s - s2)/s,2);
	
	if(m_cell->getID() == 1)
	{
	  //CFLog(NOTICE,"error ht in sol: " << errorht << "\n");
	}
	if (errorht > 0.00001)
	{
	  //CFLog(NOTICE,"high error: " << errorht << " in " << m_cell->getID() << "\n");
	}

        eht += quadPntWeights[iSol]*errorht;
        es += quadPntWeights[iSol]*errors;

        // L 1
        const CFreal errorhtL1 = fabs(ht - ht2)/ht;
        const CFreal errorsL1 = fabs(s - s2)/s;

        ehtL1 += quadPntWeights[iSol]*errorhtL1;
        esL1 += quadPntWeights[iSol]*errorsL1;
      }
      
      // loop over the states
      for (CFuint iState = 0; iState < m_cellStates->size(); ++iState)
      {
        // set the physical data starting from the inner state
        m_eulerVarSet->computePhysicalData((*((*m_cellStates)[iState])),m_solPhysData);

        const CFreal M2 = m_solPhysData[EulerTerm::V]/m_solPhysData[EulerTerm::A];
        const CFreal p2 = m_solPhysData[EulerTerm::P];
        const CFreal T2 = m_solPhysData[EulerTerm::T];
        const CFreal rho2 = m_solPhysData[EulerTerm::RHO];
        const CFreal V2 = m_solPhysData[EulerTerm::V];
        const CFreal s2 = p2*pow(rho2,-gamma);
        const CFreal ht2 = gamma/(gamma-1.0)*R*T2+V2*V2/2.0;
    
        const CFreal errorhtLi = fabs(ht - ht2)/ht;
        const CFreal errorsLi = fabs(s - s2)/s;
        totaleht += errorhtLi;
	totales += errorsLi;
	if(errorhtLi > ehtLi)
	{
	  ehtLi = errorhtLi;
	}
	if(errorsLi > esLi)
	{
	  esLi = errorsLi;
	}
	
	//(*((*m_cellStates)[iState]))[0] = errorhtLi;
	//(*((*m_cellStates)[iState]))[1] = errorsLi;
      }
      
      for (CFuint iPnt = 0; iPnt < nbrQPnts; ++iPnt)
      {
        deletePtr(quadStates[iPnt]);
      }
      }
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  // L2
  const CFreal errorEnthalpy = pow(eht,0.5);
  const CFreal errorEntropy = pow(es,1);
  
  //L1
  const CFreal errorEnthalpyL1 = ehtL1;
  const CFreal errorEntropyL1 = esL1;
  
  //Li
  const CFreal errorEnthalpyLi = ehtLi;
  const CFreal errorEntropyLi = esLi;
  
  // print out errors
  CFLog(NOTICE, "error ht L2: " << errorEnthalpy << "\n");
  CFLog(NOTICE, "error s L2: " << errorEntropy << "\n");
  CFLog(NOTICE, "error ht L1: " << errorEnthalpyL1 << "\n");
  CFLog(NOTICE, "error s L1: " << errorEntropyL1 << "\n");
  CFLog(NOTICE, "error ht Linf: " << errorEnthalpyLi << "\n");
  CFLog(NOTICE, "error s Linf: " << errorEntropyLi << "\n");
  CFLog(NOTICE, "av error ht: " << totaleht/nbStates << "\n");
  CFLog(NOTICE, "av error s: " << totales/nbStates << "\n");
  //CFLog(NOTICE, "ht: " << ht << "\n");
  //CFLog(NOTICE, "s: " << s << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorEuler::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  FluxReconstructionSolverCom::setup();
  
  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();//3D
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in ComputeErrorEuler!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData  );
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // create TensorProductGaussIntegrator
  TensorProductGaussIntegrator tempGI(DIM_2D,static_cast<CFPolyOrder::Type>((frLocalData[0]->getPolyOrder())*2-1));//3D
  m_tpIntegrator = tempGI;

  // get quadrature point coordinates and wheights
  m_quadPntCoords = m_tpIntegrator.getQuadPntsMappedCoords();
  
  // get solution point local coordinates
  m_quadCoefs = frLocalData[0]->getSolPolyValsAtNode(m_quadPntCoords);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

