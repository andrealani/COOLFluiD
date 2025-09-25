#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHDTerm.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/ComputeErrorMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/QuadFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/HexaFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TriagFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"

#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodCommandProvider<
    ComputeErrorMHD,FluxReconstructionSolverData,FluxReconstructionMHDModule >
  ComputeErrorMHDProvider("ComputeErrorMHD");

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<CFuint>
    ("QuadratureOrder", "Order of the quadrature rule used for error integration (default: 3)");
  
  options.addConfigOption<std::string>
    ("RefSolutionType", "Type of reference solution: 'Manufactured' (default), 'AlfvenWave', or 'MHDVortex'");
  
  options.addConfigOption<CFuint>
    ("ShowRate", "Show rate of the error information (default: 1)");
  
  options.addConfigOption<CFuint>
    ("MonitoredVar", "ID of the variable to monitor for error computation (default: 0 - rho)");
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorMHD::ComputeErrorMHD(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_states("states"),
  m_MHDVarSet(CFNULL),
  m_cellAvgSolCoefs(),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_tpIntegrator(),
  m_quadPntCoords(),
  m_quadSolPntsLocalCoords(),
  m_quadCellAvgSolCoefs(),
  m_solPolyValsAtQuadPnts(),
  m_nbrQuadSolPnts(),
  m_nbrEqs(),
  m_solPntsLocalCoords(),
  m_frLocalDataQuad(CFNULL),
  m_error_rho(),
  m_error_rho_Li(),
  m_error_rho_L1(),
  m_error_B(),
  m_total_error_rho(),
  m_total_error_rho_Li(),
  m_total_avg_e_rho(),
  m_total_error_rho_L1(),
  m_total_error_B(),
  m_quadCoefs(),
  m_quadOrder(),
  m_refSolutionType(MANUFACTURED),
  m_refSolutionTypeStr("Manufactured"),
  m_showRate(),
  m_monitoredVar()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_quadOrder = 3;
  setParameter("QuadratureOrder", &m_quadOrder);
  
  m_refSolutionTypeStr = "Manufactured";
  setParameter("RefSolutionType", &m_refSolutionTypeStr);
  
  m_showRate = 1;
  setParameter("ShowRate", &m_showRate);
  
  m_monitoredVar = 0;
  setParameter("MonitoredVar", &m_monitoredVar);
}

//////////////////////////////////////////////////////////////////////////////

ComputeErrorMHD::~ComputeErrorMHD()
{
  CFAUTOTRACE;
  
  // Clean up quadrature element data
  if (m_frLocalDataQuad != CFNULL)
  {
    deletePtr(m_frLocalDataQuad);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  ComputeErrorMHD::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_states);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::execute()
{
  CFTRACEBEGIN

  CFreal e_rho = 0.0;
  CFreal e_rhoLi = 0.0;
  CFreal e_rhoL1 = 0.0;
  CFreal e_B = 0.0;
  CFreal e_rhoLicell = 0.0;
  CFreal total_e_rho = 0.0;
  CFuint nbStates = 0;
  
  m_error_rho = 0.;
  m_error_rho_Li = 0.;
  m_error_rho_L1 = 0.;
  m_error_B = 0.;
  m_total_error_rho = 0.;
  m_total_error_rho_Li = 0.;
  m_total_avg_e_rho = 0.;
  m_total_error_rho_L1 = 0.;
  m_total_error_B = 0.;
  m_total_nbStates = 0;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoDataCell = m_cellBuilder->getDataGE();
  geoDataCell.trs = cells;

  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  
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

      // Count quadrature solution points for consistent counting across orders
      nbStates += m_nbrQuadSolPnts;

      // Compute Jacobian determinant at quadrature solution points for fixed integration
      const std::valarray<CFreal> jacobDetQuad =
                        m_cell->computeGeometricShapeFunctionJacobianDeterminant(*m_quadSolPntsLocalCoords);

      // Loop over quadrature solution points for fixed integration
      // This provides consistent integration accuracy regardless of simulation polynomial order
      for (CFuint iQuad = 0; iQuad < m_nbrQuadSolPnts; ++iQuad)
      {
        // Extrapolate solution from simulation solution points to quadrature point
        RealVector stateAtQuad(m_nbrEqs);
        stateAtQuad = 0.0;
        
        for (CFuint iSol = 0; iSol < (*m_cellStates).size(); ++iSol)
        {
          stateAtQuad += m_solPolyValsAtQuadPnts[iQuad][iSol] * (*((*m_cellStates)[iSol]));
        }

        // Get physical coordinates at quadrature point using cell geometric mapping
        const RealVector& quadLocalCoord = (*m_quadSolPntsLocalCoords)[iQuad];
        RealVector coordAtQuad = m_cell->computeCoordFromMappedCoord(quadLocalCoord);
        
        const CFreal x = coordAtQuad[XX];
        const CFreal y = coordAtQuad[YY];
        const CFreal z = (dim==DIM_3D) ? coordAtQuad[ZZ] : 0.0;

        // Compute analytical solution at quadrature point
        RealVector analyticalSol(m_nbrEqs);
        computeAnalyticalSolution(x, y, z, analyticalSol);

        const CFreal p = stateAtQuad[7];
        const CFreal rho = stateAtQuad[0];
        const CFreal u = stateAtQuad[1];
        const CFreal v = stateAtQuad[2];
        const CFreal w = stateAtQuad[3];

        const CFreal Bx = stateAtQuad[4];
        const CFreal By = stateAtQuad[5];
        const CFreal Bz = stateAtQuad[6];

        // L2 norm computation using FIXED quadrature integration
        CFreal errorho = pow((analyticalSol[m_monitoredVar] - stateAtQuad[m_monitoredVar]),2); // /analyticalSol[m_monitoredVar]
        e_rho += (*m_quadCellAvgSolCoefs)[iQuad]*errorho*jacobDetQuad[iQuad];

        // L∞ norm computation - pointwise maximum (no integration)
        const CFreal error_rhoLi = fabs(stateAtQuad[m_monitoredVar] - analyticalSol[m_monitoredVar]);
        e_rhoLi = std::max(e_rhoLi, error_rhoLi);

        // L1 norm computation
        e_rhoL1 += (*m_quadCellAvgSolCoefs)[iQuad] * fabs(stateAtQuad[m_monitoredVar] - analyticalSol[m_monitoredVar]) * jacobDetQuad[iQuad];

        // Special magnetic field error for MHD vortex case
        if (m_refSolutionType == MHD_VORTEX) {
          const CFreal Bx_error = fabs(stateAtQuad[4] - analyticalSol[4]);  // Bx
          const CFreal By_error = fabs(stateAtQuad[5] - analyticalSol[5]);  // By
          e_B += (*m_quadCellAvgSolCoefs)[iQuad] * (Bx_error + By_error) * jacobDetQuad[iQuad];
        }

        // Track total error for averaging
        total_e_rho += error_rhoLi;
      }
      // OLD cell-level max (was taking max of integrated values):
      // e_rhoLicell = std::max(e_rhoLicell,e_rhoLi);
      }
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  //L2
  m_error_rho = pow(e_rho,1.);
  
  //Li - CORRECTED: use the actual maximum error across all points
  m_error_rho_Li = e_rhoLi;
  
  //L1 
  m_error_rho_L1 = e_rhoL1;

  // Magnetic field error (normalized by A=20 for MHD vortex case)
  if (m_refSolutionType == MHD_VORTEX) {
    m_error_B = e_B / 20.0;  // Divide by A=20 BEFORE MPI reduction
  } else {
    m_error_B = 0.0;
  }
  
  // OLD INCORRECT Li computation (was using cell-level max of integrated values):
  // m_error_rho_Li = e_rhoLicell;

  //Av Err
  //m_avg_e_rho = total_e_rho ;
  
  // print out errors
  const std::string nsp = this->getMethodData().getNamespace();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  #ifdef CF_HAVE_MPI
      MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
      PE::GetPE().setBarrier(nsp);
      const CFuint count = 1;
      MPI_Allreduce(&m_error_rho, &m_total_error_rho, count, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(&m_error_rho_Li, &m_total_error_rho_Li, count, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(&m_error_rho_L1, &m_total_error_rho_L1, count, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(&m_error_B, &m_total_error_B, count, MPI_DOUBLE, MPI_SUM, comm);
      MPI_Allreduce(&nbStates, &m_total_nbStates, count, MPI_UNSIGNED, MPI_SUM, comm);
  #else
      // No MPI - totals equal local values
      m_total_error_rho = m_error_rho;
      m_total_error_rho_Li = m_error_rho_Li;
      m_total_error_rho_L1 = m_error_rho_L1;
      m_total_error_B = m_error_B;
      m_total_nbStates = nbStates;
  #endif
    
  if (PE::GetPE().GetRank(nsp) == 0 && iter % m_showRate == 0) 
  {
    // print errors
    CFLog(NOTICE, "error var[" << m_monitoredVar << "] L2: " << pow(m_total_error_rho,0.5) << "\n");
    CFLog(NOTICE, "error var[" << m_monitoredVar << "] Linf: " << m_total_error_rho_Li << "\n");
    CFLog(NOTICE, "error var[" << m_monitoredVar << "] L1: " << m_total_error_rho_L1 << "\n");
    
    // Special output for MHD vortex case
    if (m_refSolutionType == MHD_VORTEX) {
      CFLog(NOTICE, "error B: " << m_total_error_B << "\n");
    }
  }

  PE::GetPE().setBarrier(nsp);
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  FluxReconstructionSolverCom::setup();
  
  // Process reference solution type config option
  if (m_refSolutionTypeStr == "Manufactured") {
    m_refSolutionType = MANUFACTURED;
  } else if (m_refSolutionTypeStr == "AlfvenWave") {
    m_refSolutionType = ALFVEN_WAVE;
  } else if (m_refSolutionTypeStr == "MHDVortex") {
    m_refSolutionType = MHD_VORTEX;
  } else {
    throw Common::BadValueException (FromHere(),"RefSolutionType should be 'Manufactured', 'AlfvenWave', or 'MHDVortex'");
  }
  
  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get MHD 2D varset
  //m_MHDVarSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();//3D
  //if (m_MHDVarSet.isNull())
  //{
  //  throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHDProjectionVarSet in ComputeErrorMHD!");
  //}

  // resize the physical data for internal and ghost solution points
  //m_MHDVarSet->getModel()->resizePhysicalData(m_solPhysData  );
  
  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  
  // get element shape
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  
  // Create configurable order FluxReconstructionElementData for fixed integration rule
  // This ensures consistent integration accuracy across all polynomial orders
  // for convergence studies, avoiding order-dependent integration errors
  const CFPolyOrder::Type quadOrder = static_cast<CFPolyOrder::Type>(m_quadOrder);
  
  switch(elemShape)
  {
    case CFGeoShape::LINE:
    {
      throw Common::NotImplementedException (FromHere(),"FR has not been implemented for 1D");
    } break;
    case CFGeoShape::TRIAG:
    {
      m_frLocalDataQuad = new TriagFluxReconstructionElementData(quadOrder);
    } break;
    case CFGeoShape::QUAD:
    {
      m_frLocalDataQuad = new QuadFluxReconstructionElementData(quadOrder);
    } break;
    case CFGeoShape::TETRA:
    {
      m_frLocalDataQuad = new TetraFluxReconstructionElementData(quadOrder);      
    } break;
    case CFGeoShape::HEXA:
    {
      m_frLocalDataQuad = new HexaFluxReconstructionElementData(quadOrder);
    } break;
    case CFGeoShape::PRISM:
    {
      m_frLocalDataQuad = new PrismFluxReconstructionElementData(quadOrder);      
    } break;      
    default:
    {
      throw Common::NotImplementedException (FromHere(),"FR method not implemented for elements of type "
                                    + StringOps::to_str(elemShape) + ".");
    }
  }
  
  // Get quadrature solution points and cell averaging coefficients for fixed integration
  m_quadSolPntsLocalCoords = m_frLocalDataQuad->getSolPntsLocalCoords();
  m_quadCellAvgSolCoefs = m_frLocalDataQuad->getCellAvgSolCoefs();
  m_nbrQuadSolPnts = m_quadSolPntsLocalCoords->size();
  
  // Get polynomial values for extrapolation from simulation sol pnts to quadrature sol pnts
  m_solPolyValsAtQuadPnts = frLocalData[0]->getSolPolyValsAtNode(*m_quadSolPntsLocalCoords);

  // get the cell averaging coefficients for current simulation order (for comparison)
  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();

  // get solution point local coordinates for current simulation order
  m_solPntsLocalCoords = frLocalData[0]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::computeAnalyticalSolution(const CFreal x, const CFreal y, const CFreal z, 
                                                 RealVector& analyticalSol)
{
  switch(m_refSolutionType)
  {
    case MANUFACTURED:
      computeManufacturedSolution(x, y, z, analyticalSol);
      break;
    case ALFVEN_WAVE:
      computeAlfvenWaveSolution(x, y, z, analyticalSol);
      break;
    case MHD_VORTEX:
      computeMHDVortexSolution(x, y, z, analyticalSol);
      break;
    default:
      throw Common::ShouldNotBeHereException (FromHere(),"Unknown reference solution type");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::computeManufacturedSolution(const CFreal x, const CFreal y, const CFreal z,
                                                   RealVector& analyticalSol)
{
  // Ensure output vector has correct size
  analyticalSol.resize(m_nbrEqs);
  
  const CFreal r = sqrt(x*x + y*y + z*z);
  
  // MHD Manufactured solution
  analyticalSol[0] = pow(r, -5./2.);        // rho
  analyticalSol[1] = x / sqrt(r);           // u
  analyticalSol[2] = y / sqrt(r);           // v
  analyticalSol[3] = z / sqrt(r);           // w
  analyticalSol[4] = x / (r*r*r);           // Bx
  analyticalSol[5] = y / (r*r*r);           // By
  analyticalSol[6] = z / (r*r*r) + 0.017;   // Bz
  analyticalSol[7] = pow(r, -5./2.);        // p
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::computeAlfvenWaveSolution(const CFreal x, const CFreal y, const CFreal z,
                                                 RealVector& analyticalSol)
{
  // Ensure output vector has correct size
  analyticalSol.resize(m_nbrEqs);
  
  /*analyticalSol[0] = 1.0;  // rho
  analyticalSol[1] = -0.1*0.89415423683*sin(2*M_PI*(x*0.44775908783+y*0.89415423683));  // u
  analyticalSol[2] = 0.1*0.44775908783*sin(2*M_PI*(x*0.44775908783+y*0.89415423683));  // v
  analyticalSol[3] = 0.1*cos(2*M_PI*(x*0.44775908783+y*0.89415423683));  // w
  analyticalSol[4] = 0.44775908783-0.1*0.89415423683*sin(2*M_PI*(x*0.44775908783+y*0.89415423683));  // Bx
  analyticalSol[5] = 0.89415423683+0.1*0.44775908783*sin(2*M_PI*(x*0.44775908783+y*0.89415423683));  // By
  analyticalSol[6] = 0.1*cos(2*M_PI*(x*0.44775908783+y*0.89415423683));  // Bz
  analyticalSol[7] = 0.1;  // p
*/

  analyticalSol[0] = 1.0;  // rho
  analyticalSol[1] = -0.1*0.5*sin(2*M_PI*(x*0.8660254+y*0.5));  // u
  analyticalSol[2] = 0.1*0.8660254*sin(2*M_PI*(x*0.8660254+y*0.5));  // v
  analyticalSol[3] = 0.1*cos(2*M_PI*(x*0.8660254+y*0.5));  // w
  analyticalSol[4] = 0.8660254-0.1*0.5*sin(2*M_PI*(x*0.8660254+y*0.5));  // Bx
  analyticalSol[5] = 0.5+0.1*0.8660254*sin(2*M_PI*(x*0.8660254+y*0.5));  // By
  analyticalSol[6] = 0.1*cos(2*M_PI*(x*0.8660254+y*0.5));  // Bz
  analyticalSol[7] = 0.1;  // p

  /*
  analyticalSol[0] = 1.0;  // rho
  analyticalSol[1] = 0.;  // u
  analyticalSol[2] = 0.1*sin(2*M_PI*x);  // v
  analyticalSol[3] = 0.1*cos(2*M_PI*x);  // w
  analyticalSol[4] = 1.;  // Bx
  analyticalSol[5] = 0.1*sin(2*M_PI*x);  // By
  analyticalSol[6] = 0.1*cos(2*M_PI*x);  // Bz
  analyticalSol[7] = 0.1;  // p
  */
}

//////////////////////////////////////////////////////////////////////////////

void ComputeErrorMHD::computeMHDVortexSolution(const CFreal x, const CFreal y, const CFreal z,
                                                RealVector& analyticalSol)
{
  // Ensure output vector has correct size
  analyticalSol.resize(m_nbrEqs);
  
  // Translation velocity [1,1] 
  // For time t, the translated coordinates are: x_trans = x - t, y_trans = y - t
  // Since this is called during simulation, we need to get the current time
  const CFreal t = SubSystemStatusStack::getActive()->getCurrentTime();
  const CFreal x_trans = x - t;  // Translate by velocity [1,1]
  const CFreal y_trans = y - t;
  
  // Compute r for the translated coordinates
  const CFreal r = sqrt(x_trans*x_trans + y_trans*y_trans);
  
  // Constants from the initial solution
  const CFreal vortex_strength = 5.38948938512;
  const CFreal sqrt2 = sqrt(2.0);
  const CFreal pi = M_PI;
  
  // Exact solution (translation of initial conditions)
  analyticalSol[0] = 1.0;  // rho
  analyticalSol[1] = 1.0 - y_trans * (sqrt2 * vortex_strength) / (2.0 * pi) * exp(1.0 - r*r);  // u
  analyticalSol[2] = 1.0 + x_trans * (sqrt2 * vortex_strength) / (2.0 * pi) * exp(1.0 - r*r);  // v
  analyticalSol[3] = 0.0;  // w
  analyticalSol[4] = -y_trans * vortex_strength / (2.0 * pi) * exp(1.0 - r*r);  // Bx
  analyticalSol[5] = x_trans * vortex_strength / (2.0 * pi) * exp(1.0 - r*r);   // By
  analyticalSol[6] = 0.0;  // Bz
  
  // Pressure term
  const CFreal pressure_term = (vortex_strength*vortex_strength * (1.0 - r*r) - sqrt2*sqrt2 * vortex_strength*vortex_strength) / (8.0 * pi*pi);
  analyticalSol[7] = 1.0 + pressure_term * exp(1.0 - r*r);  // p
  
  analyticalSol[8] = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

