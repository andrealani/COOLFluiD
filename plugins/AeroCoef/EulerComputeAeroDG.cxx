#include "AeroCoef/AeroCoefDG.hh"
#include "AeroCoef/EulerComputeAeroDG.hh"

#include "Common/PE.hh"

#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
// #include "Framework/ParsingErrorException.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/EulerVarSet.hh"

#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "DiscontGalerkin/DiscontGalerkinSolver.hh"
#include "DiscontGalerkin/DiscontGalerkinSolverData.hh"

#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::DiscontGalerkin;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EulerComputeAeroDG, DataProcessingData, AeroCoefDGModule> eulerComputeAeroDG("EulerComputeAeroDG");

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Definition of the function defining the angle between the flow and the airfoil.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coeficients.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFreal >("uInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("rhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("pInf","Pressure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

EulerComputeAeroDG::EulerComputeAeroDG(const std::string& name) :
  DataProcessingCom(name),
  m_faceTrsGeoBuilder(),
  m_sockets(),
  m_updateVarSet(),
  m_vars("t"),
  m_eval(0.0,1),
  m_lift(),
  m_drag(),
  m_alphadeg(0.),
  m_alpharad(0.)
{
  m_fileWall = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFileWall = "Wall";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_saveRateWall = 1;
  setParameter("SaveRateWall",&m_saveRateWall);

  m_function = std::string("0.");
  setParameter("Alpha",&m_function);

  m_nameOutputFileAero = "AeroCoef.plt";
  setParameter("OutputFileAero",&m_nameOutputFileAero);

  m_saveRateAero = 1;
  setParameter("SaveRateAero",&m_saveRateAero);

  m_uInf = 0.;
  setParameter("uInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("rhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("pInf",&m_pInf);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

EulerComputeAeroDG::~EulerComputeAeroDG()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
EulerComputeAeroDG::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
EulerComputeAeroDG::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::setup()
{
 m_lift = 0.;
 m_drag = 0.;

 m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<Physics::NavierStokes::EulerVarSet>();
 m_updateVarSet->getModel()->resizePhysicalData(m_dataState);

 cf_assert(m_uInf > 0.);
 cf_assert(m_rhoInf > 0.);
 cf_assert(m_pInf > 0.);

 prepareOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // Create the socket source for the boundaryNormals for each TRS to which this command applies
  // This is because we have different normals in FiniteVolume and FluctSplit / FiniteElement
//   const CFuint nbTrs = getTrsNames().size();
//   for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs)
//   {
//     std::string socketName = getTrsName(iTrs) + "-boundaryNormals";
//     m_sockets.createSocketSource<const CFreal*>(socketName);
//   }

  try {
    setFunction();
  }
  catch (Common::Exception& e) {
    CFout << e.what() << "\n";
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"EulerComputeAeroDG::setFuntion(): alpha function wrongly defined.");

  m_functionParser.Parse(m_function, m_vars);

  if (m_functionParser.ErrorMsg() != 0) {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(m_functionParser.ErrorMsg());
    msg += " Function: " + m_function;
    msg += " Vars: "     + m_vars;
//     throw ParsingErrorException (FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute the value of the angle
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(m_eval);

  // transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;

  // Execute and save file if needed...
  if(!(iter % m_saveRateWall))
  {
      computeWall();
  }

  if(!(iter % m_saveRateAero))
  {
    computeAero();
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::computeWall()
{
//   CFAUTOTRACE;
//
//   prepareOutputFileWall(); // file handle is opened here
//
//   // unused // const CFuint nbEqs = PhysicalModel::getNbEq();
//   // unused //  const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   const CFreal R = m_updateVarSet->getModel()->getR();
//   const CFreal gamma = m_updateVarSet->getModel()->getGamma();
//   const CFreal gammaMinus1 = gamma - 1.;
//   CFreal invRho;
//   CFreal rhoK2;
//   CFreal p;
//   CFreal Cp;
//   CFreal a2;
//   CFreal Mach;
//   CFreal T;
//
//   Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> >
//      geoBuilder = getMethodData().getFaceBuilder();
//
//   // get InnerCells TopologicalRegionSet
//   SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
//
//   // get Inlet Boundary Faces = CurrentTopologicalRegionSet
//   SafePtr<TopologicalRegionSet> faces = getCurrentTRS();
//
//
//   const CFreal time = SubSystemStatusStack::getActive()->getCurrentTime();
//   const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
//
//   cf_assert(m_fileWall->isopen());
//   ofstream& fout = m_fileWall->get();
//
//   //loop over all wall boundary faces
//   for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
//   {
//     const RealVector& state = nstates[(*nodesIdx)[iNState]];
//     Node *const node = nodes[(*nodesIdx)[iNState]];
//
//     invRho = 1./(state)[0];
//     rhoK2 = 0.5*((state)[1]*(state)[1] + (state)[2]*(state)[2])*invRho;
//     p = gammaMinus1*((state)[3] - rhoK2);
//
//     Cp = (p - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);
//
//     a2 = gamma*p/(state)[0];
//     Mach = sqrt((((state)[1]*(state)[1])*invRho + ((state)[2]*(state)[2])*invRho)/a2);
//     T = a2/(gamma*R);
//
//     CFreal TDim = T * m_updateVarSet->getModel()->getTempRef();
//     CFreal pDim = p * m_updateVarSet->getModel()->getPressRef();
//
//     const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
//     const CFreal rhoDim = (state)[0] * rhoRef;
//
//     // Output to File
//     fout
//     << (*node)[XX]  << " "
//     << (*node)[YY]  << " "
//     << m_alphadeg   << " "
//     << iter         << " "
//     << time         << " "
//     << pDim         << " "
//     << Mach         << " "
//     << Cp           << " "
//     << TDim         << " "
//     << rhoDim       << "\n";
//   }
//
//   m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::computeAero()
{
  CFAUTOTRACE;

  //unused//  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  //Get the datahandle containing the boundary Normals
//   SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
//   const std::string TRSName = trs->getName();
//   const std::string socketName = TRSName + "-boundaryNormals";
//
//   DataHandle<const CFreal*> boundaryNormals = m_sockets.
//     getSocketSource<const CFreal*>(socketName)->getDataHandle();
//   DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get Inlet Boundary Faces = CurrentTopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = getCurrentTRS();

  //get number of inlet faces
  const CFuint nbFaces = faces->getLocalNbGeoEnts();

  //   get DGMethodData
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<DiscontGalerkinSolver> DGMethod = spaceMethod.d_castTo<DiscontGalerkinSolver>();
  cf_assert(DGMethod.isNotNull());
  SafePtr<DiscontGalerkinSolverData> DGData = DGMethod->getMethodData().d_castTo<DiscontGalerkinSolverData>();

  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> > geoBuilder = DGData->getFaceBuilder();

  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  RealVector normal(0.0, nbDim);
  RealVector Cp(nbDim);
  Cp = 0.;

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;

  m_lift = .0;
  m_drag = .0;
  State tempState;

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    geoData.idx = iFace;
    GeometricEntity *const currFace = geoBuilder->buildGE();
    //get nodes of face
    std::vector<Node*>&  nodes  = *currFace->getNodes();

    GeometricEntity* cellLeft  = currFace->getNeighborGeo(0);
    //get nodes and states of the cell
    std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();
    //get number of states of the cell
    const CFuint nbStatesInCellLeft = left_cell_states.size();

      CFreal detJacobi;
      if (nbDim == 2)
      {
        RealVector hlp = *nodes[1] - *nodes[0];
        detJacobi = sqrt(hlp[0]*hlp[0]+hlp[1]*hlp[1]);
      }
      else
      {
        //     detJacobi = face.computeVolume()*2;
        RealVector hlp1 = *nodes[1] - *nodes[0];
        RealVector hlp2 = *nodes[2] - *nodes[0];
        detJacobi = sqrt((hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1])*(hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1]) + (hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])*(hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])+(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0])*(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0]));//*2.0/2.0;
      }

      //*****************************************************************
      //*****************************************************************
      //SET Integrator
      //compute shape function in quadrature points
      const std::vector<RealVector>& leftShapeFunctions =  DGData->getContourIntegrator()->getSolutionIntegrator(cellLeft)->computeShapeFunctionsAtQuadraturePoints();

      //numbers of quadrature points
      CFuint m_nbKvadrPoint = DGData->getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

      //set weights for element quadrature
      const std::vector<RealVector>& leftWeight = DGData->getContourIntegrator()->getSolutionIntegrator(cellLeft)->getCoeff();

      //get coordinates of quadrature points
      const std::vector<RealVector>& leftCoord = DGData->getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();

      //compute gradient of shape functions in quadrature points
      std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);

      //*****************************************************************
      //*****************************************************************

      //find local index of face in cell (must be improved)
      CFuint m_idxFaceFromLeftCell= 10;
      if (nbDim == 2)
      {
        for (CFuint i=0; i < 3; i++)
          for (CFuint j=0; j < 2; j++)
            if((*left_cell_nodes[(i)%3]==*nodes[(j)%2])&&(*left_cell_nodes[(i+1)%3]==*nodes[(j+1)%2])) m_idxFaceFromLeftCell=i;
      }
      else
      {
        for (CFuint i=0; i < 4; i++)
          for (CFuint j=0; j < 3; j++)
            if(((*left_cell_nodes[(i)%4]==*nodes[(j)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j+2)%3]))||((*left_cell_nodes[(i)%4]==*nodes[(j+2)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j)%3])))
            {
              m_idxFaceFromLeftCell=i;
            }
      }
      cf_assert(m_idxFaceFromLeftCell!=10);

    // get the face normal
      if (nbDim == 2)
      {
        normal[0]= (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[YY] - (*left_cell_nodes[m_idxFaceFromLeftCell])[YY];
        normal[1]= (*left_cell_nodes[m_idxFaceFromLeftCell])[XX] - (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[XX];
        normal /= sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
      }
      else
      {
        normal = currFace->computeAvgCellNormal();
        RealVector HlpNormal = *left_cell_nodes[(m_idxFaceFromLeftCell+3)%4] - *left_cell_nodes[(m_idxFaceFromLeftCell)%4];
        if ((normal[0]*HlpNormal[0]+normal[1]*HlpNormal[1]+normal[2]*HlpNormal[2]) > 0)
        {
          normal *=-1;
        }
        normal /= sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
      }

      CFreal press = 0.;
      for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
      {
        CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
        {
          tempState[iEq]= 0.;
        }
        //computation of state in point of kvadrature - from previous time step
        for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
        {
          RealVector &state = *left_cell_states[iState]->getData();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
          {
            tempState[iEq] += leftShapeFunctions[leftIndex][iState]*state[iEq];
          }
        }
        cf_assert((tempState[0])>0);
        CFreal invRho = 1./(tempState[0]);
        CFreal p;
        if (nbDim == 3)
        {
          CFreal rhoK2 = 0.5*(tempState[1]*tempState[1] + tempState[2]*tempState[2] + tempState[3]*tempState[3])*invRho;
          p = gammaMinus1*(tempState[4] - rhoK2);
        }
        else
        {
          CFreal rhoK2 = 0.5*(tempState[1]*tempState[1] + tempState[2]*tempState[2])*invRho;
          p = gammaMinus1*(tempState[3] - rhoK2);
        }
//         press += (m_pInf-p);
        press += leftWeight[0][kvadrature_point]*p*detJacobi;
      }
//       press /= m_nbKvadrPoint;

      //Compute Cp (scaled by m, depending on ownership of face to processor)
      //press-pInf are switched to take into account the minus sign in normal

      Cp[XX] += press * normal[XX];
      Cp[YY] += press * normal[YY];
      if (nbDim == 3)
      {
        Cp[ZZ] += press * normal[ZZ];
      }
      // release the face
      geoBuilder->releaseGE();
  }

  // adimensionalize Cp
  Cp /= (0.5 * m_rhoInf * m_uInf * m_uInf * refLength);
  if (nbDim == 3)
  {
    m_lift += (-sin(m_alpharad)*Cp[XX] + cos(m_alpharad)*Cp[ZZ]);
    m_drag += (cos(m_alpharad)*Cp[XX] + sin(m_alpharad)*Cp[ZZ]);
  }
  else
  {
    m_lift += (-sin(m_alpharad)*Cp[XX] + cos(m_alpharad)*Cp[YY]);
    m_drag += (cos(m_alpharad)*Cp[XX] + sin(m_alpharad)*Cp[YY]);
  }

  m_xForceCoef = Cp[XX] / refLength;
  m_yForceCoef = Cp[YY] / refLength;
  if (nbDim == 3) {m_zForceCoef = Cp[ZZ] / refLength;}

  //Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::prepareOutputFileWall()
{
  CFAUTOTRACE;

  cf_assert (!m_fileWall->isopen());

  using namespace boost::filesystem;
  path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
  file = PathAppender::getInstance().appendAllInfo( file );

  ofstream& fout = m_fileWall->open(file);

  fout << "TITLE = Wall_Values" << "\n";
  fout << "VARIABLES = x y Alpha Iter Time Pressure Mach Cp Temperature Density" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::prepareOutputFileAero()
{
  CFAUTOTRACE;

  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  convergenceFile << "TITLE = Aerodynamic_Coeficients"  << "\n";
  convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef" << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void EulerComputeAeroDG::updateOutputFileAero()
{
  CFAUTOTRACE;

  boost::filesystem::path fpath (m_nameOutputFileAero);
//  fpath = PathAppender::getInstance().appendAllInfo( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  convergenceFile
  << subSysStatus->getNbIter()        << " "
  << subSysStatus->getCurrentTime()   << " "
  << m_alphadeg                       << " "
  << m_lift                           << " "
  << m_drag                           << " \n";
//   << m_xForceCoef                     << " "
//   << m_yForceCoef                     << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoeff

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

