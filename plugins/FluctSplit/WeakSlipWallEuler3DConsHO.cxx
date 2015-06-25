#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "WeakSlipWallEuler3DConsHO.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {




    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallEuler3DConsHO, FluctuationSplitData, FluctSplitNavierStokesModule> weakSlipWallEuler3DConsHOProvider("WeakSlipWallEuler3DConsHO");

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DConsHO::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("alpha","distribution coefficient");
   options.addConfigOption< bool >("ExactNormSphere","Use the exact normals for a Sphere");
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DConsHO::WeakSlipWallEuler3DConsHO(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _fluxes(),
  _flagState()
{
   addConfigOptionsTo(this);
  _alpha = 1.0;
   setParameter("alpha",&_alpha);
  m_exact_norm = false;
  setParameter("ExactNormSphere",&m_exact_norm);
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallEuler3DConsHO::~WeakSlipWallEuler3DConsHO()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakSlipWallEuler3DConsHO::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DConsHO::setup()
{

  _varSet->setup();
  // size the fluxes
//   _fluxes.resize(MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell()-1); //why?????? 
                                                                                     //may be it's only for triang and tetra..

  _fluxes.resize(6);

//    cout<<"nb states in faces============= "<< _fluxes.size() <<endl;

  for (CFuint i = 0; i < _fluxes.size(); ++i) {
    _fluxes[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  _flagState.resize(states.size());
  _flagState = false;

// //......... DEBUG ...............
//  CDfile.open("Cd_computation.dat", ios::app);
// Fluxfile << " TITLE  =  Convergence of FField " << endl;
// Fluxfile << " VARIABLES = Iter Res[0] " << endl;
// //...............................

  
 }

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DConsHO::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  RealVector normal(0.0, dim);
  const CFreal third = 1./3.;
  const CFreal halfOminAlphaThird = 0.5*(1. - _alpha)*third;
  const CFreal alphaThird = _alpha*third;
  
//   /*const CFreal oneMinusAlpha = 1. - _alpha;
//    const CFreal sixth = 1./6.;*/
 
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = getCurrentTRS();

  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();


//............................for Fluxes file....................................
   Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
//   DistributionData& distdata = getMethodData().getDistributionData();
//   CFuint CellID = distdata.cellID;
//...............................................................................


//     Pnx = 0.0;  //set to zero the output Pnx

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
     

      geoData.idx = iFace;
      GeometricEntity *const currFace = geoBuilder->buildGE();
 
      const vector<State*> *const statesInFace = currFace->getStates();
   
      const CFuint nbStatesInFace = statesInFace->size();
    
      setFaceNormal(currFace->getID(), normal);       // set the face normal


      for (CFuint is = 0; is < nbStatesInFace; ++is) 
       {
	
         if(m_exact_norm) // EXACT NORMALS
          {
           const CFreal NormOld = normal.norm2();
           Node& node0 = (*statesInFace)[is]->getCoordinates();
           const CFreal x0 = node0[XX];
           const CFreal y0 = node0[YY];
           const CFreal z0 = node0[ZZ];
           const CFreal R = sqrt(x0*x0+y0*y0+z0*z0);
           normal[0] = x0/R;
           normal[1] = y0/R;
           normal[2] = z0/R;

           normal *=  NormOld;  //to re-set the FaceArea..
          }

     // compute the normal fluxes corrections for both the statesof this cell
        computeNormalFlux(*(*statesInFace)[is], normal,_fluxes[is]);

       }
// //------------ CD computation ------------------------
//       for (CFuint is = 0; is < nbStatesInFace; ++is) 
//        {
//         computeCD(*(*statesInFace)[is],  normal);
//        }
// //....................................................




// //....................... DEBUG ......................
//          SumFlux = 0.0;
// 
//          for (CFuint is = 0; is < nbStatesInFace; ++is) 
//           {
//            SumFlux += _fluxes[is][0];
//           }
// //....................................................


    /// @todo implement case when the face is quadrilateral
if (nbStatesInFace == 6) {
      const CFuint stateID0 = (*statesInFace)[0]->getLocalID();
      const CFuint stateID1 = (*statesInFace)[1]->getLocalID();
      const CFuint stateID2 = (*statesInFace)[2]->getLocalID();
      const CFuint stateID3 = (*statesInFace)[3]->getLocalID();
      const CFuint stateID4 = (*statesInFace)[4]->getLocalID();
      const CFuint stateID5 = (*statesInFace)[5]->getLocalID();


// //....................................................
//  CFuint Lab = 0;
//  bool flag;
// 
// for (CFuint ist = 0; ist < 6; ist++)
// {
// 
//    if( (*statesInFace)[ist]->getLocalID() == 146 )
//     {
//      flag = true;
//     }
// }
// //....................................................





  /**Sous-triangle 1: 0-3-5 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID0]) {
        rhs(stateID0, iEq, nbEqs) +=
          alphaThird*_fluxes[0][iEq] + halfOminAlphaThird*
          (_fluxes[3][iEq] + _fluxes[5][iEq]);

        _flagState[stateID0] = true;
      }

      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) +=
          alphaThird*_fluxes[3][iEq] + halfOminAlphaThird*
          (_fluxes[0][iEq] + _fluxes[5][iEq]);

        _flagState[stateID3] = true;
      }

      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) +=
          alphaThird*_fluxes[5][iEq] + halfOminAlphaThird*
          (_fluxes[0][iEq] + _fluxes[3][iEq]);

        _flagState[stateID5] = true;
      }
    }
  
   /**Sous-triangle 2: 3-1-4 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) +=
          alphaThird*_fluxes[3][iEq] + halfOminAlphaThird*
          (_fluxes[1][iEq] + _fluxes[4][iEq]);

        _flagState[stateID3] = true;
      }

      if (!isUpdated[stateID1]) {
        rhs(stateID1, iEq, nbEqs) +=
          alphaThird*_fluxes[1][iEq] + halfOminAlphaThird*
          (_fluxes[3][iEq] + _fluxes[4][iEq]);

        _flagState[stateID1] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) +=
          alphaThird*_fluxes[4][iEq] + halfOminAlphaThird*
          (_fluxes[1][iEq] + _fluxes[3][iEq]);

        _flagState[stateID4] = true;
      }
    }
  
   /**Sous-triangle 3: 5-4-2 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) +=
          alphaThird*_fluxes[5][iEq] + halfOminAlphaThird*
          (_fluxes[4][iEq] + _fluxes[2][iEq]);

        _flagState[stateID5] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) +=
          alphaThird*_fluxes[4][iEq] + halfOminAlphaThird*
          (_fluxes[5][iEq] + _fluxes[2][iEq]);

        _flagState[stateID4] = true;
      }

      if (!isUpdated[stateID2]) {
        rhs(stateID2, iEq, nbEqs) +=
          alphaThird*_fluxes[2][iEq] + halfOminAlphaThird*
          (_fluxes[5][iEq] + _fluxes[4][iEq]);

        _flagState[stateID2] = true;
      }
    }
  
  /**Sous-triangle 3: 5-4-3 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) -=
          alphaThird*_fluxes[5][iEq] + halfOminAlphaThird*
          (_fluxes[4][iEq] + _fluxes[3][iEq]);

        _flagState[stateID5] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) -=
          alphaThird*_fluxes[4][iEq] + halfOminAlphaThird*
          (_fluxes[5][iEq] + _fluxes[3][iEq]);

        _flagState[stateID4] = true;
      }

      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) -=
          alphaThird*_fluxes[3][iEq] + halfOminAlphaThird*
          (_fluxes[5][iEq] + _fluxes[4][iEq]);

        _flagState[stateID3] = true;
      }
    }
/*
// cout<<"**** IN CYCLE*******"<<endl;

      // distribute contributions to the nodes
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {



//   CFreal c_div = 1.0/12.0;

  CFreal c0= 0.0;
  CFreal c1= 0.0;
  CFreal c2= 0.0;
  CFreal c3= 1.0/3.0;
  CFreal c4= 1.0/3.0;
  CFreal c5= 1.0/3.0;


if (!isUpdated[stateID0]) {
   rhs(stateID0, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[0][iEq]+oneMinusAlpha*0.2*(_fluxes[1][iEq]+_fluxes[2][iEq]+_fluxes[3][iEq]+_fluxes[4][iEq]+_fluxes[5][iEq]));
  _flagState[stateID0] = true;
	}

if (!isUpdated[stateID1]) {
   rhs(stateID1, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[1][iEq]+oneMinusAlpha*0.2*( _fluxes[0][iEq]+_fluxes[2][iEq]+_fluxes[3][iEq]+_fluxes[4][iEq]+_fluxes[5][iEq]));
  _flagState[stateID1] = true;
	}

if (!isUpdated[stateID2]) {
   rhs(stateID2, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[2][iEq]+oneMinusAlpha*0.2*( _fluxes[1][iEq]+_fluxes[0][iEq]+_fluxes[3][iEq]+_fluxes[4][iEq]+_fluxes[5][iEq]));
  _flagState[stateID2] = true;
	}

if (!isUpdated[stateID3]) {
   rhs(stateID3, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[3][iEq]+oneMinusAlpha*0.2*( _fluxes[1][iEq]+_fluxes[2][iEq]+_fluxes[0][iEq]+_fluxes[4][iEq]+_fluxes[5][iEq]));
  _flagState[stateID3] = true;
	}

if (!isUpdated[stateID4]) {
   rhs(stateID4, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[4][iEq]+oneMinusAlpha*0.2*( _fluxes[1][iEq]+_fluxes[2][iEq]+_fluxes[3][iEq]+_fluxes[0][iEq]+_fluxes[5][iEq]));
  _flagState[stateID4] = true;
	}

if (!isUpdated[stateID5]) {
   rhs(stateID5, iEq, nbEqs) -=
   sixth*(_alpha*_fluxes[5][iEq]+oneMinusAlpha*0.2*( _fluxes[1][iEq]+_fluxes[2][iEq]+_fluxes[3][iEq]+_fluxes[4][iEq]+_fluxes[0][iEq]));
  _flagState[stateID5] = true;
	}







// if (!isUpdated[stateID0]) {
//    rhs(stateID0, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID0] = true;
// 	}
// 
// if (!isUpdated[stateID1]) {
//    rhs(stateID1, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID1] = true;
// 	}
// 
// if (!isUpdated[stateID2]) {
//    rhs(stateID2, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID2] = true;
// 	}
// 
// if (!isUpdated[stateID3]) {
//    rhs(stateID3, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID3] = true;
// 	}
// 
// if (!isUpdated[stateID4]) {
//    rhs(stateID4, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID4] = true;
// 	}
// 
// if (!isUpdated[stateID5]) {
//    rhs(stateID5, iEq, nbEqs) -=
//    c_div*(c0*_fluxes[0][iEq]+c1*_fluxes[1][iEq]+c2*_fluxes[2][iEq]+c3*_fluxes[3][iEq]+c4*_fluxes[4][iEq]+c5*_fluxes[5][iEq]);
//   _flagState[stateID5] = true;
// 	}

// if (flag ==true) {cout << "rhs_146(ieq="<< iEq<<")= "<<rhs(146, iEq, nbEqs) <<endl;}

      } //end cycle on nb of equation*/






    } //enf IF nbStates ==6*/

    // release the face
    geoBuilder->releaseGE();

  }  //end cycle on TRS faces


//..................CD file writing ........................
//   Pnx /= 6.0;   //because we have 6 states in each face
//   CDfile << subSysStatus->getNbIter() << "  " << Pnx <<endl;
//..........................................................



// //..................Fluxes file writing ........................
// Fluxfile << subSysStatus->getNbIter() << "  " << SumFlux <<endl;
// //..............................................................

//   this could be done in a more efficient way ...
   const CFuint nbStates = _flagState.size();
   for (CFuint i = 0; i < nbStates; ++i) {
     if (_flagState[i] == true) {
       isUpdated[i] = true;
     }
   }
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallEuler3DConsHO::computeNormalFlux(const RealVector& state,
                                                const RealVector& normal,
                                                RealVector& flux) const
 {
   const CFreal nx = normal[0]; 
   const CFreal ny = normal[1]; 
   const CFreal nz = normal[2]; 
   const CFreal rho = state[0];
   const CFreal u = state[1]/rho;
   const CFreal v = state[2]/rho;
   const CFreal w = state[3]/rho;
   const CFreal un = u*nx + v*ny + w*nz;
   const CFreal gamma = _varSet->getModel()->getGamma();
   const CFreal halfGammaMinus1 = 0.5*(gamma - 1);
   const CFreal rhoV2 = rho*(u*u + v*v + w*w);

   flux[0] = rho*un;
   flux[1] = un*state[1];
   flux[2] = un*state[2];
   flux[3] = un*state[3];
   flux[4] = un*(gamma*state[4] - halfGammaMinus1*rhoV2);
 }

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// void WeakSlipWallEuler3DConsHO::computeCD(const RealVector& state,
//                                                 const RealVector& normal)
//  {
//    const CFreal nx = normal[0]; 
// //    const CFreal ny = normal[1]; 
// //    const CFreal nz = normal[2]; 
//    const CFreal rho = state[0];
//    const CFreal u = state[1]/rho;
//    const CFreal v = state[2]/rho;
//    const CFreal w = state[3]/rho;
// //    const CFreal un = u*nx + v*ny + w*nz;
//    const CFreal gamma = _varSet->getModel()->getGamma();
// //    const CFreal halfGammaMinus1 = 0.5*(gamma - 1);
//    const CFreal rhoV2 = rho*(u*u + v*v + w*w);
// 
// //------------------------------- p*nx*A_cell / 1/2 *rho* u^2 * D^2  (and then divided by 6..)
//  Pnx +=  - nx * (  (gamma-1.0)*(state[4]-0.5*rhoV2)  ) / (  0.5*(1.22503)*(158.30559*158.30559)*(0.09)  );
// //-------------------------------
// 
//  }

//////////////////////////////////////////////////////////////////////////////










void WeakSlipWallEuler3DConsHO::setFaceNormal(const CFuint faceID,
                                            RealVector& normal)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  normal[0] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 0);
  normal[1] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 1);
  normal[2] = -normals[cellLocalID]->getFaceNormComp(iFaceLocal, 2);
}

//////////////////////////////////////////////////////////////////////////////

// void WeakSlipWallEuler3DConsHO::configure(const Config::ConfigArgs& args) //OLD
// {
 void WeakSlipWallEuler3DConsHO::configure(Config::ConfigArgs& args) //NEW!!!!!
 {
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
