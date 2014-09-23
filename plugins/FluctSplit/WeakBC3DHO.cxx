#include "WeakBC3DHO.hh"
#include "InwardNormalsData.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

WeakBC3DHO::WeakBC3DHO(const std::string& name) :
  WeakBC(name),
  socket_isUpdated("isUpdated"),
  m_fluxes()
{
}

//////////////////////////////////////////////////////////////////////////////

WeakBC3DHO::~WeakBC3DHO()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC3DHO::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WeakBC::needsSockets();

  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3DHO::setup()
{
  WeakBC::setup();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // size the fluxes
  // triangular faces ONLY for  now
  m_fluxes.resize(6);
  for (CFuint i = 0; i < m_fluxes.size(); ++i) {
    m_fluxes[i].resize(nbEqs);
  }

// //......... DEBUG ...............
// Fluxfile.open("Fluxes_FField.dat", ios::app);
// Fluxfile << " TITLE  =  Convergence of FField " << endl;
// Fluxfile << " VARIABLES = Iter Res[0] " << endl;
// //...............................

}

//////////////////////////////////////////////////////////////////////////////

void WeakBC3DHO::executeOnTrs()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//  const CFreal fifth = 1./5.;
//   const CFreal halfO.minAlphafifth = 0.5*(1. - m_alpha)*fifth;
//   const CFreal alphafifth = m_alpha*fifth;

  // states forming the ghost cell
  // @todo this works only for triangular boundary faces
  vector<State*> states(2);


//  // get the data handle for the states
//  DataHandle<State*> allStates = socket_states.getDataHandle();
//-----------NEW-------------
  // get the data handle for the states
  DataHandle<State*,GLOBAL> allStates = socket_states.getDataHandle();
//---------------------------
 const CFreal third = 1./3.;
  const CFreal halfOEminAlphaThird = 0.5*(1. - m_alpha)*third;
  const CFreal alphaThird = m_alpha*third;
//   const CFreal oneMinusAlpha = 1. - m_alpha;
//    const CFreal third = 1./6.;

  // get the data handle for the is Updated flags
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  // get the data handle for the rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // create two ghost states (no local ID is needed for these states)
  State* gstate0 = new State();
  State* gstate1 = new State();
  State* gstate2 = new State();
  State* gstate3 = new State();
  State* gstate4 = new State();
  State* gstate5 = new State();

  Common::SafePtr<Framework::TopologicalRegionSet> const faces = getCurrentTRS();

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = faces;
  const CFuint nbFaces = faces->getLocalNbGeoEnts();


//.......for Fluxes file...........................................................
Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
//.................................................................................


//cycle on all the faces of the TRS......................
   for (CFuint iFace = 0; iFace < nbFaces; ++iFace) 
{

    // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity& currFace = *geoBuilder->buildGE();

  //  cf_assert(currFace.nbNodes() == 3);

    // set the face normal
    setFaceNormal(currFace.getID());

    State *const state0 = currFace.getState(0);
    State *const state1 = currFace.getState(1);
    State *const state2 = currFace.getState(2);
    State *const state3 = currFace.getState(3);
    State *const state4 = currFace.getState(4);
    State *const state5 = currFace.getState(5);

    const CFuint stateID0 = state0->getLocalID();
    const CFuint stateID1 = state1->getLocalID();
    const CFuint stateID2 = state2->getLocalID();
    const CFuint stateID3 = state3->getLocalID();
    const CFuint stateID4 = state4->getLocalID();
    const CFuint stateID5 = state5->getLocalID();

  
    // set the appropriate values in the ghost states
    setGhostState(*state0, *gstate0);
    setGhostState(*state1, *gstate1);
    setGhostState(*state2, *gstate2);
    setGhostState(*state3, *gstate3);
    setGhostState(*state4, *gstate4);
    setGhostState(*state5, *gstate5);

//...................here I compute the nodal fluxes
//...........................................................................
    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate0;
    states[1] = state0;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[0]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate1;
    states[1] = state1;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[1]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate2;
    states[1] = state2;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[2]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate3;
    states[1] = state3;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[3]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate4;
    states[1] = state4;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[4]);

    // set ghost state == first state
    // corresponding B-state == second state
    states[0] = gstate5;
    states[1] = state5;
    // compute the corresponding fluctuation and the update coeff
    computeFlux(states, m_fluxes[5]);
//......................................................................


// //....................... DEBUG ......................
//          SumFlux = 0.0;
// 
//          for (CFuint is = 0; is < 6; ++is) 
//           {
//            SumFlux += m_fluxes[is][0];
//           }
// //....................................................



// //..........HERE I CALCULATE THE AVERAGE PRESSURE in the face states........
// //..........and then in the ghost state.....................................
// 
//       CFreal P_av = 0.0;
//       CFreal P_av_G = 0.0;
// 
//       for (CFuint is = 0; is < 6; ++is) 
//        {
// 
//        State *const state = currFace.getState(is);
// 
//        const CFreal rho = (*state)[0];
//        const CFreal u = (*state)[1]/rho;
//        const CFreal v = (*state)[2]/rho;
//        const CFreal w = (*state)[3]/rho;
// 
//        const CFreal gamma = 1.4;
//        const CFreal halfGammaMinus1 = 0.2;
//        const CFreal rhoV2 = rho*(u*u + v*v + w*w);
// 
//        P_av += (gamma-1.)* (*state)[4] - halfGammaMinus1*rhoV2;
// 
//        }
//  
//        P_av = P_av/6.0; 
// //        cout << "P_av " << P_av <<endl;
// 
// {     //here I calculate ghost Pressure...........
//        const CFreal rho = (*gstate0)[0];
//        const CFreal u = (*gstate0)[1]/rho;
//        const CFreal v = (*gstate0)[2]/rho;
//        const CFreal w = (*gstate0)[3]/rho;
// 
//        const CFreal gamma = 1.4;
//        const CFreal halfGammaMinus1 = 0.2;
//        const CFreal rhoV2 = rho*(u*u + v*v + w*w);
// 
//        P_av_G += (gamma-1.)* (*gstate0)[4] - halfGammaMinus1*rhoV2;
// //        cout << "P_av_G " << P_av_G <<endl;
// }
// 
// //.......Now, I calculate the damping factor............
// 
//       CFreal Damp_Factor;
// 
//       if ( (P_av_G/P_av < 0.0) || (P_av_G/P_av > 2.0))
//        {
//         Damp_Factor = 0.0;
//        }
//       else
//        {
//         CFreal PI = acos(-1.0);
//         Damp_Factor = (1.0 - sin( P_av_G/P_av * PI + PI/2.0) )/2.0;
//        }
// 
// //......................................................



    for (CFuint iEq = 0; iEq < nbEqs; ++iEq) 
     {


//   CFreal c0= 1.0/6.0;
//   CFreal c1= 1.0/6.0;
//   CFreal c2= 1.0/6.0;
//   CFreal c3= 1.0/6.0;
//   CFreal c4= 1.0/6.0;
//   CFreal c5= 1.0/6.0;
// 
// 
//      if (!isUpdated[stateID0]) 
//       {
//       rhs(stateID0, iEq, nbEqs) +=
// (alphafifth*c0*m_fluxes[0][iEq] + halfO.minAlphafifth*(c1*m_fluxes[1][iEq] + c2*m_fluxes[2][iEq]+ c3*m_fluxes[3][iEq]+ c4*m_fluxes[4][iEq]+ c5*m_fluxes[5][iEq]) );
//       m_flagState[stateID0] = true;
// 	}
// 
//      if (!isUpdated[stateID1]) 
//       {
//       rhs(stateID1, iEq, nbEqs) +=
// (alphafifth*c1*m_fluxes[1][iEq] + halfO.minAlphafifth*(c0*m_fluxes[0][iEq] + c2*m_fluxes[2][iEq]+ c3*m_fluxes[3][iEq]+ c4*m_fluxes[4][iEq]+ c5*m_fluxes[5][iEq]) );
//       m_flagState[stateID1] = true;
// 	}
// 
//      if (!isUpdated[stateID2]) 
//       {
//       rhs(stateID2, iEq, nbEqs) +=
// (alphafifth*c2*m_fluxes[2][iEq] + halfO.minAlphafifth*(c1*m_fluxes[1][iEq] + c0*m_fluxes[0][iEq]+ c3*m_fluxes[3][iEq]+ c4*m_fluxes[4][iEq]+ c5*m_fluxes[5][iEq]) );
//       m_flagState[stateID2] = true;
// 	}
// 
//      if (!isUpdated[stateID3]) 
//       {
//       rhs(stateID3, iEq, nbEqs) +=
// (alphafifth*c3*m_fluxes[3][iEq] + halfO.minAlphafifth*(c1*m_fluxes[1][iEq] + c2*m_fluxes[2][iEq]+ c0*m_fluxes[0][iEq]+ c4*m_fluxes[4][iEq]+ c5*m_fluxes[5][iEq]) );
//       m_flagState[stateID3] = true;
// 	}
// 
//      if (!isUpdated[stateID4]) 
//       {
//       rhs(stateID4, iEq, nbEqs) +=
// (alphafifth*c4*m_fluxes[4][iEq] + halfO.minAlphafifth*(c1*m_fluxes[1][iEq] + c2*m_fluxes[2][iEq]+ c3*m_fluxes[3][iEq]+ c0*m_fluxes[0][iEq]+ c5*m_fluxes[5][iEq]) );
//       m_flagState[stateID4] = true;
// 	}
// 
//      if (!isUpdated[stateID5]) 
//       {
//       rhs(stateID5, iEq, nbEqs) +=
// (alphafifth*c5*m_fluxes[5][iEq] + halfO.minAlphafifth*(c1*m_fluxes[1][iEq] + c2*m_fluxes[2][iEq]+ c3*m_fluxes[3][iEq]+ c4*m_fluxes[4][iEq]+ c0*m_fluxes[0][iEq]) );
//       m_flagState[stateID5] = true;
// 	}


  CFreal c0= 0.0;
  CFreal c1= 0.0;
  CFreal c2= 0.0;
  CFreal c3= 1.0/3.0;
  CFreal c4= 1.0/3.0;
  CFreal c5= 1.0/3.0;

  /**Sous-triangle 1: 0-3-5 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID0]) {
        rhs(stateID0, iEq, nbEqs) +=
          alphaThird*m_fluxes[0][iEq] + halfOEminAlphaThird*
          (m_fluxes[3][iEq] + m_fluxes[5][iEq]);

        m_flagState[stateID0] = true;
      }

      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) +=
          alphaThird*m_fluxes[3][iEq] + halfOEminAlphaThird*
          (m_fluxes[0][iEq] + m_fluxes[5][iEq]);

        m_flagState[stateID3] = true;
      }

      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) +=
          alphaThird*m_fluxes[5][iEq] + halfOEminAlphaThird*
          (m_fluxes[0][iEq] + m_fluxes[3][iEq]);

        m_flagState[stateID5] = true;
      }
    }
  
   /**Sous-triangle 2: 3-1-4 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) +=
          alphaThird*m_fluxes[3][iEq] + halfOEminAlphaThird*
          (m_fluxes[1][iEq] + m_fluxes[4][iEq]);

        m_flagState[stateID3] = true;
      }

      if (!isUpdated[stateID1]) {
        rhs(stateID1, iEq, nbEqs) +=
          alphaThird*m_fluxes[1][iEq] + halfOEminAlphaThird*
          (m_fluxes[3][iEq] + m_fluxes[4][iEq]);

        m_flagState[stateID1] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) +=
          alphaThird*m_fluxes[4][iEq] + halfOEminAlphaThird*
          (m_fluxes[1][iEq] + m_fluxes[3][iEq]);

        m_flagState[stateID4] = true;
      }
    }
  
   /**Sous-triangle 3: 5-4-2 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) +=
          alphaThird*m_fluxes[5][iEq] + halfOEminAlphaThird*
          (m_fluxes[4][iEq] + m_fluxes[2][iEq]);

        m_flagState[stateID5] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) +=
          alphaThird*m_fluxes[4][iEq] + halfOEminAlphaThird*
          (m_fluxes[5][iEq] + m_fluxes[2][iEq]);

        m_flagState[stateID4] = true;
      }

      if (!isUpdated[stateID2]) {
        rhs(stateID2, iEq, nbEqs) +=
          alphaThird*m_fluxes[2][iEq] + halfOEminAlphaThird*
          (m_fluxes[5][iEq] + m_fluxes[4][iEq]);

        m_flagState[stateID2] = true;
      }
    }
  
  /**Sous-triangle 3: 5-4-3 **/

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
      // distribute corrections to the boundary nodes
      if (!isUpdated[stateID5]) {
        rhs(stateID5, iEq, nbEqs) -=
          alphaThird*m_fluxes[5][iEq] + halfOEminAlphaThird*
          (m_fluxes[4][iEq] + m_fluxes[3][iEq]);

        m_flagState[stateID5] = true;
      }

      if (!isUpdated[stateID4]) {
        rhs(stateID4, iEq, nbEqs) -=
          alphaThird*m_fluxes[4][iEq] + halfOEminAlphaThird*
          (m_fluxes[5][iEq] + m_fluxes[3][iEq]);

        m_flagState[stateID4] = true;
      }

      if (!isUpdated[stateID3]) {
        rhs(stateID3, iEq, nbEqs) -=
          alphaThird*m_fluxes[3][iEq] + halfOEminAlphaThird*
          (m_fluxes[5][iEq] + m_fluxes[4][iEq]);

        m_flagState[stateID3] = true;
      }
    }
     }
  /*
  
  
     if (!isUpdated[stateID0]) {
      rhs(stateID0, iEq, nbEqs) += third*(m_alpha*m_fluxes[0][iEq]+oneMinusAlpha*0.2*(m_fluxes[1][iEq]+m_fluxes[2][iEq]+m_fluxes[3][iEq]+m_fluxes[4][iEq]+m_fluxes[5][iEq]));
      m_flagState[stateID0] = true;
	}

     if (!isUpdated[stateID1]) {
      rhs(stateID1, iEq, nbEqs) += third*(m_alpha*m_fluxes[1][iEq]+oneMinusAlpha*0.2*( m_fluxes[0][iEq]+m_fluxes[2][iEq]+m_fluxes[3][iEq]+m_fluxes[4][iEq]+m_fluxes[5][iEq]));
      m_flagState[stateID1] = true;
	}

     if (!isUpdated[stateID2]) {
      rhs(stateID2, iEq, nbEqs) += third*(m_alpha*m_fluxes[2][iEq]+oneMinusAlpha*0.2*( m_fluxes[1][iEq]+m_fluxes[0][iEq]+m_fluxes[3][iEq]+m_fluxes[4][iEq]+m_fluxes[5][iEq]));
      m_flagState[stateID2] = true;
	}

     if (!isUpdated[stateID3]) {
      rhs(stateID3, iEq, nbEqs) += third*(m_alpha*m_fluxes[3][iEq]+oneMinusAlpha*0.2*( m_fluxes[1][iEq]+m_fluxes[2][iEq]+m_fluxes[0][iEq]+m_fluxes[4][iEq]+m_fluxes[5][iEq]));
      m_flagState[stateID3] = true;
	}

     if (!isUpdated[stateID4]) {
      rhs(stateID4, iEq, nbEqs) += third*(m_alpha*m_fluxes[4][iEq]+oneMinusAlpha*0.2*( m_fluxes[1][iEq]+m_fluxes[2][iEq]+m_fluxes[3][iEq]+m_fluxes[0][iEq]+m_fluxes[5][iEq]));
      m_flagState[stateID4] = true;
	}

     if (!isUpdated[stateID5]) {
      rhs(stateID5, iEq, nbEqs) += third*(m_alpha*m_fluxes[5][iEq]+oneMinusAlpha*0.2*( m_fluxes[1][iEq]+m_fluxes[2][iEq]+m_fluxes[3][iEq]+m_fluxes[4][iEq]+m_fluxes[0][iEq]));
      m_flagState[stateID5] = true;
	}*/

//   /*  }*/



  //Release the geometric entity
  geoBuilder->releaseGE();
  }

// //..................Fluxes file writing ........................
// Fluxfile << subSysStatus->getNbIter() << "  " << SumFlux <<endl;
// //..............................................................


  // this could be done in a more efficient way ...
//   const CFuint nbStates = m_flagState.size();
//   for (CFuint i = 0; i < nbStates; ++i) {
//     if (m_flagState[i] == true) {
//       isUpdated[i] = true;
//     }
//   }

  deletePtr(gstate0);
  deletePtr(gstate1);
  deletePtr(gstate2);
  deletePtr(gstate3);
  deletePtr(gstate4);
  deletePtr(gstate5);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
