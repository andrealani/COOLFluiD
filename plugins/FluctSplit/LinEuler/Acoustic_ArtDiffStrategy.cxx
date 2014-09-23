#include "Acoustic_ArtDiffStrategy.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctSplit.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<Acoustic_ArtDiffStrategy,
                       FluctuationSplitData,
                       ArtificialDiffusionStrategy,
                       FluctSplitModule>
                       acousticArtDiffStrategyProvider("Acoustic");

//////////////////////////////////////////////////////////////////////////////

void Acoustic_ArtDiffStrategy::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Artnu","Set up the cutoff frequency of the filter.");
}

//////////////////////////////////////////////////////////////////////////////

Acoustic_ArtDiffStrategy::Acoustic_ArtDiffStrategy(const std::string& name) :
  ArtificialDiffusionStrategy(name)
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_nu = 1e-5;
  setParameter("Artnu",&m_nu);
}

//////////////////////////////////////////////////////////////////////////////

Acoustic_ArtDiffStrategy::~Acoustic_ArtDiffStrategy()
{
}

////////////////////////////////////////////////////////////////////////////

void Acoustic_ArtDiffStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  ArtificialDiffusionStrategy::setup();

  m_dim = Framework::PhysicalModelStack::getActive()->getDim();
  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_stateNormal.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState) {
      m_stateNormal[iState].resize(m_dim);
  }

}

////////////////////////////////////////////////////////////////////////////

/**
 * 
 * @param residual 
 */

void Acoustic_ArtDiffStrategy::addArtificialDiff(std::vector<RealVector>& residual)
{

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbNodes = residual.size();


  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();
  vector<State*>& states = *ddata.states;
  const CFreal volume = ddata.cell->computeVolume();

//   RealVector c;
//   c.resize(nbNodes);
  RealVector r;
  r.resize(nbNodes);
  RealVector Rs;
  Rs.resize(nbNodes);
  RealVector Damping;
  Damping.resize(nbEqs);
  RealVector SmoothRes;
  SmoothRes.resize(nbEqs);

  const CFuint cell_ID = getMethodData().getDistributionData().cellID;


//    for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
//      for (CFuint iState = 0; iState < nbNodes; ++iState) {
//        m_stateNormal[iState][iDim] = m_normals[cell_ID]->getNodalNormComp(iState,iDim);
//        }
//    }


//**************************************like NS **************************************************
//  CFreal normals=0.;
//
//   for (CFuint iState = 0; iState < nbNodes; ++iState) {
//     residual[iState] = 0.;
//     Rs[iState] = m_nu/volume;
//   }
// 
// 
//    for (CFuint node = 0; node < nbNodes; ++node) {
//      if(!isBState[states[node]->getLocalID()]) {
//        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
//          Damping[iEq] = 0.;
//          SmoothRes[iEq] = 0.;
//          for (CFuint other = 0; other < nbNodes; ++other) {
//            for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
//              normals = m_stateNormal[node][iDim]*m_stateNormal[other][iDim];
//            }
//            SmoothRes[iEq] += normals*(*states[other])[iEq];
//          }
//        Damping[iEq] = Rs[node]*SmoothRes[iEq];
//        residual[node][iEq] = Rs[node]*Damping[iEq];
//        }
//      }
//      else {
//        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
//          residual[node][iEq] = 0.0;
//      }
//    }

//************************************** like NS **************************************************



//************************************** FD like **************************************************

  for (CFuint iState = 0; iState < nbNodes; ++iState) {
    residual[iState] = 0.;
    Rs[iState] = m_nu/volume;
//     c[iState] = 0.;
//     for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
//       c[iState] += m_stateNormal[iState][iDim];
//     }
  }

   for (CFuint node = 0; node < nbNodes; ++node) {
     const Node& currnode = states[node]->getCoordinates();
     if(!isBState[states[node]->getLocalID()]) {

       for (CFuint other = 0; other < nbNodes; ++other) {
         r[other] =0.;
       }
 
       for (CFuint other = 0; other < nbNodes; ++other) {
          const Node& othernode = states[other]->getCoordinates();
          for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
            r[other] += (currnode[iDim]-othernode[iDim])*(currnode[iDim]-othernode[iDim]);
         }
 	 r[other] = -sqrt(r[other]);
       }

       CFreal sum = 0.;
       for (CFuint other = 0; other < nbNodes; ++other)
         sum += r[other];
       r[node] = -sum;

       for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          Damping[iEq] = 0.;
          SmoothRes[iEq] = 0.;

          for (CFuint other = 0; other < nbNodes; ++other) {
            SmoothRes[iEq] += r[other]*((*states[other])[iEq]);
          }
//           Damping[iEq] = Rs[node]*c[node]*SmoothRes[iEq];
          Damping[iEq] = Rs[node]*SmoothRes[iEq];
          residual[node][iEq] = Damping[iEq];
        }
      }
      else {
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          residual[node][iEq] = 0.0;
      }
    }

//************************************artificial **************************************************


}

//////////////////////////////////////////////////////////////////////////////
    }// End namespace FluctSplit

}// End namespace COOLFluiD
