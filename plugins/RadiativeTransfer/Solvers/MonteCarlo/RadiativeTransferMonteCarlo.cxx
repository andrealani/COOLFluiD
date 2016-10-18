#include "RadiativeTransfer/Solvers/MonteCarlo/RadiativeTransferMonteCarlo.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTrackingAxi.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking3D.hh"
#include "LagrangianSolver/ParticleTracking/ParticleTracking2D.hh"

//////////////////////////////////////////////////////////////////////////////
//
namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<RadiativeTransferMonteCarlo<ParticleTrackingAxi>,
                      DataProcessingData,     RadiativeTransferModule>
RadiativeTransferMonteCarloAxiFVMCCProvider("RadiativeTransferMonteCarloAxiFVMCC");


MethodCommandProvider<RadiativeTransferMonteCarlo<ParticleTracking3D>,
                      DataProcessingData,    RadiativeTransferModule>
RadiativeTransferMonteCarlo3DFVMCCProvider("RadiativeTransferMonteCarlo3DFVMCC");


MethodCommandProvider<RadiativeTransferMonteCarlo<ParticleTracking2D>,
                      DataProcessingData,    RadiativeTransferModule>
RadiativeTransferMonteCarlo2DFVMCCProvider("RadiativeTransferMonteCarlo2DFVMCC");


//////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////
///
//void RadiativeTransferMonteCarlo::myReduceSpectra(const RealMatrix& indata,
//                                                  RealMatrix& outdata)
//{
//  CFLog(INFO, "RadiativeTransferMonteCarlo::myReduceSpectra() => start\n");

//  cf_always_assert(indata.size() > 0);
//  cf_always_assert(outdata.size() > 0);
//  cf_always_assert(outdata.nbRows() == indata.nbRows());

//  // here implement the spectra reduction
//  m_deltaWav = 1e-10;

//  const CFreal minWav     = m_radLibrary->getMinWavelength() *m_deltaWav;
//  const CFreal maxWav     = m_radLibrary->getMaxWavelength() *m_deltaWav;
//  const CFreal nInStride  = m_radLibrary->getWavelengthStride();
//  const CFreal nOutStride = m_reducedSpectraSize;

//  vector<CFreal> inWav(nInStride);
//  vector<CFreal> inEm(nInStride);
//  vector<CFreal> inAm(nInStride);

//  const CFreal dWav= (maxWav-minWav ) / nOutStride;

//  CFreal sum1=0.;
//  for(CFuint j=0; j<indata.nbCols();j+=3){
//    sum1 += indata(0,j+1);
//    //cout<< indata(0,j)<<' '<<indata(0,j+1)<<' '<<indata(0,j+2)<<endl;
//  }
//  cout<<"sum: "<<sum1<<endl;

//  for (CFuint s = 0; s < outdata.nbRows(); ++s) {

//    //extract the vectors from the indata Matrix
//    for(CFuint i=0; i<nInStride; ++i){
//      inWav[i] = indata(s,i*3  );
//      inEm [i] = indata(s,i*3+1);
//      inAm [i] = indata(s,i*3+2);
//    }

//    vector<CFreal> tempOutWav,tempOut Am,tempOutEm;
//    tempOutWav.reserve(nInStride/nOutStride);
//    tempOutAm.reserve(nInStride/nOutStride);
//    tempOutEm.reserve(nInStride/nOutStride);

//    for( CFuint outIdx= 0, inIdx=0; outIdx<nOutStride; ++outIdx){
//      CFreal outWavMin = minWav    + dWav*outIdx;
//      CFreal outWavMax = outWavMin + dWav;

//      if(inIdx>0){
//        //add first point at outWavMin
//        CFreal startWavEm=linearInterpol(inWav[inIdx],inEm[inIdx],inWav[inIdx+1],inEm[inIdx+1],outWavMin);
//        CFreal startWavAm=linearInterpol(inWav[inIdx],inAm[inIdx],inWav[inIdx+1],inAm[inIdx+1],outWavMin);

//        tempOutWav.push_back( outWavMin );
//        tempOutEm.push_back( startWavEm );
//        tempOutAm.push_back( startWavAm );
//        if(inWav[inIdx+1]<outWavMax){
//          ++inIdx;
//        }
//      }
//      while(inWav[inIdx+1]<outWavMax && inIdx<inWav.size()){
//        //add middle points
//        tempOutWav.push_back( inWav[inIdx] );
//        tempOutAm.push_back(  inAm[inIdx]  );
//        tempOutEm.push_back(  inEm[inIdx]  );
//        ++inIdx;
//      }

//      if(inIdx<inWav.size() ){
//        //add last point at outWavMax
//        CFreal endWavEm=linearInterpol(inWav[inIdx],inEm[inIdx],inWav[inIdx+1],inEm[inIdx+1],outWavMax);
//        CFreal endWavAm=linearInterpol(inWav[inIdx],inAm[inIdx],inWav[inIdx+1],inAm[inIdx+1],outWavMax);

//        tempOutWav.push_back( outWavMax );
//        tempOutEm.push_back( endWavEm );
//        tempOutAm.push_back( endWavAm );
//      }

//      //calculate the integral and add it in the reduced coeff vectors
//      m_wavReduced(s,outIdx) = ( tempOutWav.back()+tempOutWav.front() )/2.;
//      CFreal sumEm=0.,sumAm=0.;
//      for(CFuint i=0; i<tempOutWav.size()-1; ++i){
//        sumEm+= ( tempOutEm[i+1]+tempOutEm[i] )/2.*( tempOutWav[i+1]-tempOutWav[i] );
//        sumAm+= ( tempOutAm[i+1]+tempOutAm[i] )/2.*( tempOutWav[i+1]-tempOutWav[i] );
//      }

//      m_emReduced(s,outIdx) = sumEm;
//      m_amReduced(s,outIdx) = sumAm/( tempOutWav.back()-tempOutWav.front() );

//      //clear the temporary vectors for new calculation
//      tempOutWav.clear();
//      tempOutEm.clear();
//      tempOutAm.clear();
//    }
//  }

//  sum1=0.;
//  for(CFuint j=0; j<m_emReduced.nbCols();++j){
//      //cout<< m_wavReduced(0,j)<<' '<<m_emReduced(0,j)<<endl;
//      sum1+=m_emReduced(0,j);
//  }

//  cout<<"Total radPower "<<3.1415*sum1<<endl;

//  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
//  for (CFuint s = 0; s < m_emReduced.nbRows(); ++s) {
//    m_pdfEm[s]= new vector<CFreal>;
//    m_pdfEm[s]->resize(m_emReduced.nbCols(),0.);

//    //first pass to get the total emission coefficient
//    //use this information to get the cell's total Radiative Power

//    //still serial
//    CFreal emSum1=0.;
//    for (CFuint j=0; j < m_emReduced.nbCols(); ++j ){
//      emSum1+=m_emReduced(s,j);
//    }
//    CFreal Volume= (m_isAxi) ? m_axiVolumes[s] : volumes[s];
//    CFreal Q_out =4.*MathConsts::CFrealPi()*emSum1*Volume;
//    m_gOutRadPowers[s]=Q_out;

//    //second pass to get the [0->1] comulative probably distribution or [0,1] energy fraction
//    if(m_mcSampling){
//      CFreal emSum2=0;
//      for (CFuint j=0; j < m_emReduced.nbCols(); ++j ){
//        emSum2+=m_emReduced(s,j);
//        (*(m_pdfEm[s]))[j]=emSum2/emSum1;
//      }
//    }
//    else{

//      for (CFuint j=0; j < m_emReduced.nbCols(); ++j ){
//        (*(m_pdfEm[s]))[j]=m_emReduced(s,j)/emSum1;
//      }
//    }
//  }
////  for(CFuint j=0; j<(*(m_pdfEm[0])).size();++j){
////      cout<< (*(m_pdfEm[0]))[j]<<endl;
////  }
//  CFLog(INFO, "RadiativeTransferMonteCarlo::myReduceSpectra() => END\n");
//}

} // namespace RadiativeTransfer

} // namespace COOLFluiD

