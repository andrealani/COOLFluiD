// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethodMHD.hh"
#include "NewtonMethod/StdUpdateSolPP.hh"

#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Common/BadValueException.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "Framework/ConvergenceMethod.hh"
#include "Framework/ConvergenceMethodData.hh"
#include "Framework/ConvergenceStatus.hh"
#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonIterator.hh"

#include "Framework/PathAppender.hh"
#include "Framework/GeometricEntity.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUpdateSolPP, NewtonIteratorData, NewtonMethodMHDModule> 
stdUpdateSolPPProvider("StdUpdateSolPP");

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSolPP::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<vector<CFreal>,Config::DynamicOption<> >("Relaxation","Relaxation factor");
  options.addConfigOption< bool >("Validate","Check that each update creates variables with physical meaning");
  options.addConfigOption< CFreal >("pressureBoundaryValue", "the boundary-pressure value for the pressure.");
  options.addConfigOption< CFreal >("densityBoundaryValue", "the boundary-density value for the pressure.");
}

//////////////////////////////////////////////////////////////////////////////

StdUpdateSolPP::StdUpdateSolPP(const std::string& name) : 
  NewtonIteratorCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_invalidStates("invalidStates")
{
  addConfigOptionsTo(this);
   
  m_alpha = vector<CFreal>();
  setParameter("Relaxation",&m_alpha);

  m_validate = false;
  setParameter("Validate",&m_validate);

  _pBC = MathTools::MathConsts::CFrealEps();
  setParameter("pressureBoundaryValue", &_pBC);

  _rhoBC = MathTools::MathConsts::CFrealEps();
  setParameter("densityBoundaryValue", &_rhoBC);
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSolPP::setup()
{
  const CFuint nbStates = socket_states.getDataHandle().size();

  if (m_validate) {
    socket_invalidStates.getDataHandle().resize(nbStates);
  }
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  if (m_alpha.size() == 0) {
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = 1.;
    }
  }

  if (m_alpha.size() == 1) {
    const CFreal value = m_alpha[0];
    m_alpha.resize(nbEqs);
    for (CFuint i = 0; i < nbEqs; ++i) {
      m_alpha[i] = value;
    }
  }

  if (m_alpha.size() != nbEqs) {
    throw BadValueException (FromHere(),"StdUpdateSolPP::setup() : m_alpha.size() != nbEqs");
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdUpdateSolPP::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  // rhs is the temporary placeholder for the dU
  DataHandle<CFreal>& dU = rhs;

  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint states_size = states.size();

  SafePtr<FilterState> filterState = getMethodData().getFilterState();
  SafePtr<FilterRHS> filterRHS     = getMethodData().getFilterRHS();
  //>> mark 2024.06.17
  CFuint nbStatesWithNegPressure = 0;
  CFuint nbStatesWithNegrho = 0; 
  //<< mark 2024.06.17  
  //>> mark 2024.10.02
  CFreal T_Sun = 1.5e6;   // K
  CFreal mu_cor = 1.27;   // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal Bref = 2.2e-4;   // T
  CFreal mu0 = 1.2566e-6;
  CFreal rhoref = 1.67e-13;
  CFreal kB = 1.38e-23;
  CFreal pref = 0.03851;
  CFreal vref = 4.80e5;
  CFreal smooth_factor = 0.0;
  //CFreal _rhoBC=2.0;
  //CFreal _pBC = 0.25801090625;
  CFreal Bmagmax=1.0;
  //<< mark 2024.10.02
  const CFreal RSun = 6.9551e8; // m
  const CFreal GMsun = 1.327474512e20; // SI value
  const CFreal g0 = pow(vref, 2) / RSun;
  const CFreal g = -(GMsun / RSun / RSun) / g0; // gravity term dimensionless


  //------------->> Mark 2025.03.11-----------------------
  for (CFuint iState = 0; iState < states_size; ++iState) 
  {
	  State& cur_state = *states[iState];
	  if (cur_state.isParUpdatable() ) 
      //CFreal time_code = SubSystemStatusStack::getActive()->getCurrentTimeDim();	
      //CFreal time_hours = time_code*1447.2/3600.0;
      //CFreal min_decimalhour=(int(time_hours*100.0)-int(time_hours)*100)*0.6;  
	  //if ((cur_state.isParUpdatable() && std::abs(min_decimalhour)<7.0) || (SubSystemStatusStack::getActive()->getNbIter()<2))  
	  {
		  Bmagmax = max(std::sqrt(cur_state[4] * cur_state[4] + cur_state[5] * cur_state[5] + cur_state[6] * cur_state[6]),Bmagmax);
	  }
  } 
  //--------------- Mark 2025.03.11<<---------------------
  for (CFuint iState = 0; iState < states_size; ++iState)
  {
    State& cur_state = *states[iState];
    // do the update only if the state is parallel updatable
    if (cur_state.isParUpdatable())
    {

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	filterRHS->filter(iEq, dU(iState, iEq, nbEqs));
	cur_state[iEq] += m_alpha[iEq] * dU(iState, iEq, nbEqs);
      }
	  //if(plasmaBeta<0.01)
	  //{
	//	cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
	 // }

	  //>> mark 2024.10.02
	  
	  if ((cur_state[7] < (_pBC*1.e-8)) || (cur_state[7] > (_pBC*1.e2))){
		  cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
	  }
	  if ((cur_state[0] < (_rhoBC*1.e-8)) || (cur_state[0] > (_rhoBC*1.e1))){
		  cur_state[0] -= m_alpha[0] * dU(iState, 0, nbEqs);
	  }
	  

	  CFreal x = cur_state.getCoordinates()[XX];
	  CFreal y = cur_state.getCoordinates()[YY];
	  CFreal z = cur_state.getCoordinates()[ZZ];
	  CFreal r2 = x*x + y*y + z*z;
	  CFreal r = std::sqrt(r2);
	  CFreal maxvA = 2.0*Bmagmax*Bref / pow((_rhoBC * rhoref * mu0), 0.5); 
	  //if ((r < 1.2) && (SubSystemStatusStack::getActive()->getNbIter()>=0)){
	  if (SubSystemStatusStack::getActive()->getNbIter()>=0){
		  CFreal Bmag = std::sqrt(cur_state[4] * cur_state[4] + cur_state[5] * cur_state[5] + cur_state[6] * cur_state[6]) + 1.0e-14;
		  //CFreal minPlasmaBeta = 0.001; // 0.02;  
		  CFreal minPlasmaBeta = 0.001;
		  CFreal maxPlasmaBeta = 1.0;
		  CFreal pmag = (Bmag*Bref)*(Bmag*Bref) / 2. / mu0;

		  CFreal PressureBoundary_dimless = cur_state[7];
		  CFreal pth = PressureBoundary_dimless * pref;
		  CFreal plasmaBeta = pth / pmag;

		  CFreal densityBoundary_dimless = cur_state[0];
		  //---->> Mark 2025.03.28---------------------------------------------
		  CFreal Temperature_ratio = PressureBoundary_dimless / densityBoundary_dimless;
		  //------ Mark 2025.03.28<<-------------------------------------------
		  //CFreal maxvA = 3.0e6; //6.0e6; // 
		  CFreal vA = Bmag*Bref / pow((densityBoundary_dimless * rhoref * mu0), 0.5);

		  smooth_factor = 0.0;
		  CFreal dist = (vA - maxvA) / 2.0e3;
		  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

		  CFreal rho_temp = (smooth_factor * (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0 + densityBoundary_dimless * rhoref * (1.0 - smooth_factor)) / rhoref;
		  cur_state[0] = rho_temp;
		  //---->> Mark 2025.03.28---------------------------------------------
		  //cur_state[7] = rho_temp*Temperature_ratio;
		  //------ Mark 2025.03.28<<-------------------------------------------
		  CFreal p_temp = cur_state[7];
		  //(pth < 1.e-9)
		  /*if ((plasmaBeta < 0.0001) || (pth < 1.e-5)){
			  //if (p_temp < 1.e-3){
			  cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
			  p_temp = cur_state[7];
			  nbStatesWithNegrho += 1;
		  }
		  */
		 
		  //if (r < 1.01){
			  //cur_state[0] = max(_rhoBC*std::pow(1.0 / r, 5), _rhoBC*std::pow(1.0 / 22.0, 5));
		//	  cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
		 // }
		   
          //if(r < 1.05){
		  smooth_factor = 0.0;
		  pth = p_temp*pref;
		  plasmaBeta = pth / pmag;
		  dist = (minPlasmaBeta - plasmaBeta) / 2.0e-6; //if current beta is 0.019, then dist is positive, smooth factor is close to 1
		  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);  //if the current beta is  0.3, then dist is negative and so smooth factor is 0
		  p_temp = (smooth_factor * pmag * minPlasmaBeta + pth * (1.0 - smooth_factor)) / pref;
		  //---->> Mark 2025.03.28---------------------------------------------
		  //cur_state[0] = p_temp/Temperature_ratio;
		  //------ Mark 2025.03.28<<-------------------------------------------
		  cur_state[7] = p_temp;
		  //}


		  //>> mark 2024.10.02
		  /*
		  if ((cur_state[7] < (_pBC*1.e-8)) || (cur_state[7] > (_pBC*1.e2))){
			  cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
		  }
		  if ((cur_state[0] < (_rhoBC*1.e-8)) || (cur_state[0] > (_rhoBC*1.e1))){
			  cur_state[0] -= m_alpha[0] * dU(iState, 0, nbEqs);
		  }
		  */

		  
		  /*
		  if (r < 2.0){
			  CFreal vr = x / r*cur_state[1] + y / r*cur_state[2] + z / r*cur_state[3];
			  if (vr < -0.3*std::sqrt(-2.0*g*(1.25 - min(r, 1.25)))){
				  cur_state[1] -= m_alpha[1] * dU(iState, 1, nbEqs);
				  cur_state[2] -= m_alpha[2] * dU(iState, 2, nbEqs);
				  cur_state[3] -= m_alpha[3] * dU(iState, 3, nbEqs);
			  }
		  }
		  */
		  

		  if(r < 1.1+0.002){
		    CFreal V_local=std::sqrt(cur_state[1] * cur_state[1] + cur_state[2] * cur_state[2] + cur_state[3] * cur_state[3])*vref;
	        CFreal C_sound=std::sqrt(1.667*p_temp*pref/(rho_temp*rhoref));
			//rho_temp=cur_state[0];
			//vA = Bmag*Bref / pow((rho_temp * rhoref * mu0), 0.5);
			//C_sound=max(C_sound,vA); 
			CFreal r_fac = 1.0;
			
			dist = 0.3 + std::tanh((r - 1.0)*8.68*r_fac);
	        dist=min(dist,1.0);
	        smooth_factor = (C_sound/V_local)*dist; ///2.d0	
            smooth_factor=min(smooth_factor,1.0); 
			
	        //smooth_factor = C_sound/V_local/2.0;
            //smooth_factor=min(smooth_factor,1.0);
	        cur_state[1]=cur_state[1]*smooth_factor;
	        cur_state[2]=cur_state[2]*smooth_factor;
	        cur_state[3]=cur_state[3]*smooth_factor;
		  }



	  }

	   
	  
	  
	  
//<< mark 2024.10.02
  
//>> mark 2024.06.17
/*
    CFreal Bref = 2.2e-4;   
	CFreal pref = 0.03851;
    CFreal rhoref = 1.67e-13;
    CFreal mu0 = 1.2566e-6;
	if(SubSystemStatusStack::getActive()->getNbIter()>20){
	  CFreal Bmag = std::sqrt(cur_state[4]*cur_state[4] + cur_state[5]*cur_state[5] + cur_state[6]*cur_state[6]); 
	  CFreal pmag = (Bmag*Bref)*(Bmag*Bref)/2./mu0;
	  CFreal PressureBoundary_dimless = cur_state[7]; 
      CFreal pth = PressureBoundary_dimless * pref;
	  CFreal plasmaBeta = pth/pmag;
      CFreal pressure_CorrectionVal=1.e-7;
	  //if (cur_state[7] < 1.e-10) { 
	//  if ((plasmaBeta < 0.02) || (pth < 1.e-9)) {
	//	cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
	//	nbStatesWithNegPressure += 1;
    //  }  	
  CFreal smooth_factor = 0.0; 
  CFreal densityBoundary_dimless = cur_state[0];
  CFreal maxvA = 2.0e6;
  CFreal vA = Bmag*Bref / pow((densityBoundary_dimless * rhoref * mu0),0.5);

  smooth_factor = 0.0;
  CFreal dist = (vA - maxvA)/2.0e3;
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

  CFreal rho_temp = (smooth_factor * (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0 + densityBoundary_dimless * rhoref * (1.0 - smooth_factor)) / rhoref;
      if ((cur_state[0] < pressure_CorrectionVal) || (plasmaBeta < 0.0001) || (pth < 1.e-9)) {
          cur_state[0] -= m_alpha[0] * dU(iState, 0, nbEqs);
		  cur_state[7] -= m_alpha[7] * dU(iState, 7, nbEqs);
		  rho_temp = cur_state[0];
		  nbStatesWithNegrho += 1;
	  }	
	  cur_state[0] = rho_temp;
	}
*/
//<< mark 2024.06.17	
      // apply a polymorphic filter to the state
      filterState->filter(cur_state);
    }
    else {
      // reset to 0 the RHS for ghost states in order to avoid 
      // inconsistencies in the parallel L2 norm computation
      cf_assert(!cur_state.isParUpdatable());
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
	dU(iState, iEq, nbEqs) = 0.;
      }
    }
  }
//>> mark 2024.06.17	
  if (nbStatesWithNegPressure > 0)
    cout << "There were " << nbStatesWithNegPressure << " states with negative pressure." << endl;
  if (nbStatesWithNegrho > 0)
    cout << "There were " << nbStatesWithNegrho << " states with negative rho." << endl;  
//<< mark 2024.06.17	  
  ///This vector will temporally contain the ID of any vector considered not valid.
  std::vector<CFuint> badStatesIDs;

  ///Name of the file where to store invalid states in.
  boost::filesystem::path filepath ( "unphysical_states.plt");

  filepath = PathAppender::getInstance().appendParallel( filepath );

  // loop for unphysicalness check:
  if ( m_validate )
  {
    Common::SafePtr<SpaceMethod> theSpaceMethod = getMethodData().getCollaborator<SpaceMethod>();
    Common::SafePtr<SpaceMethodData> theSpaceMethodData = theSpaceMethod->getSpaceMethodData();
    Common::SafePtr<Framework::ConvectiveVarSet> theVarSet = theSpaceMethodData->getUpdateVar();

    ///Reset the vector:
    badStatesIDs.resize(0);

    for (CFuint iState = 0; iState < states_size; ++iState)
    {
      State& cur_state = *states[iState];
      if ( !( theVarSet->isValid(cur_state) ) )
      {
        badStatesIDs.push_back( iState );
      }
    } // for all states
  } // if validate

  if (badStatesIDs.size() > 0)
  {
    cf_assert(m_validate);
    // sort the invalid states in order to be able to search later on 
    sort(badStatesIDs.begin(), badStatesIDs.end());
    
    DataHandle<CFreal> invalidStates = socket_invalidStates.getDataHandle();
    cf_assert(invalidStates.size() == states_size);
    invalidStates = 0.;
    for (CFuint i = 0; i < badStatesIDs.size(); ++i) {
      invalidStates[badStatesIDs[i]] = 1.;
    }
    
    CFLog(VERBOSE, "StdUpdateSolPP::execute() => [" << badStatesIDs.size() << "] invalid states detected\n");
    // correct all unphysical states 
    correctUnphysicalStates(badStatesIDs);
    
    ///Gets the global iteration
    CFuint cur_global_iter = SubSystemStatusStack::getActive()->getNbIter();
    ///Gets the Newton method procedure iteration
    CFuint cur_iter = getMethodData().getConvergenceStatus().iter;

    ///opens a file to put the states which have unphysical variables
    ofstream fout ( filepath.string().c_str(), ios::app );

    if ( !fout.is_open() )
    {
      cout << endl << "File [" << filepath.string().c_str() << "]  is not open.\n";
    }
    ///
    const std::vector<std::string>& theVarNames=getMethodData().getCollaborator<SpaceMethod>()->getSpaceMethodData()->getUpdateVar()->getVarNames();
    CFuint dim = PhysicalModelStack::getActive()->getDim();
    ///

    ///Zone header:
    fout << "VARIABLES = \"X\"\n ";
    fout << "\"Y\"" << endl;
    if (dim == 3 ){
        fout << "\"Z\"" << endl;
    }
    for (CFuint i = 0; i < theVarNames.size(); i++){
      fout << "\"" << theVarNames[i] << "\"" << endl;
    }
    ///


    CFuint badStates_size = badStatesIDs.size();

    ///Line  
    fout << "ZONE T = \"Global iteration "<< cur_global_iter << ", Newton iteration " << cur_iter << "\"" << endl;
    ///Line
    fout << "I="<< badStates_size << ", J=1, K=1, ZONETYPE=Ordered"<< endl;
    ///Line
    fout << "DATAPACKING=POINT" << endl;
    ///Line
    fout << "DT=(";
    for ( CFuint i = 0 ; i < dim + theVarNames.size(); i++){
      fout << "DOUBLE ";
    }

    fout << ")" << endl;
    ///

    for (CFuint index = 0; index < badStates_size; ++index)
    {
     State& invalid_state  = *states[ badStatesIDs[index] ];
     fout << invalid_state.getCoordinates() << " " << invalid_state << std::endl;
    } // for all invalid states

    fout.close();

  } // writing badStates

  // reset to 0 the update coefficient
  updateCoeff = 0.0;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdUpdateSolPP::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> > StdUpdateSolPP::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_invalidStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
