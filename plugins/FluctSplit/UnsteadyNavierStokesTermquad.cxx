#include "Environment/ObjectProvider.hh"

#include "Common/BadValueException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "FluctSplit/UnsteadyNavierStokesTermquad.hh"
#include "FluctSplit/InwardNormalsData.hh"

#include "FluctSplit/NavierStokesTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<UnsteadyNavierStokesTermquad,
		       FluctuationSplitData,
		       ComputeDiffusiveTerm,
		       FluctSplitSpaceTimeNavierStokesModule>
UnsteadyNavierStokesquadDiffusiveTermProvider("UnsteadyNavierStokesquad");

//////////////////////////////////////////////////////////////////////////////

UnsteadyNavierStokesTermquad::UnsteadyNavierStokesTermquad(const std::string& name) :
  ComputeDiffusiveTerm(name),
  _diffVar(CFNULL),
  _updateVar(CFNULL),
  _radius(0.),
  _states(),
  _values(),
  _avValues(CFNULL),
  _avState(CFNULL),
  _normal(), 
  _gradRho(),
  _gradRhoU(),
  _gradRhoV(),
  _edata(),
  socket_pastStates("pastStates")
{
}

//////////////////////////////////////////////////////////////////////////////

UnsteadyNavierStokesTermquad::~UnsteadyNavierStokesTermquad()
{
  deletePtr(_avValues);
  deletePtr(_avState);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::setDiffusiveVarSet(SafePtr<DiffusiveVarSet> diffVar)
{
  _diffVar = diffVar.d_castTo<NavierStokesVarSet>();
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::setUpdateVarSet(SafePtr<ConvectiveVarSet> updateVar)
{
  _updateVar = updateVar.d_castTo<EulerVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::computeDiffusiveTerm
(GeometricEntity *const geo, vector<RealVector>& result, bool updateCoeffFlag)
{ 

  vector<State*> *const cellStates = geo->getStates();
  nbCellStates = cellStates->size();
  ovNbStatesInCell = 1.0/nbCellStates;
  nbEqs = PhysicalModelStack::getActive()->getNbEq();
  cellID = geo->getID();
 DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
 volume = socket_volumes.getDataHandle()[cellID];
 for (CFuint i = 0; i < nbCellStates; ++i) {
    for (CFuint j = 0; j < nbEqs; ++j) {
      result[i][j] = 0.0;
    }
  }


DistributionData& ddata = getMethodData().getDistributionData();
 // Start with the part from the past
  DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
  // get the paststates in this cell
  for (CFuint i = 0; i < nbCellStates; ++i) {
    const CFuint stateID = (*ddata.states)[i]->getLocalID();
    _states[i] = pastStatesStorage[stateID];
  }

  
  computeDiffRes(false,result,false);



   for (CFuint i = 0; i < nbCellStates; ++i) {
     _states[i] = (*cellStates)[i];
   }

   computeDiffRes(true,result,updateCoeffFlag);


}
/////////////////////////////////////////////////////////////////////////////////


void UnsteadyNavierStokesTermquad::computeDiffRes(bool isPast, vector<RealVector>& result, bool updateCoeffFlag)
{
  CFreal dt = SubSystemStatusStack::getActive()->getDT();
  dtCoeff = -dt/2.0;

  DistributionData& ddata = getMethodData().getDistributionData();
  const CFuint nbStatesInCell = ddata.states->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  NSTerm& model = _diffVar->getModel();
 // CF_DEBUG_POINT;
  
  if (nbStatesInCell ==4){
// CF_DEBUG_POINT;
  const Node& node0 = (*ddata.states)[0]->getCoordinates();
  const Node& node1 = (*ddata.states)[1]->getCoordinates();
  const Node& node2 = (*ddata.states)[2]->getCoordinates();
  const Node& node3 = (*ddata.states)[3]->getCoordinates();
 //CF_DEBUG_POINT;
  CFreal detJ = (1.0/8.0)*((node3[YY]-node1[YY])*(node2[XX]-node0[XX])-(node2[YY]-node0[YY])*(node3[XX]-node1[XX]));
  CFreal ovdetJ = 1.0/detJ;
  CFreal J11 = 2.0*(node1[XX]-node0[XX]);
  CFreal J12 = 2.0*(node1[YY]-node0[YY]);
  CFreal J21 = 2.0*(node3[XX]-node0[XX]);
  CFreal J22 = 2.0*(node3[YY]-node0[YY]);
  CFreal onethird = 1.0/3.0;
  CFreal onesixth = 1.0/6.0;


   computeAverageState(nbStatesInCell, ovNbStatesInCell, nbEqs, _states, *_avState);
  this->getMethodData().getDistributionData().avState = *_avState;
  getMethodData().getUpdateVar()->computePhysicalData(*_avState, _edata);
  CFreal rhobar = _edata[EulerTerm::RHO];
  CFreal ubar = _edata[EulerTerm::VX];
  CFreal vbar = _edata[EulerTerm::VY];
 
  const CFreal R =_updateVar->getModel()->getR();

  const CFreal Rinv = 1.0/R;
  
  CFreal T0 = (0.4*Rinv/(*_states[0])[0])*((*_states[0])[3]-0.5*((*_states[0])[1]*(*_states[0])[1]+(*_states[0])[2]*(*_states[0])[2])/(*_states[0])[0]);
  CFreal T1 = (0.4*Rinv/(*_states[1])[0])*((*_states[1])[3]-0.5*((*_states[1])[1]*(*_states[1])[1]+(*_states[1])[2]*(*_states[1])[2])/(*_states[0])[0]);
  CFreal T2 = (0.4*Rinv/(*_states[2])[0])*((*_states[2])[3]-0.5*((*_states[2])[1]*(*_states[2])[1]+(*_states[2])[2]*(*_states[2])[2])/(*_states[0])[0]);
  CFreal T3 = (0.4*Rinv/(*_states[3])[0])*((*_states[3])[3]-0.5*((*_states[3])[1]*(*_states[3])[1]+(*_states[3])[2]*(*_states[3])[2])/(*_states[0])[0]);

  // CF_DEBUG_OBJ(T0);
  // CF_DEBUG_OBJ((*_states[0])[3]);
  // CF_DEBUG_OBJ((*_states[0])[2]);
  // CF_DEBUG_OBJ((*_states[0])[1]);
  // CF_DEBUG_OBJ((*_states[0])[0]);

  // CF_DEBUG_OBJ(T1);
  // CF_DEBUG_OBJ(T2);
  // CF_DEBUG_OBJ(T3);


  CFreal u0 = (((*_states[0])[1] - (*_states[0])[0]*ubar)/rhobar);
  CFreal u1 = (((*_states[1])[1] - (*_states[1])[0]*ubar)/rhobar);
  CFreal u2 = (((*_states[2])[1] - (*_states[2])[0]*ubar)/rhobar);
  CFreal u3 = (((*_states[3])[1] - (*_states[3])[0]*ubar)/rhobar);

  CFreal v0 = (((*_states[0])[2] - (*_states[0])[0]*vbar)/rhobar);
  CFreal v1 = (((*_states[1])[2] - (*_states[1])[0]*vbar)/rhobar);
  CFreal v2 = (((*_states[2])[2] - (*_states[2])[0]*vbar)/rhobar);
  CFreal v3 = (((*_states[3])[2] - (*_states[3])[0]*vbar)/rhobar);

// CF_DEBUG_POINT;
  /******************************* Flux*dl0 ********************************************************************/

  CFreal dudxdlidx = (u0*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
		      u1*(-J22*J22*onethird+J12*J12*onesixth) +
		      u2*(-J22*J22*onesixth+0.5*J12*J22-J12*J12*onesixth) +
		      u3*(J22*J22*onesixth-J12*J12*onethird))*ovdetJ;
 
  CFreal dudydlidx = (u0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      u1*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
		      u2*(J22*J21*onesixth -J12*J21*0.25 -J11*J22*0.25 +J12*J11*onesixth) +
		      u3*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird))*ovdetJ;
 

  CFreal dvdxdlidx = (v0*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
		      v1*(-J22*J22*onethird+J12*J12*onesixth) +
		      v2*(-J22*J22*onesixth+0.5*J12*J22-J12*J12*onesixth) +
		      v3*(J22*J22*onesixth-J12*J12*onethird))*ovdetJ;
 
  CFreal dvdydlidx = (v0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      v1*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
		      v2*(J22*J21*onesixth -J12*J21*0.25 -J11*J22*0.25 +J12*J11*onesixth) +
		      v3*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird))*ovdetJ;

  CFreal dudxdlidy = (u0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      u1*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
		      u2*(J22*J21*onesixth-J21*J12*0.25-J22*J11*0.25+J12*J11*onesixth) +
		      u3*(-J21*J22*onesixth+J11*J22*0.25-J21*J12*0.25+J11*J12*onethird))*ovdetJ;

  CFreal dudydlidy = (u0*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
		      u1*(-J21*J21*onethird+J11*J11*onesixth) +
		      u2*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
		      u3*(J21*J21*onesixth-J11*J11*onethird))*ovdetJ;

  CFreal dvdxdlidy = (v0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      v1*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
		      v2*(J22*J21*onesixth-J21*J12*0.25-J22*J11*0.25+J12*J11*onesixth) +
		      v3*(-J21*J22*onesixth+J11*J22*0.25-J21*J12*0.25+J11*J12*onethird))*ovdetJ;

  CFreal dvdydlidy = (v0*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
		      v1*(-J21*J21*onethird+J11*J11*onesixth) +
		      v2*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
		      v3*(J21*J21*onesixth-J11*J11*onethird))*ovdetJ;

  CFreal dTdxdlidy = (T0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      T1*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
		      T2*(J22*J21*onesixth-J21*J21*0.25-J22*J11*0.25+J12*J11*onesixth) +
		      T3*(-J21*J22*onesixth+J11*J22*0.25-J21*J12*0.25+J11*J12*onethird))*ovdetJ;

  CFreal dTdydlidy = (T0*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
		      T1*(-J21*J21*onethird+J11*J11*onesixth) +
		      T2*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
		      T3*(J21*J21*onesixth-J11*J11*onethird))*ovdetJ;

  CFreal dTdxdlidx = (T0*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
		      T1*(-J22*J22*onethird+J12*J12*onesixth) +
		      T2*(-J22*J22*onesixth+0.5*J12*J22-J12*J12*onesixth) +
		      T3*(J22*J22*onesixth-J12*J12*onethird))*ovdetJ;

  CFreal dTdydlidx = (T0*(-J22*J21*onethird+J12*J21*0.25+J22*J11*0.25-J12*J11*onethird) +
		      T1*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
		      T2*(J22*J21*onesixth -J12*J21*0.25 -J11*J22*0.25 +J12*J11*onesixth) +
		      T3*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird))*ovdetJ;

  vector<RealVector*>& gradients =
    getMethodData().getDistributionData().gradients;

  
  // unit normal
 
  
  CFreal mu = _diffVar->getDynViscosity(*_avState, gradients);
  CFreal lambda =  _diffVar->getThermConductivity(*_avState, mu);;
   // CF_DEBUG_OBJ(lambda);

  result[0][0] += 0.0;
  result[0][1] += (dtCoeff)*mu*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)+dudydlidy+dvdxdlidy);
  result[0][2] += (dtCoeff)*mu*(dudydlidx+dvdxdlidx+(4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy));
  result[0][3] += (dtCoeff)*(mu*(ubar*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)) + vbar*(dudydlidx+dvdxdlidx)) - lambda*dTdxdlidx +
			     mu*(ubar*(dudydlidy+dvdxdlidy) + vbar*((4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy))) - lambda*dTdydlidy);     
 
  //CF_DEBUG_OBJ(result[0]);

// CF_DEBUG_POINT;
 /******************************* Flux*dl1 ********************************************************************/

  dudxdlidx = (u0*(-J22*J22*onethird+J12*J12*onesixth) +
	       u1*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5) +
	       u2*(J22*J22*onesixth-J12*J12*onethird) +
	       u3*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth))*ovdetJ;

  dudydlidx = (u0*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
	       u1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       u2*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       u3*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth))*ovdetJ;

  dvdxdlidx = (v0*(-J22*J22*onethird+J12*J12*onesixth) +
	       v1*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5) +
	       v2*(J22*J22*onesixth-J12*J12*onethird) +
	       v3*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth))*ovdetJ;

  dvdydlidx = (v0*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
	       v1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       v2*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       v3*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth))*ovdetJ;


  dudxdlidy = (u0*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
	       u1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       u2*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       u3*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth))*ovdetJ;


  dudydlidy = (u0*(-J21*J21*onethird+J11*J11*onesixth) +
	       u1*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5) +
	       u2*(J21*J21*onesixth-J11*J11*onethird) +
	       u3*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth))*ovdetJ;


  dvdxdlidy = (v0*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
	       v1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       v2*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       v3*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth))*ovdetJ;


  dvdydlidy = (v0*(-J21*J21*onethird+J11*J11*onesixth) +
	       v1*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5) +
	       v2*(J21*J21*onesixth-J11*J11*onethird) +
	       v3*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth))*ovdetJ;   

  dTdxdlidy = (T0*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J12*J11*onesixth) +
	       T1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       T2*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       T3*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth))*ovdetJ;

  dTdydlidy = (T0*(-J21*J21*onethird+J11*J11*onesixth) +
	       T1*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5) +
	       T2*(J21*J21*onesixth-J11*J11*onethird) +
	       T3*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth))*ovdetJ;  

  dTdxdlidx = (T0*(-J22*J22*onethird+J12*J12*onesixth) +
	       T1*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5) +
	       T2*(J22*J22*onesixth-J12*J12*onethird) +
	       T3*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth))*ovdetJ;

  dTdydlidx = (T0*(J21*J22*onethird-J11*J22*0.25+J12*J21*0.25-J12*J11*onesixth) +
	       T1*(-J22*J21*onethird-J12*J21*0.25-J22*J11*0.25-J11*J12*onethird) +
	       T2*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       T3*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth))*ovdetJ;
/*check3 pas de probleme*/

	 //CF_DEBUG_POINT;	     
  
  result[1][0] += 0.0;
  result[1][1] += (dtCoeff)*mu*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)+dudydlidy+dvdxdlidy);
  result[1][2] += (dtCoeff)*mu*(dudydlidx+dvdxdlidx+(4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy));
  result[1][3] += (dtCoeff)*(mu*(ubar*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)) + vbar*(dudydlidx+dvdxdlidx)) - lambda*dTdxdlidx +
			     mu*(ubar*(dudydlidy+dvdxdlidy) + vbar*((4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy))) - lambda*dTdydlidy);     

    //    CF_DEBUG_OBJ(result[1]);

 /******************************* Flux*dl2 ********************************************************************/

 

  dudxdlidx = (u0*(-J22*J22*onesixth+J12*J22*0.5-J12*J12*onesixth) +
	       u1*(J22*J22*onesixth-J12*J12*onethird) +
	       u2*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
	       u3*(-J22*J22*onethird+J12*J12*onesixth))*ovdetJ;

  dudydlidx = (u0*(J22*J21*onesixth-J21*J12*0.25-J22*J11*0.25+J12*J11*onesixth) +
	       u1*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       u2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       u3*(J21*J22*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth))*ovdetJ;

  dvdxdlidx = (v0*(-J22*J22*onesixth+J12*J22*0.5-J12*J12*onesixth) +
	       v1*(J22*J22*onesixth-J12*J12*onethird) +
	       v2*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
	       v3*(-J22*J22*onethird+J12*J12*onesixth))*ovdetJ;

  dvdydlidx = (v0*(J22*J21*onesixth-J21*J12*0.25-J22*J11*0.25+J12*J11*onesixth) +
	       v1*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       v2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       v3*(J21*J22*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth))*ovdetJ;

  dudxdlidy = (u0*(J22*J21*onesixth-J12*J21*0.25-J11*J22*0.25+J12*J11*onesixth) +
	       u1*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       u2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       u3*(J21*J22*onethird-J11*J22*0.25+J21*J12*0.25-J11*J12*onesixth))*ovdetJ;

  dudydlidy = (u0*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
	       u1*(J21*J21*onesixth-J11*J11*onethird) +
	       u2*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
	       u3*(-J21*J21*onethird+J11*J11*onesixth))*ovdetJ;

  dvdxdlidy = (v0*(J22*J21*onesixth-J12*J21*0.25-J11*J22*0.25+J12*J11*onesixth) +
	       v1*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       v2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       v3*(J21*J22*onethird-J11*J22*0.25+J21*J12*0.25-J11*J12*onesixth))*ovdetJ;

  dvdydlidy = (v0*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
	       v1*(J21*J21*onesixth-J11*J11*onethird) +
	       v2*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
	       v3*(-J21*J21*onethird+J11*J11*onesixth))*ovdetJ;   

  dTdxdlidy = (T0*(J22*J21*onesixth-J12*J21*0.25-J11*J22*0.25+J12*J11*onesixth) +
	       T1*(-J22*J21*onesixth-J12*J21*0.25+J11*J22*0.25+J12*J11*onethird) +
	       T2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       T3*(J21*J22*onethird-J11*J22*0.25+J21*J12*0.25-J11*J12*onesixth))*ovdetJ;
  
  dTdydlidy = (T0*(-J21*J21*onesixth+J21*J11*0.5-J11*J11*onesixth) +
	       T1*(J21*J21*onesixth-J11*J11*onethird) +
	       T2*(J21*J21*onethird+J11*J11*onethird-J11*J21*0.5) +
	       T3*(-J21*J21*onethird+J11*J11*onesixth))*ovdetJ;

  dTdxdlidx = (T0*(-J22*J22*onesixth+J12*J22*0.5-J12*J12*onesixth) +
	       T1*(J22*J22*onesixth-J12*J12*onethird) +
	       T2*(J22*J22*onethird+J12*J12*onethird-J22*J12*0.5) +
	       T3*(-J22*J22*onethird+J12*J12*onesixth))*ovdetJ;

  dTdydlidx = (T0*(J22*J21*onesixth-J21*J12*0.25-J22*J11*0.25+J12*J11*onesixth) +
	       T1*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J11*J12*onethird) + /*check3 probleme*/
	       T2*(-J22*J21*onethird-J11*J12*onethird+J12*J21*0.25+J11*J22*0.25) +
	       T3*(J21*J22*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth))*ovdetJ;

  result[2][0] += 0.0;
  result[2][1] += (dtCoeff)*mu*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)+dudydlidy+dvdxdlidy);
  result[2][2] += (dtCoeff)*mu*(dudydlidx+dvdxdlidx+(4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy));
  result[2][3] += (dtCoeff)*(mu*(ubar*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)) + vbar*(dudydlidx+dvdxdlidx)) - lambda*dTdxdlidx +
			     mu*(ubar*(dudydlidy+dvdxdlidy) + vbar*((4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy))) - lambda*dTdydlidy);     

  //  CF_DEBUG_OBJ(result[2]);
 //CF_DEBUG_POINT;
 /******************************* Flux*dl3 ********************************************************************/

  dudxdlidx = (u0*(J22*J22*onesixth-J12*J12*onethird) +
	       u1*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth) +
	       u2*(-J22*J22*onethird+J12*J12*onesixth) +
	       u3*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5))*ovdetJ;

  dudydlidx = (u0*(-J21*J22*onesixth+J11*J22*0.25-J12*J21*0.25+J11*J12*onethird) +
	       u1*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth) +
	       u2*(J21*J22*onethird-J11*J22*0.25+J21*J12*0.25-J11*J12*onesixth) +
	       u3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;

  dvdxdlidx = (v0*(J22*J22*onesixth-J12*J12*onethird) +
	       v1*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth) +
	       v2*(-J22*J22*onethird+J12*J12*onesixth) +
	       v3*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5))*ovdetJ;

  dvdydlidx = (v0*(-J21*J22*onesixth+J11*J22*0.25-J12*J21*0.25+J11*J12*onethird) +
	       v1*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth) +
	       v2*(J21*J22*onethird-J11*J22*0.25+J21*J22*0.25-J11*J12*onesixth) +
	       v3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;

  dudxdlidy = (u0*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird) +
	       u1*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth) +
	       u2*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth) +
	       u3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;
 
  dudydlidy = (u0*(J21*J21*onesixth-J11*J11*onethird) +
	       u1*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth) +
	       u2*(-J21*J21*onethird+J11*J11*onesixth) +
	       u3*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5))*ovdetJ;

  dvdxdlidy = (v0*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird) +
	       v1*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth) +
	       v2*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth) +
	       v3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;

  dvdydlidy = (v0*(J21*J21*onesixth-J11*J11*onethird) +
	       v1*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth) +
	       v2*(-J21*J21*onethird+J11*J11*onesixth) +
	       v3*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5))*ovdetJ;


  dTdxdlidy = (T0*(-J22*J21*onesixth+J12*J21*0.25-J11*J22*0.25+J12*J11*onethird) +
	       T1*(J22*J21*onesixth+J12*J21*0.25+J11*J22*0.25+J12*J11*onesixth) +
	       T2*(J22*J21*onethird-J12*J21*0.25+J11*J22*0.25-J11*J12*onesixth) +
	       T3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;

  dTdydlidy = T0*((J21*J21*onesixth-J11*J11*onethird) +
		  T1*(-J21*J21*onesixth-J21*J11*0.5-J11*J11*onesixth) +
		  T2*(-J21*J21*onethird+J11*J11*onesixth) +
		  T3*(J21*J21*onethird+J11*J11*onethird+J21*J11*0.5))*ovdetJ;

  dTdxdlidx = (T0*(J22*J22*onesixth-J12*J12*onethird) +
	       T1*(-J22*J22*onesixth-J12*J22*0.5-J12*J12*onesixth) +
	       T2*(-J22*J22*onethird+J12*J12*onesixth) +
	       T3*(J22*J22*onethird+J12*J12*onethird+J22*J12*0.5))*ovdetJ;

/*check3 pas de probleme*/

  dTdydlidx = (T0*(-J21*J22*onesixth+J11*J22*0.25-J12*J21*0.25+J11*J12*onethird) +
	       T1*(J21*J22*onesixth+J11*J22*0.25+J12*J21*0.25+J12*J11*onesixth) +
	       T2*(J21*J22*onethird-J11*J22*0.25+J21*J12*0.25-J11*J12*onesixth) +
	       T3*(-J22*J21*onethird-J12*J11*onethird-J11*J22*0.25-J12*J21*0.25))*ovdetJ;

/*check2*/
  /*************CHECK UNTILL THERE *****************/

  result[3][0] += 0.0;
  result[3][1] += (dtCoeff)*mu*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)+dudydlidy+dvdxdlidy);
  result[3][2] += (dtCoeff)*mu*(dudydlidx+dvdxdlidx+(4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy));
  result[3][3] += (dtCoeff)*(mu*(ubar*((4.0/3.0)*dudxdlidx-(2.0/3.0)*(dvdydlidx)) + vbar*(dudydlidx+dvdxdlidx)) - lambda*dTdxdlidx +
			     mu*(ubar*(dudydlidy+dvdxdlidy) + vbar*((4.0/3.0)*dvdydlidy-(2.0/3.0)*(dudxdlidy))) - lambda*dTdydlidy);  

  // CF_DEBUG_OBJ(result[3]);
 //CF_DEBUG_POINT;
  if (isPast) {
    DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
    const CFreal dcoeff = (updateCoeffFlag) ? mu/rhobar*detJ : 0.0;
     if (updateCoeffFlag) {
       for (CFuint i = 0; i < nbCellStates; ++i) {
	 updateCoeff[(*ddata.states)[i]->getLocalID()] += sqrt(volume)*dcoeff;
       }
     }
}

  }

  else {
    //CF_DEBUG_OBJ(nbStatesInCell);
    ovNbStatesInCell = 1.0/nbStatesInCell;
    computeAverageState(nbStatesInCell, ovNbStatesInCell, nbEqs, _states, *_avState);
//CF_DEBUG_POINT;
    this->getMethodData().getDistributionData().avState = *_avState;
//CF_DEBUG_POINT;
    getMethodData().getUpdateVar()->computePhysicalData(*_avState, _edata);
//CF_DEBUG_POINT;
    computeCellGradientsAndAverageStatequad( _edata);
  //CF_DEBUG_POINT;
    DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
    DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  //CF_DEBUG_POINT;
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    NSTerm& model = _diffVar->getModel();
    vector<RealVector*>& gradients = ddata.gradients;
  //CF_DEBUG_POINT;
    // set the diffusive term
    for (CFuint i = 0; i < nbCellStates; ++i) {
      // this is not the unit normal !!
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	_normal[iDim] = normals[cellID]->getNodalNormComp(i,iDim);
      }
      const RealVector& flux = _diffVar->getFlux(*_avState, gradients, _normal, 0.0);
    result[i] = (dtCoeff)*flux;
    //CF_DEBUG_POINT;
  const CFreal dcoeff = (updateCoeffFlag) ? (model.getPhysicalData())[NSTerm::MU]/
    (_edata[EulerTerm::RHO]*this->socket_volumes.getDataHandle()[cellID]*dim*dim) : 0.0;
      if (updateCoeffFlag) {
      const CFreal faceArea = normals[cellID]->getAreaNode(i);
      updateCoeff[(*ddata.states)[i]->getLocalID()] += faceArea*faceArea*dcoeff;
    }
  }
  }
  //CF_DEBUG_POINT;
} 
//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::setup()
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodesInControlVolume =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  
  _states.resize(nbNodesInControlVolume);
  _values.resize(nbEqs, nbNodesInControlVolume);
  _normal.resize(dim); 
  _gradRho.resize(dim);
  _gradRhoU.resize(dim);
  _gradRhoV.resize(dim);
  
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
    resizePhysicalData(_edata); 
  _avValues = new State();
  _avState = new State();
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::configure ( Config::ConfigArgs& args )
{
  ComputeDiffusiveTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadyNavierStokesTermquad::computeCellGradientsAndAverageStatequad( const RealVector& pdata)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  const CFuint nbCellStates = 3;
  //CF_DEBUG_POINT;
  // compute vars that will be used to compute the gradients
  _diffVar->setGradientVars(_states, _values, nbCellStates);
  //CF_DEBUG_POINT;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dimCoeff = 1./dim;
  const CFreal coeffGrad = dimCoeff/volume;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  //CF_DEBUG_POINT;
  vector<RealVector*>& gradients =
    getMethodData().getDistributionData().gradients;
  //CF_DEBUG_POINT;
  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    RealVector& grad = *gradients[iEq];
    grad = 0.0;
    
    for (CFuint is = 0; is < nbCellStates; ++is) {
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
	grad[iDim] += _values(iEq,is)*normals[cellID]->getNodalNormComp(is,iDim);
      }
    }
    grad *= coeffGrad;
  }
  // CF_DEBUG_POINT;
  // new velocity gradients
  _gradRho = 0.0;
  // CF_DEBUG_POINT;
  for (CFuint is = 0; is < nbCellStates; ++is) {
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      _gradRho[iDim] += (*_states[is])[0]*normals[cellID]->getNodalNormComp(is,iDim);
    }
  }
  _gradRho *= coeffGrad;
  
  _gradRhoU = 0.0;
  _gradRhoV = 0.0;
  const CFreal avRho = pdata[EulerTerm::RHO];
  const CFreal avU = pdata[EulerTerm::VX];
  const CFreal avV = pdata[EulerTerm::VY];
  
  for (CFuint is = 0; is < nbCellStates; ++is) {
    for (CFuint iDim = 0; iDim < dim; ++iDim) {    
      _gradRhoU[iDim] += (*_states[is])[1]*normals[cellID]->getNodalNormComp(is,iDim);
      _gradRhoV[iDim] += (*_states[is])[2]*normals[cellID]->getNodalNormComp(is,iDim);	
    }
  }
  _gradRhoU *= coeffGrad;
  _gradRhoV *= coeffGrad;
  
  RealVector& gradU = *gradients[1];
  RealVector& gradV = *gradients[2];    
  gradU = (_gradRhoU-avU*_gradRho)/avRho;
  gradV = (_gradRhoV-avV*_gradRho)/avRho;
  
  SafePtr<EulerTerm> eulerModel = _updateVar->getModel();
  if(avRho < 0.0) {
    cout << "negative avRho = " << avRho << " in cell "<< cellID << endl;
    cout << "pdata = " << pdata << endl;
    throw BadValueException (FromHere(),"Negative average density in cell");
  }
  //CF_DEBUG_POINT;
  if (pdata[EulerTerm::A] < 0.0) {
    cout << "negative a = " << pdata[EulerTerm::A] << " in cell "<< cellID << endl;
    cout << "pdata = " << pdata << endl;
    throw BadValueException (FromHere(),"Negative average sound speed in cell");
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
