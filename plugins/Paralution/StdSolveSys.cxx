// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/LSSIdxMapping.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"

#include "Paralution/Paralution.hh"
#include "Paralution/StdSolveSys.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Paralution {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSolveSys, ParalutionLSSData, ParalutionModule> 
stdSolveSysParalutionProvider("StdSolveSys");

//////////////////////////////////////////////////////////////////////////////

StdSolveSys::StdSolveSys(const std::string& name) :
  ParalutionLSSCom(name),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_rhs("rhs"), 
  _upLocalIDs(),
  _upStatesGlobalIDs(),
  IterCounter(0)
{
}

//////////////////////////////////////////////////////////////////////////////

StdSolveSys::~StdSolveSys()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSolveSys::execute()
{
  using namespace paralution;
  CFAUTOTRACE;

  CFLog(VERBOSE, "StdSolveSys::execute() \n");

Stopwatch<WallTime> stopTimer;
stopTimer.start();

   ParalutionLSSData& MethodData = getMethodData();

   DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
   cf_assert(_upLocalIDs.size() == _upStatesGlobalIDs.size());
   const CFuint vecSize = _upLocalIDs.size();
   const CFuint nbEqs = MethodData.getNbSysEquations();
  
   ParalutionMatrix& mat = MethodData.getMatrix();      
   ParalutionVector& rhsVec = MethodData.getRhsVector();
   ParalutionVector& solVec = MethodData.getSolVector();
 //For re-use of the ls
   IterativeLinearSolver<LocalMatrix<CFreal>, LocalVector<CFreal>, CFreal >& ls = MethodData.getKSP();
   Preconditioner<LocalMatrix<CFreal>, LocalVector<CFreal>, CFreal >& p = MethodData.getPreconditioner();  
   bool firstIter = MethodData.getFirstIter();
//Maybe this can be stored in this object, so there is no need to call it every time
   bool useGPU = MethodData.getUseGPU();
   CFuint reBuildRatio = MethodData.getreBuildRatio();


   // assemble the matrix in the HOST
   mat.finalAssembly(rhs.size());

   CFLog(VERBOSE, "StdSolveSys::execute() ==> Matrix assembled! \n");
  // mat.convertToCSR();
  // CFLog(NOTICE, "StdSolveSys::execute() ==> Matrix converted to CSR! \n");
//   // after the first iteration freeze the non zero structure
//   // of the matrix
//   if (SubSystemStatusStack::getActive()->getNbIter() == 1) {
//     mat.freezeNonZeroStructure();
//   }

//   // the rhs is copied into the ParalutionVector for the rhs
//   // _upStatesGlobalIDs[i] is different from _upLocalIDs[i]
//   // in case multiple LSS are used
    for (CFuint i = 0; i < vecSize; ++i) {
    //  cout << "nbSysEq = " << nbEqs << ", upLocalIDs = " << _upLocalIDs[i]
    //  	 << ", upStatesGlobalIDs = " << _upStatesGlobalIDs[i] << endl;
      rhsVec.setValue(_upStatesGlobalIDs[i], rhs[_upLocalIDs[i]]);
    }

//   // assemble the rhs vector in HOST
    rhsVec.Assembly();
    CFLog(VERBOSE, "StdSolveSys::execute() ==> RHS assembled! \n");


    solVec.setValue(0.0);
    solVec.Assembly();

    CFLog(VERBOSE, "StdSolveSys::execute() ==> Solution vector assembled! \n");

    if (useGPU){
      mat.moveToGPU();
      rhsVec.moveToGPU();
      solVec.moveToGPU();
    }

    CFLog(VERBOSE, "StdSolveSys::execute() ==> moveToGPU succesful! \n");


IterCounter++;
  //FOR RE-USE THE LS
if (firstIter){
  ls.Build();
  getMethodData().setFirstIter(false);
}
if (IterCounter%reBuildRatio == 0){
  p.ReBuildNumeric();
  IterCounter=0;
}




rhsVec.Solve(ls, solVec.getVecPtr());
//ls.Clear(); //Only if not reusing the ls

if(useGPU){
solVec.moveToCPU();
mat.moveToCPU();
rhsVec.moveToCPU();
}
//solVec.printToFile("Sol.txt");

  solVec.copy2(&rhs[0], &_upLocalIDs[0], vecSize);
  //mat.resetToZeroEntries(); //needed?


CFLog(NOTICE, "StdSolveSys::execute() took " << stopTimer << "s, with " << ls.GetIterationCount() << " iterations \n");


/*
x.Allocate("x", mat.get_ncol());
x.Zeros();

mat.info();
mat.WriteFileMTX("matrix.mtx");

// Linear Solver
GMRES<LocalMatrix<double>, LocalVector<double>, double > ls;
//preconditioner
ILU<LocalMatrix<double>,LocalVector<double>,double> p;
  p.Set(0);

ls.SetOperator(mat);
 ls.SetPreconditioner(p);

ls.Build();
ls.Verbose(2) ;

ls.Solve(rhs, &x);

ls.Clear();

x.MoveToHost();
for(int counter=0;counter<6;counter++)printf(" x[%i]=%f \n",counter,x[counter]);

*/








//   const CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
//   if (getMethodData().getSaveRate() > 0) {
//     if (getMethodData().isSaveSystemToFile() || (nbIter%getMethodData().getSaveRate() == 0)) { 
//       const string mFile = "mat-iter" + StringOps::to_str(nbIter) + ".dat";
//       mat.printToFile(mFile.c_str());
//       // mat.printToScreen();
      
//       const string vFile = "rhs-iter" + StringOps::to_str(nbIter) + ".dat";
//       rhsVec.printToFile(vFile.c_str());
//       // rhsVec.printToScreen();
//     }
//   }
 
//   // reuse te preconditioner
// #if PETSC_VERSION_MINOR==7
//   ParalutionBool reusePC = ((nbIter-1)%getMethodData().getPreconditionerRate() == 0) ?
//      PETSC_FALSE : PETSC_TRUE;
//   CFLog(VERBOSE, "StdSolveSys::execute() => reusePC [" << reusePC <<"]\n");
//   CHKERRCONTINUE(KSPSetReusePreconditioner(ksp,reusePC));
//   PC& pc = getMethodData().getPreconditioner();
//   CHKERRCONTINUE(PCSetReusePreconditioner(pc,reusePC));
// #endif
 
// #if PETSC_VERSION_MINOR==6 || PETSC_VERSION_MINOR==7
//   CFuint ierr = KSPSetOperators(ksp, mat.getMat(), mat.getMat());
// #else
//   CFuint ierr = KSPSetOperators
//     (ksp, mat.getMat(), mat.getMat(),DIFFERENT_NONZERO_PATTERN);
// #endif
  
//   //This is to allow viewing the matrix structure in X windows
//   //Works only if Paralution is compiled with X
//   if(getMethodData().isShowMatrixStructure())
//   {
//     ierr = MatView(mat.getMat(), PETSC_VIEWER_DRAW_WORLD);
//     CHKERRCONTINUE(ierr);
//   }

//   ierr = KSPSetUp(ksp);
//   CHKERRCONTINUE(ierr);

// #if PETSC_VERSION_MINOR==7
//   if (getMethodData().getPCType()==PCBJACOBI && getMethodData().useGPU()) {
//     ParalutionInt its,nlocal,first;
//     KSP* subksp;
//     PC   subpc;  
//     PC& pc = getMethodData().getPreconditioner();
//     PCBJacobiGetSubKSP(pc,&nlocal,&first,&subksp);
//     for (CFuint i=0; i<nlocal; i++) {
//       KSPGetPC(subksp[i],&subpc);
//       PCSetType(subpc,PCILU);
//     }
//   }
// #endif

//   ierr = KSPSolve(ksp, rhsVec.getVec(), solVec.getVec());
//   CHKERRCONTINUE(ierr);

//   CFint iter = 0;
//   ierr = KSPGetIterationNumber(ksp, &iter);
//   CHKERRCONTINUE(ierr);
  
//   if (nbIter%getMethodData().getKSPConvergenceShowRate() == 0) {
//     CFLog(INFO, "KSP convergence reached at iteration: " << iter << "\n");
//   }
  
//   solVec.copy(&rhs[0], &_upLocalIDs[0], vecSize);
}

//////////////////////////////////////////////////////////////////////////////

void StdSolveSys::setup()
{
  CFAUTOTRACE;

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  const bool useNodeBased = getMethodData().useNodeBased();
  const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
  const CFuint nbEqs = getMethodData().getNbSysEquations();
  
  // count the number of updatable states
  CFuint nbUpdatableStates = 0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if ((!useNodeBased && states[i]->isParUpdatable()) ||
	(useNodeBased && nodes[i]->isParUpdatable())) {
      ++nbUpdatableStates;
    }
  }

  const CFuint vecSize = nbUpdatableStates * nbEqs;
  // indexes for the insertion of elements in a ParalutionVector
  _upStatesGlobalIDs.reserve(vecSize);
  _upLocalIDs.reserve(vecSize);

  const LSSIdxMapping& idxMapping = getMethodData().getLocalToGlobalMapping();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();
  const std::valarray<bool>& maskArray = *getMethodData().getMaskArray();

  // indexes are set to global ID for updatable states
  for (CFuint i = 0; i < nbStates; ++i) {
    const CFuint localID = i*totalNbEqs;
    if ((!useNodeBased && states[i]->isParUpdatable()) ||
	(useNodeBased && nodes[i]->isParUpdatable())) {
      const CFuint sID = (!useNodeBased) ? states[i]->getLocalID() : nodes[i]->getLocalID();
      CFint globalID = static_cast<CFint>(idxMapping.getColID(sID))*nbEqs;
      
      for (CFuint iEq = 0; iEq < totalNbEqs; ++iEq) {
	if (maskArray[iEq]) {
	  _upStatesGlobalIDs.push_back(globalID++);
	  _upLocalIDs.push_back(static_cast<CFint>(localID + iEq));
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > StdSolveSys::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_rhs);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Paralution  

} // namespace COOLFluiD
