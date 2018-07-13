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

   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
   DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
   const bool useNodeBased = getMethodData().useNodeBased();
   const CFuint nbStates = (!useNodeBased) ? states.size() : nodes.size();
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

#ifdef CF_HAVE_CUDA
   if (getMethodData().getBuildOnGPU()) {
     mat.updateDiagBlocks(nbStates, nbEqs);
   }
#endif
   
   // assemble the matrix in the HOST
   mat.finalAssembly(rhs.size());
//mat.moveToCPU();
//mat.printToFile("ParalutionMatrix.txt");
//abort();


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
    CFLog(VERBOSE, "StdSolveSys: useGPU " << useGPU << " ; buildOnGPU " << getMethodData().getBuildOnGPU() << " \n");
    if (useGPU){
      if (!getMethodData().getBuildOnGPU()){mat.moveToGPU();}
      rhsVec.moveToGPU();
      solVec.moveToGPU();
    }




    IterCounter++;
   //FOR RE-USE THE LS
   if (firstIter){
      mat.AssignToSolver(ls);
      ls.SetPreconditioner(p);
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
      if (!getMethodData().getBuildOnGPU()){mat.moveToCPU();}
      solVec.moveToCPU();
      rhsVec.moveToCPU();
   }
   mat.LeavePointers();

//mat.resetToZeroEntriesGPU();
//mat.finalAssembly(rhs.size());
//mat.moveToCPU();
//mat.printToFile("ParalutionMatrix.txt");
//abort();





//solVec.printToFile("Sol.txt");
//abort();
   solVec.copy2(&rhs[0], &_upLocalIDs[0], vecSize);


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
