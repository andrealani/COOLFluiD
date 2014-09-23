// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "Trilinos/Trilinos.hh"
#include "StdParSolveSys.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSIdxMapping.hh"

#include "ml_MultiLevelPreconditioner.h"

//#define TB_TEST_1DPOISSON  // test problem for multigrid hacked into this code

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace Teuchos;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdParSolveSys, TrilinosLSSData, TrilinosModule> stdParSolveSysProvider("StdParSolveSys");

//////////////////////////////////////////////////////////////////////////////

StdParSolveSys::StdParSolveSys(const std::string& name) :
  TrilinosLSSCom(name),
  socket_states("states"),
  socket_rhs("rhs"),
  _lid(NULL),
  _gid(NULL)
{
}

//////////////////////////////////////////////////////////////////////////////

StdParSolveSys::~StdParSolveSys()
{
  delete[] _lid;
  delete[] _gid;
}

//////////////////////////////////////////////////////////////////////////////

void StdParSolveSys::execute()
{
  using Teuchos::rcp;
  using Teuchos::RCP;

  Common::Stopwatch<Common::WallTime> stopTimer;
  stopTimer.start();

  TrilinosMatrix* mat = getMethodData().getMatrix();
  TrilinosVector* rhsVec = getMethodData().getRhsVector();
  TrilinosVector* solVec = getMethodData().getSolVector();

  // assemble the matrix: TODO is this needed ?
  // assemble the matrix:
  if (!(mat->isAssembled())) {
    mat->beginAssembly();
    mat->endAssembly();
  }

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  RCP<Epetra_CrsMatrix> epetra_A=rcp(new Epetra_CrsMatrix(*mat->getMat()));
  RCP<Epetra_Vector>    epetra_x=rcp(new Epetra_Vector(*solVec->getVec()));
  RCP<Epetra_Vector>    epetra_b=rcp(new Epetra_Vector(*rhsVec->getVec()));

  // filling the rhs: TODO make more efficient
  const CFuint nrhs=rhsVec->getLocalSize();
  for (CFuint i=0; i<nrhs; i++) {
    (*epetra_b)[i]=rhs[_lid[i]];
    //rhsVec->setValue(_gid[i],rhs[_lid[i]]);
  }
  const CFuint nsol=solVec->getLocalSize();
  for (CFuint i=0; i<nsol; i++) {
    (*epetra_x)[i]=0.;
    //solVec->setValue(_gid[i],0.);
  }

  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder(getMethodData().getOptionsXML()); // the most important in general setup

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream(); // TODO: decouple from fancyostream to ostream or to C stdout when possible
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  Teuchos::CommandLineProcessor  clp(false); // false: don't throw exceptions

  RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( epetra_A );
  RCP<Thyra::VectorBase<double> >         x = Thyra::create_Vector( epetra_x, A->domain() );
  RCP<const Thyra::VectorBase<double> >   b = Thyra::create_Vector( epetra_b, A->range() );

  // r = b - A*x, initial L2 norm
  double nrm_r=0.;
  CFreal systemResidual=-1.;
  {
    Epetra_Vector epetra_r(*epetra_b);
    Epetra_Vector epetra_A_x(epetra_A->OperatorRangeMap());
    epetra_A->Apply(*epetra_x,epetra_A_x);
    epetra_r.Update(-1.0,epetra_A_x,1.0);
    epetra_r.Norm2(&nrm_r);
  }

  // Reading in the solver parameters from the parameters file and/or from
  // the command line.  This was setup by the command-line options
  // set by the setupCLP(...) function above.
  linearSolverBuilder.readParameters(0); /* out.get() if want confirmation about the xml file within trilinos */
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = linearSolverBuilder.createLinearSolveStrategy(""); // create linear solver strategy
  lowsFactory->setVerbLevel((Teuchos::EVerbosityLevel)(getMethodData().getOutputLevel())); // set verbosity

/*
  // print back default and current settings
  if (opts->trilinos.dumpDefault!=0) {
    fflush(stdout); cout << flush;
    _MMESSAGE_(0,1,"Dumping Trilinos/Stratimikos solver defaults to files: 'trilinos_default.txt' and 'trilinos_default.xml'...\n");
    fflush(stdout); cout << flush;
    std::ofstream ofs("./trilinos_default.txt");
    linearSolverBuilder.getValidParameters()->print(ofs,PLPrintOptions().indent(2).showTypes(true).showDoc(true)); // the last true-false is the deciding about whether printing documentation to option or not
    ofs.flush();ofs.close();
    ofs.open("./trilinos_default.xml");
    Teuchos::writeParameterListToXmlOStream(*linearSolverBuilder.getValidParameters(),ofs);
    ofs.flush();ofs.close();
  }
  if (opts->trilinos.dumpCurrXML!=0) {
    fflush(stdout); cout << flush;
    _MMESSAGE_(0,1,"Dumping Trilinos/Stratimikos current settings to file: 'trilinos_current.xml'...\n");
    fflush(stdout); cout << flush;
    linearSolverBuilder.writeParamsFile(*lowsFactory,"./trilinos_current.xml");
  }
*/
  // solve the matrix
  RCP<Thyra::LinearOpWithSolveBase<double> > lows = Thyra::linearOpWithSolve(*lowsFactory, A); // create solver
  Thyra::solve(*lows, Thyra::NOTRANS, *b, &*x); // solve

  // r = b - A*x, final L2 norm
  {
    Epetra_Vector epetra_r(*epetra_b);
    Epetra_Vector epetra_A_x(epetra_A->OperatorRangeMap());
    epetra_A->Apply(*epetra_x,epetra_A_x);
    epetra_r.Update(-1.0,epetra_A_x,1.0);
    systemResidual=1./nrm_r;
    nrm_r=0.;
    epetra_r.Norm2(&nrm_r);
    systemResidual*=nrm_r;
  }

  // copying the solution into the data handle for the rhs: TODO make more efficient
  rhs = 0.0;
  for (int i=0; i<nrhs; i++) {
    //rhs[_lid[i]] = solVec->getValue(_gid[i]);
    rhs[_lid[i]] = (*epetra_x)[i];
  }

  stopTimer.stop();

  // print relative residual
  CFLog(INFO,"Linsolve time: " << stopTimer << "s\n");
  CFLog(INFO,"Rel. residual: " << systemResidual << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void StdParSolveSys::setup()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  // sizes
  int nbEqs = getMethodData().getNbSysEquations();
  int nbStates = states.size();

  // check (determining the number of locally owned unknowns --> TODO replace by the content of the cf_assertion)
  int nbMyStates = 0;
  for (int i = 0; i < nbStates; ++i) {
    if (states[i]->isParUpdatable()) {
      ++nbMyStates;
    }
  }
  const CFuint nbMyUnks = nbMyStates * nbEqs;
  cf_assert(nbMyUnks == getMethodData().getEpetraMap()->NumMyElements());

  // allocation
  _lid = new int[nbMyUnks];
  _gid = new int[nbMyUnks];

  // filling _lid and _gid
  int counter = 0;
  const LSSIdxMapping& idxMapping = getMethodData().getLocalToGlobalMapping();
  for (int i=0; i<nbStates; i++) {         // loop over the local state ids
    if (states[i]->isParUpdatable()) {
      for (int e=0; e<nbEqs; e++) {
        _lid[counter] = nbEqs*i + e;
        _gid[counter] = nbEqs*(idxMapping.getRowID(i)) + e; // TODO: check wether states[i]->getLocalID()==i
        counter++;
      }
    }
  }
  cf_assert(counter == nbMyUnks);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdParSolveSys::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD
