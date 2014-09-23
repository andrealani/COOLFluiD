#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/LSSVector.hh"
#include "MathTools/RealVector.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/DirichletBC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< DirichletBC,UFEMSolverData,UFEMPlugin > DirichletBCProvider("DirichletBC");

//////////////////////////////////////////////////////////////////////////////

void DirichletBC::defineConfigOptions(Config::OptionList& options)
{
   CFAUTOTRACE;

   options.addConfigOption< std::string >( "Symmetry","Keep matrix symmetry by 'AdjustColumn' or 'ScaleDiagonal' methods or not (default is 'No')." );
   options.addConfigOption< CFreal >( "ScaleDiagonal","ScaleDiagonal symmetry method coefficient (default 1.e20)." );
   options.addConfigOption< bool >( "Implicit", "Apply the BC implicitly (solving for deltaFI)? (default is true)" );
}

//////////////////////////////////////////////////////////////////////////////

DirichletBC::DirichletBC(const std::string& name) :
  BaseBC(name),
  m_symmetryStr(""),
  socket_rhs("rhs"),
  socket_isUpdated("isUpdated"),
  socket_bStatesNeighbors("bStatesNeighbors")
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_symmetryStr = "No";
  m_scale = 1.e20;
  m_isImplicit = true;

  setParameter( "Symmetry",      &m_symmetryStr);
  setParameter( "ScaleDiagonal", &m_scale);
  setParameter( "Implicit",      &m_isImplicit);
}

//////////////////////////////////////////////////////////////////////////////

void DirichletBC::setup()
{
  CFAUTOTRACE;
  BaseBC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void DirichletBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  BaseBC::configure(args);

  CFLog(INFO, getClassName() << ": Symmetry: "       << m_symmetryStr << "\n");
  CFLog(INFO, getClassName() << ": ScaleDiagonal: "  << m_scale       << "\n");
  CFLog(INFO, getClassName() << ": Implicit: "       << m_isImplicit  << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void DirichletBC::executeOnTrs()
{
  CFAUTOTRACE;
  BaseBC::executeOnTrs();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  DataHandle< bool > isUpdated = socket_isUpdated.getDataHandle();
  DataHandle< State*, GLOBAL > states = socket_states.getDataHandle();
  DataHandle< std::valarray< State* > > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  Common::SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "DirichletBC::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  // get system matrix and global index mapping
  const Common::SafePtr< LinearSystemSolver > lss =  getMethodData().getLinearSystemSolver()[0];
  Common::SafePtr< LSSMatrix > sysMat = lss->getMatrix();
  const LSSIdxMapping& idxMapping = lss->getLocalToGlobalMapping();

  // PhysicalModel properties and auxiliary variables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const bool isAdjust = (m_symmetryStr=="AdjustColumn");
  const bool isScale  = (m_symmetryStr=="ScaleDiagonal");
  if (!isScale) m_scale = 1.;
  CFreal implicit = (m_isImplicit? 1.:0.);

  // this should be an intermediate lightweight assembly! it is needed because
  // here you SET values while elsewhere you ADD values
  sysMat->flushAssembly();
  sysMat->finalAssembly();

  // cycle all the states in the TRS
  Common::SafePtr< std::vector< CFuint > > trsStates = trs->getStatesInTrs();
  std::vector< CFuint >::iterator itd;
  CFuint istate=0;
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd, istate++) {
    const CFuint nLocalID = *itd;
    const State *currState = states[*itd];

    //if (isUpdated[nLocalID]) continue;
    //for (CFuint iEq=0; iEq<nbEqs; ++iEq) updated[iEq]=isUpdated[]
    if (currState->isParUpdatable()) {

      // global position of node, which must have at least one neighbour
      const CFuint nGlobalID = idxMapping.getColID(nLocalID)*nbEqs;
      const CFuint nbNeigh = bStatesNeighbors[nLocalID].size();

      //cf_assert(nbNeigh>0);

      // evaluate boundary condition at space, time and states,
      // and calculate the boundary condition enforced value
      computeStateValuesDirichletBC(currState);

      // loop over equations
      for(CFuint iEq=0; iEq<nbEqs; ++iEq) if ((!isUpdated[nLocalID*nbEqs+iEq])&&(m_applyFlags[istate][iEq])) {

        // setting value to be set, depending on implicit it is either Unknown or deltaUnknown
        const CFreal applyval = m_applyVars[iEq] - implicit*((*currState)[iEq]);

        // zeroing line if not asked by the options elsehow
        if (!isScale){

          // zeroing line if scaling wasnt applied
/*
          for (CFuint j=0; j<nbNeigh; ++j) {
            const CFuint jGlobalID = idxMapping.getColID(bStatesNeighbors[nLocalID][j]->getLocalID())*nbEqs;
            for (CFuint jEq=0; jEq<nbEqs; ++jEq) sysMat->setValue(nGlobalID+iEq,jGlobalID+jEq, 0.);
          }
/*/
          sysMat->finalAssembly();
          sysMat->setRow(nGlobalID+iEq,0.,0.);
/**/

          // zeroing column by summing to rhs, to ensure symmetry of the matrix
          if (isAdjust){
            for (CFuint j=0; j<nbNeigh; ++j) {
              const CFuint nbLocalID=bStatesNeighbors[nLocalID][j]->getLocalID();
              if (bStatesNeighbors[nLocalID][j]->isParUpdatable()) {
                const CFuint jGlobalID = idxMapping.getColID(bStatesNeighbors[nLocalID][j]->getLocalID())*nbEqs;
                for (CFuint jEq=0; jEq<nbEqs; ++jEq){
                  CFreal ijval=0.;
                  sysMat->finalAssembly();
                  sysMat->getValue(jGlobalID+jEq,nGlobalID+iEq,ijval);
                  sysMat->setValue(jGlobalID+jEq,nGlobalID+iEq, 0.);
                  rhs(nbLocalID,jEq,nbEqs) -= ijval * applyval;
                }
              }
            }
          } // for isAdjust

        } // for isScale

        // set rhs and system matrix diagonal terms (scaled)
        sysMat->setValue(nGlobalID+iEq,nGlobalID+iEq, m_scale);
        rhs(nLocalID,iEq,nbEqs) = applyval * m_scale;

        // flagging is important!
        isUpdated[nLocalID*nbEqs+iEq] = true;

      }

    } // isParUpdatable?

  } // cycle all the states in the TRS

  // this should be an intermediate lightweight assembly! it is needed because
  // here you SET values while elsewhere you ADD values
  sysMat->flushAssembly();

}

//////////////////////////////////////////////////////////////////////////////

void DirichletBC::computeStateValuesDirichletBC(const Framework::State* currState)
{
  computeStateValuesBaseBC(currState);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > DirichletBC::needsSockets()
{
  CFAUTOTRACE;

  std::vector<Common::SafePtr<BaseDataSocketSink> > result=BaseBC::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_bStatesNeighbors);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

