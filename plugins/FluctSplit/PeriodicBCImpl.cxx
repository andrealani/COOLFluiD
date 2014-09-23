#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Common/BadValueException.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/PeriodicBCImpl.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PeriodicBCImpl,
                      FluctuationSplitData,
                      FluctSplitModule>
PeriodicBCImplProvider("PeriodicBCImpl");

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImpl::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCImpl::PeriodicBCImpl(const std::string& name) : PeriodicBC(name),
	  socket_bStatesNeighbors("bStatesNeighbors"),
  _jacobElem(4,4),
  _jacobElema(4,4),
  _in(4),
  _irc(4),
  _ira(4),
  _block(4,4)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

PeriodicBCImpl::~PeriodicBCImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
PeriodicBCImpl::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = PeriodicBC::needsSockets();
  result.push_back(&socket_bStatesNeighbors);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImpl::setup()
{

  PeriodicBC::setup();
}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImpl::executeOnTrs()
{
  CFAUTOTRACE;

  Common::SafePtr<LinearSystemSolver> lss =
  getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  jacobMatrix->finalAssembly();

  DataHandle < Framework::State*, Framework::GLOBAL > states      = socket_states.getDataHandle();
  DataHandle<bool>    isUpdated   = socket_isUpdated.getDataHandle();
  DataHandle<CFreal>  rhs         = socket_rhs.getDataHandle();
  DataHandle<CFreal>  updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle< std::valarray<Framework::State*> > bStatesNeighbors = socket_bStatesNeighbors.getDataHandle();

  Common::SafePtr< vector<CFuint> > const applied_trs_statesIdx = getCurrentTRS()->getStatesInTrs();


  const CFuint nbEqs = m_tmp_rhs.size();

  const CFreal cfl = getMethodData().getCFL()->getCFLValue();
  const CFreal invcfl = 1.0 / cfl;

  for (CFuint is = 0; is < applied_trs_statesIdx->size(); ++is)
  {
//    CF_DEBUG_OBJ(is);
    const CFuint applied_sid = (*applied_trs_statesIdx)[is];
    const CFuint coupled_sid = m_match_states_idx[is];
// unused // State& applied_state = *(states[applied_sid]);
// unused // State& coupled_state = *(states[coupled_sid]);

    // residuals

    for (CFuint j = 0; j < nbEqs; ++j)
        m_tmp_rhs[j] =  rhs(applied_sid,j,nbEqs);

    for (CFuint j = 0; j < nbEqs; ++j)
        m_tmp_rhs[j] += rhs(coupled_sid,j,nbEqs);

    for (CFuint j = 0; j < nbEqs; ++j)
    {
      rhs(coupled_sid,j,nbEqs) = m_tmp_rhs[j];
      rhs(applied_sid,j,nbEqs) = m_tmp_rhs[j];
    }

    /// @todo this BC is not working and needs to be finished

    // update coefficient

    const CFreal tmp_upcoeff = updateCoeff[applied_sid] + updateCoeff[coupled_sid];
    updateCoeff[applied_sid] = tmp_upcoeff;
    updateCoeff[coupled_sid] = tmp_upcoeff;



    CFuint nbVarsInvolved = 4;
    CFuint nbEntries = bStatesNeighbors[coupled_sid].size();

    cf_assert(nbEntries > 0);
    
    // set the row idsfor
     for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
      _irc[iVar] = coupled_sid*nbEqs + iVar;
      _ira[iVar] = applied_sid*nbEqs + iVar;
     }
     CFuint entryID;
     CFuint nStart;

     for (CFuint i = 0; i < nbEntries; ++i)
     {
       entryID = bStatesNeighbors[coupled_sid][i]->getLocalID();
       if (entryID != coupled_sid){
       nStart = nbEqs*entryID;


     // set the column ids
    for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
         _in[iVar] = nStart + iVar;    
    }
	  // get the elements of the jacobian matrix
    	    jacobMatrix->getValues(nbVarsInvolved,
    	  			 &_irc[0],
    	  			 nbEqs,
    	  			 &_in[0],
    	  			 &_jacobElem[0]);

     	    jacobMatrix->finalAssembly();
	   

	   // get the elements of the jacobian matrix
    	    jacobMatrix->addValues(nbVarsInvolved,
    	  			 &_ira[0],
    	  			 nbEqs,
    	  			 &_in[0],
    	  			 &_jacobElem[0]);
	    jacobMatrix->finalAssembly();
       }
     }

     nbEntries = bStatesNeighbors[applied_sid].size();

    // cf_assert(nbEntries > 0);
    
     for (CFuint i = 0; i < nbEntries; ++i)
     {
     entryID = bStatesNeighbors[applied_sid][i]->getLocalID();
      if (entryID != applied_sid){
     nStart = nbEqs*entryID;


    // set the column ids
     for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
       _in[iVar] = nStart + iVar;
	    
     }
    // get the elements of the jacobian matrix
  	    jacobMatrix->getValues(nbVarsInvolved,
   	  			 &_ira[0],
   	  			 nbEqs,
	  			 &_in[0],
	  			 &_jacobElem[0]);
	    jacobMatrix->finalAssembly();

	   // get the elements of the jacobian matrix
    	    jacobMatrix->addValues(nbVarsInvolved,
    	  			 &_irc[0],
    	  			 nbEqs,
    	  			 &_in[0],
    	  			 &_jacobElem[0]);
    	   jacobMatrix->finalAssembly();
    }
     }

     jacobMatrix->getValues(nbVarsInvolved,
    	  			 &_irc[0],
    	  			 nbEqs,
    	  			 &_irc[0],
    	  			 &_jacobElem[0]);


     jacobMatrix->getValues(nbVarsInvolved,
    	  			 &_ira[0],
    	  			 nbEqs,
    			    &_ira[0],
    				 &_jacobElema[0]);

    jacobMatrix->finalAssembly();

     for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
    	for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
    
    	  _block(iVar, jVar) = _jacobElem(iVar, jVar) + _jacobElema(iVar,jVar);
           }
        }

      jacobMatrix->setValues(nbVarsInvolved,
			     &_irc[0],
			     nbEqs,
			     &_irc[0],
			     &_block[0]);

     
     jacobMatrix->setValues(nbVarsInvolved,
			    &_ira[0],
			    nbEqs,
			    &_ira[0],
			    &_block[0]);

     jacobMatrix->finalAssembly();
     

  } // end loop states in trs

}

//////////////////////////////////////////////////////////////////////////////

void PeriodicBCImpl::configure ( Config::ConfigArgs& args )
{
  PeriodicBC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
