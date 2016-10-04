// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MeshData.hh"
#include "Framework/DataHandleOutput.hh"
#include "Common/BadValueException.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("SocketNames","The names of the sockets to output based on states.");
  options.addConfigOption< std::vector<std::string> >("VariableNames","The names of the variables.");
  options.addConfigOption< std::vector<std::string> >("CCSocketNames","The names of the sockets to output cell-centered.");
  options.addConfigOption< std::vector<std::string> >("CCVariableNames","The names of the variables for cell-centered output.");
  options.addConfigOption< std::vector<CFuint> >     ("CCBlockSize","The block size for each socket. Default is one.");
  options.addConfigOption< std::vector<std::string> >("CCTrs","The trs with the cells to output the cell-centered data.");
  options.addConfigOption< bool >("isNodal","Tells if the variables to print are nodal.");
}

//////////////////////////////////////////////////////////////////////////////

DataHandleOutput::DataHandleOutput(const std::string& name) :
  NumericalStrategy(name),
  NamespaceMember(),
  socket_states("states"),
  socket_nodes("nodes"),
  m_cached_datahandles(false)
{
  addConfigOptionsTo(this);

  m_socket_names = std::vector<std::string> ();
  setParameter("SocketNames",&m_socket_names);

  m_varnames = std::vector<std::string> ();
  setParameter("VariableNames",&m_varnames);

  m_ccsocket_names = std::vector<std::string> ();
  setParameter("CCSocketNames",&m_ccsocket_names);

  m_ccvarnames = std::vector<std::string> ();
  setParameter("CCVariableNames",&m_ccvarnames);

  m_ccvartrs = std::vector<std::string> ();

  m_ccsockettrs = std::vector<std::string> ();
  setParameter("CCTrs",&m_ccsockettrs);

  m_ccblocksize = std::vector<CFuint> ();
  setParameter("CCBlockSize",&m_ccblocksize);
  
  m_isNodal = false;
  setParameter("isNodal",&m_isNodal);
}

//////////////////////////////////////////////////////////////////////////////

DataHandleOutput::~DataHandleOutput()
{
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::configure ( Config::ConfigArgs& args )
{
  NumericalStrategy::configure(args);

  for (std::vector<std::string>::const_iterator it = m_socket_names.begin();
       it != m_socket_names.end();
       ++it)
  {
    m_sockets.createSocketSink<CFreal>(*it);
  }

  for (std::vector<std::string>::const_iterator it = m_ccsocket_names.begin();
       it != m_ccsocket_names.end();
       ++it)
  {
    m_sockets.createSocketSink<CFreal>(*it);
  }

}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::setupStateVars()
{
  CFAUTOTRACE;
  
  // setup of the state based data output
  const CFuint nbDofs = (!m_isNodal) ? 
    socket_states.getDataHandle().size() : socket_nodes.getDataHandle().size();
  m_nbvars.clear();

  CFuint total_nbvars = 0;

  // if user did not provide the varnames we will construct
  // them based on the socket names
  bool fill_varnames = m_varnames.empty();

  // for each socket compute the number of variables per socket
  // assuming that it is a multiple of the number of states
  for (std::vector<std::string>::const_iterator it = m_socket_names.begin();
       it != m_socket_names.end();
       ++it)
  {
    Common::SafePtr<DataSocketSink<CFreal> > this_socket =
      m_sockets.getSocketSink<CFreal>(*it);

    DataHandle<CFreal> dh = this_socket->getDataHandle();
    
    CFuint nbvars = dh.size() / nbDofs;
    m_nbvars.push_back(nbvars);
    
    total_nbvars += nbvars;
    
    if ( fill_varnames )
    {
      if (nbvars == 1) {
        m_varnames.push_back( *it );
      }
      else {
        for ( CFuint iv = 0; iv < nbvars; ++iv) {
          m_varnames.push_back( *it + Common::StringOps::to_str(iv) );
	}
      }
    }
  }
  
  // check the total nb of variables matches the total nb of variable names
  if ( total_nbvars !=  m_varnames.size() ) {
    CFLog(ERROR, "total_nbvars [" << total_nbvars<< "] !=  m_varnames.size() [" << m_varnames.size()<< "] \n");
    throw BadValueException (FromHere(),"Number of variable names does not match number of variables");
  }
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::setupCCVars()
{
  CFAUTOTRACE;
  m_ccnbvars.clear();
  m_ccvartrs.clear();

  CFuint total_ccnbvars = 0;

  std::vector<SafePtr<TopologicalRegionSet> > trsList =
    MeshDataStack::getActive()->getTrsList();

  //  create the block size vector is user did not do it
  if ( m_ccblocksize.empty() )
  {
    m_ccblocksize.resize(m_ccsocket_names.size(),1);
  }

  // throw error if user gave wrong number of parameters
  if ( m_ccblocksize.size() != m_ccsocket_names.size() )
    throw BadValueException (FromHere(),"Vector with block sizes doesnt have the same size has the cell-centered sockets");

  // if user does not supply the TRS assume it to be named InnerCells
  if ( m_ccsockettrs.empty() )
  {
    m_ccsockettrs.resize(m_ccsocket_names.size(),"InnerCells");
  }

  if ( m_ccsockettrs.size() != m_ccsocket_names.size())
    throw BadValueException (FromHere(),"Must define one TRS to each cell-centered socket to output");

  // if user did not provide the varnames we will construct
  // them based on the socket names
  bool fill_ccvarnames = m_ccvarnames.empty();

  // for each socket compute the number of variables per socket
  // assuming that it is a multiple of the number of states
  cf_assert(m_ccsockettrs.size() == m_ccsocket_names.size());
  for ( CFuint i = 0; i < m_ccsocket_names.size(); ++i)
  {
    Common::SafePtr<DataSocketSink<CFreal> > this_socket =
      m_sockets.getSocketSink<CFreal>( m_ccsocket_names[i] );

    // this will check that the TRS exists
    SafePtr<TopologicalRegionSet> trs =
      MeshDataStack::getActive()->getTrs( m_ccsockettrs[i] );

    // lets get the element type description
    SafePtr<vector<ElementTypeData> > elementType =
      MeshDataStack::getActive()->getElementTypeData(trs->getName());

    // compute the number of cells in teh TRS
    CFuint nbCellsInTRS = 0;
    for (CFuint iType = 0; iType < elementType->size(); ++iType)
    {
      const CFuint nbCellsInType  = (*elementType)[iType].getNbElems();
      nbCellsInTRS += nbCellsInType;
    }

    DataHandle<CFreal> dh = this_socket->getDataHandle();

    const CFuint nbvars = dh.size() / ( nbCellsInTRS * m_ccblocksize[i] );
    m_ccnbvars.push_back(nbvars);

    total_ccnbvars += nbvars;

    // attach the TRS name to each of the variables independently
    for ( CFuint iv = 0; iv < nbvars; ++iv)
      m_ccvartrs.push_back(m_ccsockettrs[i]);

    if ( fill_ccvarnames )
    {
      if (nbvars == 1)
        m_ccvarnames.push_back( m_ccsocket_names[i] );
      else
        for ( CFuint iv = 0; iv < nbvars; ++iv)
          m_ccvarnames.push_back( m_ccsocket_names[i] + Common::StringOps::to_str(iv) );
    }
  }


  // check the total nb of variables matches the total nb of variable names
  if ( total_ccnbvars !=  m_ccvarnames.size() )
    throw BadValueException (FromHere(),"Number of cell-centered variable names does not match number of variables");

}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::setup()
{
  CFAUTOTRACE;
  NumericalStrategy::setup();
  setupStateVars();
  setupCCVars();
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::unsetup()
{
  // other stuff comes here

  NumericalStrategy::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::printVarNames(std::ostream& out) const
{
  for (std::vector<std::string>::const_iterator it = m_varnames.begin();
       it != m_varnames.end();
       ++it)
  {
    out << " \"" << *it << "\"";
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<std::string> DataHandleOutput::getVarNames () const
{
  return m_varnames;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<std::string> DataHandleOutput::getCCVarNames () const
{
  return m_ccvarnames;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<std::string> DataHandleOutput::getCCVarTrs () const
{
  return m_ccvartrs;
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::getDataHandles()
{
  m_datahandles.clear();
  m_ccdatahandles.clear();
  m_cached_datahandles = false;

  for (CFuint is = 0; is < m_socket_names.size(); ++is)
  {
    DataHandle<CFreal> dh =
      m_sockets.getSocketSink<CFreal>(m_socket_names[is])->getDataHandle();

    m_datahandles.push_back(dh);
  }

  for (CFuint is = 0; is < m_ccsocket_names.size(); ++is)
  {
    DataHandle<CFreal> dh =
      m_sockets.getSocketSink<CFreal>(m_ccsocket_names[is])->getDataHandle();

    m_ccdatahandles.push_back(dh);
  }

  m_cached_datahandles = true;
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::printStateData(std::ostream& out, CFuint state_id) const
{
  cf_assert(m_cached_datahandles);
  for (CFuint is = 0; is < m_socket_names.size(); ++is)
  {
    const DataHandle<CFreal>& dh = m_datahandles[is];
    CFLog(VERBOSE, state_id+1 << "/" << dh.size() << "\n");
    const CFuint dh_nbvars = m_nbvars[is];
    for (CFuint iv = 0; iv < dh_nbvars; ++iv)
    {
      out << dh(state_id,iv,dh_nbvars) << " ";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::fillStateData(CFreal* out, CFuint state_id, CFuint& counter, 
				     int iVar) const
{
  cf_assert(m_cached_datahandles);
  
  if (iVar < 0) {
    for (CFuint is = 0; is < m_socket_names.size(); ++is) {
      const DataHandle<CFreal>& dh = m_datahandles[is];
      CFLog(VERBOSE, state_id+1 << "/" << dh.size() << "\n");
      const CFuint dh_nbvars = m_nbvars[is];
      for (CFuint iv = 0; iv < dh_nbvars; ++iv) {
	out[counter++] = dh(state_id,iv,dh_nbvars);
      }
    }
  }
  else {
    CFuint countVar = 0;
    for (CFuint is = 0; is < m_socket_names.size(); ++is) {
      const DataHandle<CFreal>& dh = m_datahandles[is];
      CFLog(VERBOSE, state_id+1 << "/" << dh.size() << "\n");
      const CFuint dh_nbvars = m_nbvars[is];
      for (CFuint iv = 0; iv < dh_nbvars; ++iv, ++countVar) {
	if (countVar == iVar) {
	  out[counter++] = dh(state_id,iv,dh_nbvars);
	  return;
	}
      }
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

void DataHandleOutput::printCCData(std::ostream& out, CFuint cell_id) const
{
  cf_assert(m_cached_datahandles);
  for (CFuint is = 0; is < m_ccsocket_names.size(); ++is)
  {
    const DataHandle<CFreal>& dh = m_ccdatahandles[is];
    const CFuint dh_nbvars = m_ccnbvars[is];
    for (CFuint iv = 0; iv < dh_nbvars; ++iv)
    {
      out << dh(cell_id,iv,dh_nbvars) << " ";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
    
void DataHandleOutput::fillStateDataCC(CFreal* out, CFuint cell_id, CFuint& counter, 
				       int iVar) const
{ 
  cf_assert(m_cached_datahandles);
  if (iVar < 0) {
    for (CFuint is = 0; is < m_ccsocket_names.size(); ++is) {
      const DataHandle<CFreal>& dh = m_ccdatahandles[is];
      const CFuint dh_nbvars = m_ccnbvars[is];
      for (CFuint iv = 0; iv < dh_nbvars; ++iv) {
	out[counter++] = dh(cell_id,iv,dh_nbvars);
      }
    }
  }
  else {
    CFuint countVar = 0;
    for (CFuint is = 0; is < m_ccsocket_names.size(); ++is) {
      const DataHandle<CFreal>& dh = m_ccdatahandles[is];
      const CFuint dh_nbvars = m_ccnbvars[is];
      for (CFuint iv = 0; iv < dh_nbvars; ++iv, ++countVar) {
	if (countVar == iVar) {
	  out[counter++] = dh(cell_id,iv,dh_nbvars);
	  return;
	}
      }
    }
  }
}
  
//////////////////////////////////////////////////////////////////////////////
    
DataHandleOutput::DataHandleInfo
DataHandleOutput::getStateData(CFuint var_id) const
{
  cf_assert(m_cached_datahandles);

  CFuint dh_var_id = 0;
  CFuint varcount  = 0;
  for ( ; dh_var_id < m_nbvars.size(); ++dh_var_id )
  {
    varcount += m_nbvars[dh_var_id];
    if ( varcount > var_id )
      break;
  }
  return make_Trio( m_nbvars[dh_var_id] - ( varcount - var_id ),
                    m_nbvars[dh_var_id],
                    m_datahandles[ dh_var_id ]);
}

//////////////////////////////////////////////////////////////////////////////

DataHandleOutput::DataHandleInfo
DataHandleOutput::getCCData(CFuint var_id) const
{
  cf_assert(m_cached_datahandles);

  CFuint dh_var_id = 0;
  CFuint varcount  = 0;
  for ( ; dh_var_id < m_ccnbvars.size(); ++dh_var_id )
  {
    varcount += m_ccnbvars[dh_var_id];
    if ( varcount > var_id )
      break;
  }
  return Common::make_Trio( m_ccnbvars[dh_var_id] - ( varcount - var_id ),
                           m_ccnbvars[dh_var_id],
                           m_ccdatahandles[ dh_var_id ]);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
DataHandleOutput::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
