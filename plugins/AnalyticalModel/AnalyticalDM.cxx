// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "AnalyticalModel/AnalyticalDM.hh"
#include "AnalyticalModel/AnalyticalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AnalyticalModel {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<AnalyticalDM,
                            Framework::DomainModel,
                            AnalyticalModelModule, 1>
aAnalyticalDMProvider("AnalyticalDM");

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("ModelDim","Define the dimensionality of the model space.");
   options.addConfigOption< std::vector<std::string> >("FuncX","Definition of the x coordinates, one per topology region.");
   options.addConfigOption< std::vector<std::string> >("FuncY","Definition of the y coordinates, one per topology region.");
   options.addConfigOption< std::vector<std::string> >("FuncZ","Definition of the z coordinates, one per topology region.");
   options.addConfigOption< std::vector<std::string> >("DFuncX","Definition of the x 1st derivatives, dim-1 * per topology region.");
   options.addConfigOption< std::vector<std::string> >("DFuncY","Definition of the y 1st derivatives, dim-1 * per topology region.");
   options.addConfigOption< std::vector<std::string> >("DFuncZ","Definition of the z 1st derivatives, dim-1 * per topology region.");
}

//////////////////////////////////////////////////////////////////////////////

AnalyticalDM::AnalyticalDM(const std::string& name) :
  Framework::DomainModel(name),
  m_modeldim(0),
  m_pardim(0)
{
  addConfigOptionsTo(this);

  setParameter("ModelDim",&m_modeldim);

  setParameter("FuncX",&m_x_func_def);
  setParameter("FuncY",&m_y_func_def);
  setParameter("FuncZ",&m_z_func_def);

  setParameter("DFuncX",&m_x_dfunc_def);
  setParameter("DFuncY",&m_y_dfunc_def);
  setParameter("DFuncZ",&m_z_dfunc_def);
}

//////////////////////////////////////////////////////////////////////////////

AnalyticalDM::~AnalyticalDM()
{
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  Framework::DomainModel::configure(args);

  if ( m_modeldim != DIM_3D && m_modeldim != DIM_2D )
    throw Common::BadValueException (FromHere(),"AnalyticalModel : wrong model dimension " + Common::StringOps::to_str(m_modeldim));

  // parametric space is always one dimension less than the model space
  m_pardim = static_cast < CFuint > ( m_modeldim - 1 );

  // check the x functions at take this as the number of topologies
  m_nb_topo = m_x_func_def.size();

  // check that the other definitions also have same size
  if ( m_x_func_def.size() != m_y_func_def.size() )
    throw Common::BadValueException (FromHere(),"AnalyticalModel : number of definitions for model in Y differs form X");

  if ( m_x_func_def.size() != m_z_func_def.size() && m_modeldim == DIM_3D)
    throw Common::BadValueException (FromHere(),"AnalyticalModel : number of definitions for model in Z differs form X");

  std::vector< std::string > vars;

  if ( m_modeldim >= DIM_2D ) vars.push_back("u");
  if ( m_modeldim == DIM_3D ) vars.push_back("v");

  m_func.resize(m_nb_topo);
  for (CFuint i = 0; i < m_nb_topo; ++i )
  {
    CFLog ( INFO , "Parsing functions for topological region [" << i << "]\n");

    std::vector< std::string > topofunc;

    topofunc.push_back(m_x_func_def[i]);
    topofunc.push_back(m_y_func_def[i]);
    if ( m_modeldim == DIM_3D ) topofunc.push_back(m_z_func_def[i]);

    m_func[i].setVariables(vars);
    m_func[i].setFunctions(topofunc);
    try
    {
      m_func[i].parse();
    }
    catch (Common::ParserException& e)
    {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }

  if ( !m_x_dfunc_def.empty() )
    if ( m_x_dfunc_def.size() == m_pardim*m_x_func_def.size() )
    {
      if ( m_x_dfunc_def.size() != m_y_dfunc_def.size() )
        throw Common::BadValueException (FromHere(),"AnalyticalModel : number of 1st deriv definitions for model in Y differs form X");

      if ( m_x_dfunc_def.size() != m_z_dfunc_def.size() && m_modeldim == DIM_3D)
        throw Common::BadValueException (FromHere(),"AnalyticalModel : number of 1st deriv definitions for model in Z differs form X");

      m_dfunc.resize(m_pardim*m_nb_topo);
      for (CFuint k = 0; k < m_nb_topo; ++k )
      {
        CFLog ( INFO , "Parsing derivative functions for topological region [" << k << "]\n");

        for (CFuint j = 0; j < m_pardim; ++j ) {

        CFuint i = k+j;
        std::vector< std::string > topofunc;

        topofunc.push_back(m_x_dfunc_def[i]);
        topofunc.push_back(m_y_dfunc_def[i]);
        if ( m_modeldim == DIM_3D ) topofunc.push_back(m_z_dfunc_def[i]);

        m_dfunc[i].setVariables(vars);
        m_dfunc[i].setFunctions(topofunc);
        try
        {
          m_dfunc[i].parse();
        }
        catch (Common::ParserException& e)
        {
          CFout << e.what() << "\n";
          throw; // retrow the exception to signal the error to the user
        }

      } // loop paramdim
      } // loop topo
    }
    else
    {
      throw Common::BadValueException (FromHere(),"AnalyticalModel : number of derivatives differs from parameter space dim times number of functions");
    }
}

//////////////////////////////////////////////////////////////////////////////

Framework::DomainModel::TRidx AnalyticalDM::getNbTopoDefs () const
{
  cf_assert(m_nb_topo == m_func.size());
  return m_nb_topo;
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::computeParamCoord (const TRidx idx, const XVector& coord , PVector& pcoord) const
{
  cf_assert(idx < m_func.size());

  /// @todo only for curves
  if ( pcoord.size() != DIM_1D )
    throw Common::NotImplementedException (FromHere(),"AnalyticalDM::computeParamCoord() only supports pardim equals 1");

  /// @todo input coord vector should be projected onto the entity

  // coords that will vary
  XVector v_coord     (coord);
  XVector v_coord_eps (coord);
  XVector error       (coord);
  CFreal  fu;
  CFreal  fu_eps;
  CFreal  dgu;
  CFreal  gu;
  CFreal  gu_eps;
  PVector guess    (pcoord);
  PVector pu       (pcoord);
  PVector du       (pcoord);

// CF_DEBUG_OBJ(coord);

  const CFreal eps = 1E-6;
  const CFreal error_target = 1E-5;
  const CFuint max_iter = 300;

  // use the supplied pcoord value as first guess
  for ( CFuint k = 0; k < max_iter; ++k)
  {
    computeCoord(idx,pcoord,v_coord);
    error = coord - v_coord;
    fu = error.sqrNorm();

    // search is based on the first derivateive of the distance function being zero
    // not the distance itself because of bad convergence
    pu = pcoord + eps;
    computeCoord(idx,pu,v_coord_eps);
    error = coord - v_coord_eps;
    fu_eps = error.sqrNorm();

    gu = ( fu_eps - fu ) / eps ;

    if ( std::abs(gu) < error_target)
    {
      CFLogWarn ("coord [" << coord << "] converged at iter [" << k << "] to value  [" << pcoord << "] w computed coord [" << v_coord << "]\n");
      break;
    }

    pu = pcoord - eps;
    computeCoord(idx,pu,v_coord_eps);
    error = coord - v_coord_eps;
    fu_eps = error.sqrNorm();

    gu_eps = ( fu - fu_eps ) / eps ;

    dgu = ( gu_eps - gu ) / eps ;

    du = ( gu / dgu );

    pcoord += du;

    // make sure that the parametric values stay inside the parametric space
    if (pcoord[XX] < 0.) pcoord[XX] = 0.;
    if (pcoord[XX] > 1.) pcoord[XX] = 1.;
  }

  // warnd if didnt converge
  if (std::abs(gu) >= error_target)
  {
    CFLogWarn (" +++ AnalyticalDM was not able to compute parametric coord from coordinates\n" <<
               " +++ guess [" << guess << "] coord [" << coord << "] for TR  [" << idx << "]\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::computeCoord (const TRidx idx, const PVector& pcoord, XVector& coord) const
{
  cf_assert(idx < m_func.size());
  cf_assert(pcoord.size() == m_pardim);
  m_func[idx].evaluate(pcoord,coord);
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const
{
  cf_assert(idx < m_func.size());
  cf_assert(pcoord.size() == m_pardim);
  cf_assert(deriv1.size() == m_pardim);
  if (!m_dfunc.empty())
  {
    for (CFuint j = 0; j < m_pardim; ++j)
    {
      cf_assert(deriv1[j].size() == m_modeldim);
      m_dfunc[idx+j].evaluate(pcoord,deriv1[j]);
    }
  }
  else
  {
    throw Common::NotImplementedException (FromHere(),"AnalyticalDM::compute1stDeriv() does not implement numerical derivation yet");
  }
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const
{
  cf_assert(idx < m_func.size());
  throw Common::NotImplementedException (FromHere(),"AnalyticalDM::compute2ndDeriv");
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalDM::computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const
{
  cf_assert(idx < m_func.size());
  throw Common::NotImplementedException (FromHere(),"AnalyticalDM::computeAll");
}

//////////////////////////////////////////////////////////////////////////////

 } // namespace AnalyticalModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
