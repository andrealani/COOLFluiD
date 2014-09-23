// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Trilinos/Trilinos.hh"


#include "TrilinosLSSData.hh"
#include "Framework/MethodCommandProvider.hh"

#include "ml_MultiLevelPreconditioner.h"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<TrilinosLSSData>, TrilinosLSSData, TrilinosModule> nullTrilinosLSSComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSSData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< int >("OutputLevel","Output level: 0=none, 5=extreme (above 3 its printing matrix).");
   options.addConfigOption< string >("OptionsXML","Xml file containing the options for trilinos.");
}

//////////////////////////////////////////////////////////////////////////////

TrilinosLSSData::TrilinosLSSData(SafePtr<std::valarray<bool> > maskArray,
				 CFuint& nbSysEquations,
				 Common::SafePtr<Framework::Method> owner) :
  LSSData(maskArray, nbSysEquations, owner),
  _comm(NULL),
  _map(NULL),
  _xVec(),
  _bVec(),
  _aMat(),
  _optionsxml("trilinosoptions.xml")
{
  addConfigOptionsTo(this);

  _outputlevel = 0;
  setParameter("OutputLevel",&_outputlevel);

  //_optionsxml = "trilinosoptions.xml";
  setParameter("OptionsXML",&_optionsxml);

  TrilinosOptions::setAllOptions();
}

//////////////////////////////////////////////////////////////////////////////

TrilinosLSSData::~TrilinosLSSData()
{
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSSData::configure ( Config::ConfigArgs& args )
{
  LSSData::configure(args);

  CFLog(INFO, "Trilinos Options:\n");
  CFLog(INFO, "OutputLevel = " << _outputlevel << "\n");
  CFLog(INFO, "OptionsXML = " << _optionsxml << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSSData::setup()
{
  LSSData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void TrilinosLSSData::unsetup()
{
  LSSData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

