// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/StdUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  using namespace Framework;

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdUnSetup,
                      CFmeshReaderData,
                      CFmeshFileReaderPlugin> stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) : CFmeshReaderCom(name)
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::~StdUnSetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD
