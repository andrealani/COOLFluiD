// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Petsc/PetscHeaders.hh" // must come before any header

#include "Petsc/PetscOptions.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Petsc {

//////////////////////////////////////////////////////////////////////////////

std::map<std::string,
         PCType,
         std::less<std::string> > PetscOptions::_pcType;

std::map<std::string,
         KSPType,
         std::less<std::string> > PetscOptions::_kspType;

std::map<std::string,
         MatOrderingType,
         std::less<std::string> > PetscOptions::_matOrderType;

//////////////////////////////////////////////////////////////////////////////

void PetscOptions::setAllOptions()
{
  setPCTypes();
  setKSPTypes();
  setMatOrderTypes();
}

//////////////////////////////////////////////////////////////////////////////

void PetscOptions::setPCTypes()
{
  _pcType["PCJACOBI"]    = PCJACOBI;
  _pcType["PCBJACOBI"]   = PCBJACOBI;
  _pcType["PCSOR"]       = PCSOR;
  _pcType["PCEISENSTAT"] = PCEISENSTAT;
  _pcType["PCICC"]       = PCICC;
  _pcType["PCILU"]       = PCILU;
  _pcType["PCASM"]       = PCASM;
  _pcType["PCCOMPOSITE"] = PCCOMPOSITE;
  _pcType["PCLU"]        = PCLU;
  _pcType["PCCHOLESKY"]  = PCCHOLESKY;
  _pcType["PCNONE"]      = PCNONE;
  _pcType["PCSHELL"]     = PCSHELL;
  _pcType["PCSACUSP"]    = PCSACUSP;
  _pcType["PCSACUSPPOLY"] = PCSACUSPPOLY;
  _pcType["PCBICGSTABCUSP"] = PCBICGSTABCUSP;
  _pcType["PCGASM"] = PCGASM;
  _pcType["PCKSP"] = PCKSP;
  _pcType["PCGAMG"] = PCGAMG;
#ifdef CF_HAVE_VIENNACL
  _pcType["PCSAVIENNACL"] = PCSAVIENNACL;
  _pcType["PCCHOWILUVIENNACL"] = PCCHOWILUVIENNACL;
#endif
}

//////////////////////////////////////////////////////////////////////////////

void PetscOptions::setMatOrderTypes()
{
  _matOrderType["MATORDERING_NATURAL"]   = MATORDERINGNATURAL;
  _matOrderType["MATORDERING_ND"]        = MATORDERINGND;
  _matOrderType["MATORDERING_1WD"]       = MATORDERING1WD;
  _matOrderType["MATORDERING_RCM"]       = MATORDERINGRCM;
  _matOrderType["MATORDERING_QMD"]       = MATORDERINGQMD;
  _matOrderType["MATORDERING_ROWLENGTH"] = MATORDERINGROWLENGTH;
  _matOrderType["MATORDERING_AMD"]       = MATORDERINGAMD;
}

//////////////////////////////////////////////////////////////////////////////

void PetscOptions::setKSPTypes()
{
  _kspType["KSPRICHARDSON"] = KSPRICHARDSON;
  // _kspType["KSPCHEBYCHEV"]  = KSPCHEBYCHEV;
  _kspType["KSPCG"]         = KSPCG;
  _kspType["KSPBICG"]       = KSPBICG;
  _kspType["KSPGMRES"]      = KSPGMRES;
  _kspType["KSPBCGS"]       = KSPBCGS;
  _kspType["KSPCGS"]        = KSPCGS;
  _kspType["KSPTFQMR"]      = KSPTFQMR;
  _kspType["KSPTCQMR"]      = KSPTCQMR;
  _kspType["KSPCR"]         = KSPCR;
  _kspType["KSPLSQR"]       = KSPLSQR;
  _kspType["KSPPREONLY"]    = KSPPREONLY;
  _kspType["KSPLGMRES"]     = KSPLGMRES;
  _kspType["KSPFGMRES"]     = KSPFGMRES;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Petsc

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
