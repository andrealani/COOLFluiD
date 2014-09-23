// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"

#include "CellQualityRemeshCondition.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CellQualityRemeshCondition,
                           SimpleMeshAdapterData,
                   SimpleGlobalMeshAdapterModule>
cellQualityRemeshConditionProvider("CellQualityRemeshCondition");

//////////////////////////////////////////////////////////////////////////////

void CellQualityRemeshCondition::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("MinValue","Minimum quality acceptable.");
}

//////////////////////////////////////////////////////////////////////////////

CellQualityRemeshCondition::CellQualityRemeshCondition(const std::string& name)  :
  SimpleMeshAdapterCom(name),
  socket_qualityCell("qualityCell")
{
   addConfigOptionsTo(this);

   _minValue = 1.0;
   setParameter("MinValue",&_minValue);

}

//////////////////////////////////////////////////////////////////////////////

void CellQualityRemeshCondition::execute()
{
  CFAUTOTRACE;

  DataHandle<CFreal> qualityCell = socket_qualityCell.getDataHandle();

  CFreal minQuality = MathTools::MathConsts::CFrealMax();
  const CFuint nbCells = qualityCell.size();

  for(CFuint iCell = 0; iCell < nbCells; iCell++)
  {
    if(qualityCell[iCell] < minQuality) minQuality = qualityCell[iCell];
  }

  if(minQuality < _minValue) getMethodData().setNeedRemeshing(true);
  else getMethodData().setNeedRemeshing(false);


}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CellQualityRemeshCondition::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_qualityCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
