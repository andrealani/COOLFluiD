// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_PhysicalModelDummy_CouplingModelDummy_hh
#define COOLFluiD_PhysicalModelDummy_CouplingModelDummy_hh

//////////////////////////////////////////////////////////////////////////////

#include "PhysicalModelDummy/PhysicalModelDummy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PhysicalModelDummy {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a CouplingModelDummy
class CouplingModelDummy : public PhysicalModelDummy {

public:

  /// Constructor without arguments
  CouplingModelDummy(const std::string& name);

  /// Default destructor
  virtual ~CouplingModelDummy();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// @return name of the Physical Model type (std::string)
  virtual std::string getTypeName() const
  {
    return std::string("CouplingModelDummy");
  }

  /// Configures this object by complementing implementation in ConfigObject
  virtual void configure ( Config::ConfigArgs& args );

  /// @return vector storing IDs of the variables when sent
  const std::vector<CFuint>& getSendIDs() const {return m_sendIDs;}
  
  /// @return vector storing IDs of the variables when received
  const std::vector<CFuint>& getRecvIDs() const {return m_recvIDs;}
  
protected:  // data
  
  /// Vector storing IDs of the variables when sent
  std::vector<CFuint> m_sendIDs;
  
  /// Vector storing IDs of the variables when received
  std::vector<CFuint> m_recvIDs;
  
}; // end of class PhysicalModelDummy

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PhysicalModelDummy
}  // namespace COOLFluiD

#endif // COOLFluiD_PhysicalModelDummy_CouplingModelDummy_hh

