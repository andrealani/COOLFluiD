// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_QualifiedName_hh
#define COOLFluiD_Framework_QualifiedName_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NamedObject.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class describes the full qualified name of an entity within the correct Namespace.
///  @author Tiago Quintino
class Framework_API QualifiedName : public Common::NamedObject
{
public:

  /// Constructor accepts the namespace and name
  /// @param namesp Namespace for full qualified name
  /// @param name Name of entity within the Namespace
  QualifiedName(const std::string& namesp, const std::string& name);

  /// Destructor
  ~QualifiedName();

  /// @return the full qualified name as a string
  const std::string& str() const { return m_qname; }

  /// Accessor to the Namespace string
  /// @return CFstring with namespace
  std::string getNamespace() const
  {
    return m_namespace;
  }

  /// Operator for comparison between full qualified names
  bool operator== (const QualifiedName& qnr) const;
  
  /// @return a string with the separator for namespaces
  static CFchar separator ();

private: // data

  /// storage of the namespace
  std::string m_namespace;

  /// full qualified name
  std::string m_qname;

}; // end class QualifiedName

/// Class to compare between two QualifiedName's
struct QNameComp
{
  bool operator()(const QualifiedName& q1, const QualifiedName& q2) const
  {
//   return strcmp(s1, s2) < 0;
    return ( q1.str() < q2.str() );
  }
};

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_QualifiedName_hh
