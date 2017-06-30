// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFLUID_Common_SelfRegistPtr_hh
#define COOLFLUID_Common_SelfRegistPtr_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

#include "Common/ProviderBase.hh"
#include "Common/FailedCastException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Common {

//////////////////////////////////////////////////////////////////////////////

/// @brief Smart pointer for use to keep ownership of the self registrable objects created by the Factory's.
/// This class holds a pointer to an object that was created through
/// a Factory. This means it could possibly be allocated by a different
/// module, and whoever allocated it should also free it (to avoid undefined
/// behaviour --- mostly with overloading delete/new)
/// This smart pointer deletes the object through its provider
/// the moment the smart pointer is destroyed.
/// This pointer can hold a null pointer.
/// The pointer allows safe copying of objects thanks to an intrusive
/// reference counting technique, that requires the pointee object to revide from OwnedObject.
/// @see OwnedObject
/// @author Dries Kimpe
/// @author Andrea Lani
/// @author Tiago Quintino
template <typename TYPE>
class SelfRegistPtr {

public: // methods

  /// Default empty constructor
  SelfRegistPtr () : m_contents(CFNULL), m_provider(CFNULL) {}

  /// Constructor
  /// @param obj raw pointer to the object
  /// @param prov registered provider that created the object
  SelfRegistPtr (TYPE* obj, Common::ProviderBase* prov) :
    m_contents(obj), m_provider(prov)
  {
    cf_assert(m_contents != CFNULL);
#ifndef CF_HAVE_SINGLE_EXEC  
  cf_assert(m_provider != CFNULL);
#endif
    if (m_contents != CFNULL) { m_contents->addOwner(); }
  }

  /// Copy Constructor
  SelfRegistPtr (const SelfRegistPtr& source) :
    m_contents(source.m_contents),
    m_provider(source.m_provider)
  {
    if (m_contents != CFNULL) {
      m_contents->addOwner();
    }
  }

  /// Assignment operator
  /// @post releases ownership of the current contents
  const SelfRegistPtr<TYPE>& operator= (const SelfRegistPtr<TYPE>& source)
  {
    reset(source);
    return *this;
  }

  /// Destructor.
  /// The embedded object is freed by the correct provider.
  ~SelfRegistPtr() { release(); }

  /// Get the internal pointer
  TYPE* getPtr() const
  {
    cf_assert(m_contents != CFNULL);
    return m_contents;
  }

  /// Check if pointer is  CFNULL
  /// @return true if contents are CFNULL
  bool isNull() const {  return (m_contents == CFNULL); }

  /// Check if pointer is not CFNULL
  /// @return true if contents are not CFNULL
  bool isNotNull() const { return !isNull(); }

  /// Access to the object inside
  TYPE* operator-> () const
  {
    cf_assert (m_contents != CFNULL);
    return m_contents;
  }

  /// Access to the object inside
  TYPE & operator* () const
  {
    cf_assert (m_contents != CFNULL);
    return *m_contents;
  }

  /// Comparison operator
  bool operator== (const SelfRegistPtr<TYPE>& other) const
  {
    return (m_contents == other.m_contents && m_provider == other.m_provider);
  }

  /// Release the ownership of over the pointer contents
  void release()
  {
    if (m_contents != CFNULL)
    {
      m_contents->removeOwner();
      // ask m_provider to remove the pointed m_contents
      // object only if the later is not owned by anybody anymore
      if (m_contents->hasNoOwner())
      {
#ifndef CF_HAVE_SINGLE_EXEC
        cf_assert (m_provider != CFNULL);
#endif
        if (m_provider != CFNULL) {
          m_provider->freeInstance(m_contents);
          m_provider = CFNULL;
        } 
        m_contents = CFNULL;
      }
    }
  }

  ///  Reset the contents of the pointer to a another pointer
  void reset(const SelfRegistPtr<TYPE>& source)
  {
    release();
    m_contents = source.m_contents;
    m_provider = source.m_provider;
    if (m_contents != CFNULL) { m_contents->addOwner(); }
  }

  /// Dynamic casts the contents to the DTYPE
  /// @return another SelfRegistPtr of the DTYPE
  template <class DTYPE>
  SelfRegistPtr<DTYPE> d_castTo() const
  {
    DTYPE* rPtr = dynamic_cast<DTYPE*>(m_contents);
    if (rPtr == CFNULL)
    {
      std::string msg ("SelfRegistPtr failed dynamic cast from ");
      msg += DEMANGLED_TYPEID(TYPE);
      msg += " to ";
      msg += DEMANGLED_TYPEID(DTYPE);
      throw Common::FailedCastException (FromHere(),msg);
    }
    return SelfRegistPtr<DTYPE>(rPtr,m_provider);
  }

  /// Gets the provider that constructed the object
  /// @return pointer to the provider
  Common::ProviderBase* getProviderBase() { return m_provider; }

private: // data

  /// pointee object aggregated by this smart pointer
  TYPE* m_contents;

  /// provider object that created the m_contents
  Common::ProviderBase* m_provider;

}; // end SelfRegistPtr

//////////////////////////////////////////////////////////////////////////////

    } // namespace Common

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOFluiD_Common_SelfRegistPtr_hh
