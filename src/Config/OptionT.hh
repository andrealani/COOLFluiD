// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_OptionT_hh
#define COOLFluiD_Config_OptionT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"

#include "Config/Option.hh"
#include "Config/BadMatchException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

  template < typename TYPE > class OptionT; // forward declaration

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
class Policy {
public:

    typedef TYPE value_type;

    static void setValue(const std::string& value, OptionT< TYPE >* opt )
    {
      *(opt->m_value) = Common::StringOps::from_str < TYPE > ( value );
    }

    static std::string description (const OptionT< TYPE >* const opt)
    {
      std::ostringstream os;
      os << opt->description() << " (default = " << opt->defaultValue() << ")";
      return os.str();
    }

    static std::string getDefaultValueAsString(const OptionT< TYPE >* const opt)
    {
      return Common::StringOps::to_str < TYPE > ( opt->m_defvalue );
    }

    static std::string getValueAsString(const OptionT< TYPE >* const opt)
    {
      return Common::StringOps::to_str < TYPE >(*(opt->m_value));
    }

    static std::string debugInfo(const OptionT< TYPE >* const opt)
    {
      std::ostringstream os;
      os << "{ " << opt->getLabel() << ", " << opt->getHierarchyOwner() << ", " << opt->description() << ", " << opt->value() << " }";
      return os.str();
    }

    static std::string getXML(const OptionT< TYPE >* const opt)
    {
      std::ostringstream os;
      os << "<" << opt->getLabel() << "\t";
      os << "tree=\"" << "option" << "\" ";
      os << "type=\"" << DEMANGLED_TYPEID(TYPE) << "\" ";
      os << "mode=\"" << ( opt->isBasic() ? "basic" : "advanced" ) << "\" ";
      os << "dynamic=\"" << ( opt->isDynamic() ? "dynamic" : "static" ) << "\" ";
      os << "description=\"" << opt->description() << "\" ";
      os << ">\n";
      os << opt->getValueAsString() << "\n";
      os << "</" << opt->getLabel() << ">\n";
      return os.str();
    }

    static CFuint setCommandLineValues(
      const std::string& arg,
      CFuint position,
      int argc,
      char** argv,
      OptionT< TYPE >* opt)
    {
      int nArgs = 1;
      ++position;

      if (position < CFuint(argc)) {
          // there needs to be a following argument, which is the value
          std::string argstr(argv[position]);
          opt->checkArgument(argstr);
          opt->setValue(argstr);
          position++;         // and now consume the value
          nArgs++;
      }
      else {
            throw BadMatchException (FromHere(),"argument expected for " + opt->getLabel());
      }
      return nArgs;
    }

};

//////////////////////////////////////////////////////////////////////////////

template < typename TYPE >
class Policy < std::vector<TYPE> > {
public:

    typedef TYPE value_type;

private:

    typedef std::vector<TYPE> Container;

public:

    static void setValue(const std::string& value, OptionT<  std::vector < TYPE > >* opt)
    {
      std::vector<std::string> words;
      words = Common::StringOps::getWords(value);
      Policy < std::vector < TYPE > >::replaceValues(words,opt);
    }

    static std::string description (const OptionT<  std::vector<TYPE> >* const opt)
    {
      std::ostringstream os;
      os << opt->description();
      if (opt->m_defvalue.size() > 0) {
          os << " (default = ";
          for (typename Container::const_iterator it = opt->m_defvalue.begin(); it != opt->m_defvalue.end(); ++it) {
              if (it != opt->m_defvalue.begin()) {
                  os << ", ";
              }
              os << *it;
          }
          os << ")";
      }
      return os.str();
    }

    static std::string getDefaultValueAsString(const OptionT<  std::vector < TYPE > >* const opt)
    {
      std::string str;
      for (typename Container::const_iterator it = opt->m_defvalue.begin(); it != opt->m_defvalue.end(); ++it) {
          if (it != opt->m_defvalue.begin()) { str += " "; }
          str += Common::StringOps::to_str(*it);
      }
      return str;
    }

    static std::string getValueAsString(const OptionT<  std::vector < TYPE > >* const opt)
    {
      std::string str;
      for (typename Container::const_iterator it = opt->m_value->begin(); it != opt->m_value->end(); ++it) {
          if (it != opt->m_value->begin()) { str += " "; }
          str += Common::StringOps::to_str(*it);
      }
      return str;
    }

    static std::string debugInfo (const OptionT<  std::vector < TYPE > >* const opt)
    {
      std::ostringstream os;
      os << "{ " << opt->getLabel() << ", " << opt->getHierarchyOwner() << ", " << opt->description() << ", [";
      typename Container::const_iterator it   = opt->m_value->begin();
      typename Container::const_iterator stop = opt->m_value->end();
      for (; it != stop; ++it, os << ";")
        os << *it;
      os << "] }";
      return os.str();
    }

    static std::string getXML(const OptionT<  std::vector < TYPE > >* const opt)
    {
      std::ostringstream os;
      os << "<" << opt->getLabel() << "\t";
      os << "tree=\"" << "option" << "\" ";
      os << "type=\"" << DEMANGLED_TYPEID(TYPE) << "\" ";
      os << "mode=\"" << ( opt->isBasic() ? "basic" : "advanced" ) << "\" ";
      os << "dynamic=\"" << ( opt->isDynamic() ? "dynamic" : "static" ) << "\" ";
      os << "description=\"" << opt->description() << "\" ";
      os << ">\n";
      opt->getValueAsString();
      os << "</" << opt->getLabel() << ">\n";
      return os.str();
    }


    static CFuint setCommandLineValues(
      const std::string& arg,
      CFuint position,
      int argc,
      char** argv,
      OptionT<  std::vector < TYPE > >* opt)
    {
        int nArgs = 1;
        std::string endtag = Policy< std::vector<TYPE> >::endTag(arg);

        // collect the values until we get the end tag
        std::vector<std::string> values;
        for (int p = position + 1; p < argc && argv[p] != endtag; ++p, ++nArgs) {
            values.push_back(std::string(argv[p]));
        }

        Policy< std::vector<TYPE> >::replaceValues(values,opt);

        // advance through the end of list tag
        ++nArgs;
        return nArgs;
    }

private:

    static void replaceValues(const std::vector<std::string>& values, OptionT<  std::vector<TYPE> >* opt)
    {
     //!!! add a recurse guard for linked options?

      // empty the current values, since we'll be pushing all the new ones
      // onto the collection
      opt->m_value->erase(opt->m_value->begin(), opt->m_value->end());

      std::vector<std::string>::const_iterator vit   = values.begin();
      std::vector<std::string>::const_iterator vstop = values.end();
      for (; vit != vstop; ++vit) {
          std::string val = *vit;
          Policy< std::vector<TYPE> >::addValue(val,opt);
      }

      // any linked value gets the same treatment
      std::vector<Option*>::iterator lit   = opt->getLinkedOptions().begin();
      std::vector<Option*>::iterator lstop = opt->getLinkedOptions().end();
      for (; lit != lstop; ++lit) {
          Option* lopt = *lit;
          // the linked option should be of the same type as this one:
          OptionT< std::vector<TYPE> >* other = dynamic_cast<OptionT< std::vector < TYPE > >*>(lopt);
          Policy< std::vector<TYPE> >::replaceValues(values,other);
      }
    }

  static std::string endTag(const std::string& tag)
  {
    std::string sw("start-of-");
    std::string t;
    if(Common::StringOps::startsWith(tag, sw)) {
      t = tag.substr(sw.length());
    }
    else {
      t = tag;
    }
    return std::string("--end-of-") + t;
  }

  static void addValue(const std::string& value, OptionT<  std::vector<TYPE> >* opt)
  {
    opt->m_value->push_back(Common::StringOps::from_str<typename std::vector<TYPE>::value_type>(value));
  }

};

//////////////////////////////////////////////////////////////////////////////

///  Handles generic data types as a user-specified option
///  All of the C++ built-in data types are supported,
///  as well as the STL string class.
///  This class was shamelessly copied and adapted form Yagol.
///  All credits go to Jeff Pace <jpace@incava.org>. Long live OpenSource!!!
///  @author Tiago Quintino
template <class Type>
class OptionT : public Option {
public:

  typedef Policy<Type> StoragePolicy;
  friend class Policy<Type>;

public:

  /// Constructor
  OptionT(const std::string& name, const std::string& description);

  /// Copy constructor
  OptionT(const OptionT<Type>& rhs);

  /// Assignement operator
  OptionT<Type>& operator=(const OptionT<Type>& rhs);

  /// Virtual destructor.
  virtual ~OptionT();

  /// Links an external owned value to the Option
  /// @param value pointer to the parameter to link
  virtual void linkToValue(void* value, const std::string& typestr);

  /// Sets the default value
  /// @param defvalue pointer to the value to set
  virtual void setDefaultValue(void* defvalue, const std::string& typestr);

  /// Sets this option value from a std::string
  /// @pre label has already been checked to ensure that it matches the value's tag
  /// @param value string with value to be converted to actual type
  virtual void setValue(const std::string& value);

  /// Returns the value of the tied object.
  virtual Type value() const;

  /// Gets the option description in a --help style format.
  /// @return string with the description
  virtual std::string helpDescription() const;

  /// Returns the default value.
  virtual Type defaultValue() const;

  /// Returns the default value, as a std::string, suitable for a configuration  file.
  virtual std::string getDefaultValueAsString() const;

  /// Returns the value, as a std::string
  virtual std::string getValueAsString() const;

  /// Called when the option is successfully processed. This is delegated
  /// to the arg callback object, if any; otherwise it propagates to the
  /// method in the base class.
  virtual void onProcessed();

  /// Called after all options have been successfully processed. This is
  /// delegated to the arg callback object, if any; otherwise it propagates
  /// to the method in the base class.
  virtual void onComplete(const OptionList& opts);

  /// Returns the low-level details about this object
  /// @return string with details
  virtual std::string debugInfo() const;

  /// Returns options in XML format
  /// @return string with description
  virtual std::string getXML() const;

  /// Sets values from command line parameters
  virtual CFuint setCommandLineValues(const std::string& arg,
                                      CFuint position,
                                      int argc,
                                      char** argv);

  /// Creates a new Option* with the same values as this
  Option* clone() const;

private: // data

  /// The default value.
  Type m_defvalue;

  /// The configured value.
  Type* m_value;

}; // class OptionT

//////////////////////////////////////////////////////////////////////////////

template < class Type >
OptionT<Type>::OptionT(const std::string& name,
                                   const std::string& description) :
Option(name, description),
m_defvalue(),
m_value(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
OptionT<Type>::OptionT(const OptionT& rhs) :
Option(rhs)
{
  operator=(rhs);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
OptionT<Type>& OptionT<Type>::operator=(const OptionT& rhs)
{
  Option::operator=(rhs);
  m_defvalue = rhs.m_defvalue;
  m_value = rhs.m_value;

  return *this;
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
OptionT<Type>::~OptionT()
{
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
Option* OptionT<Type>::clone() const
{
  return new OptionT<Type>(*this);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
void OptionT<Type>::linkToValue(void* value, const std::string& typestr)
{
  cf_assert( value != CFNULL );

  std::string thistype = DEMANGLED_TYPEID(Type);
  if( thistype != typestr )
    throw BadMatchException
     (FromHere(), "A type [" + typestr + "] was passed to option [" +
       getLabel()  + "] which requires a value of type [" + thistype + "]" );

  Type * ptr_value  = static_cast<Type*>(value);
  m_value = ptr_value;
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
void OptionT<Type>::setDefaultValue(void* defvalue, const std::string& typestr)
{
  cf_assert( defvalue != CFNULL );

  std::string thistype = DEMANGLED_TYPEID(Type);
  if( thistype != typestr )
    throw BadMatchException
     (FromHere(), "A type [" + typestr + "] was passed as default to option [" +
       getLabel()  + "] which requires a value of type [" + thistype + "]" );

  Type * ptr_value  = static_cast<Type*>(defvalue);
  m_defvalue = *ptr_value;
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
void OptionT<Type>::setValue(const std::string& value)
{
  StoragePolicy::setValue(value,this);
  this->validate();
  setLinkedOptions(value);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
Type OptionT<Type>::value() const
{
  return *m_value;
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
std::string OptionT<Type>::debugInfo() const
{
  return StoragePolicy::debugInfo(this);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
std::string OptionT<Type>::getXML() const
{
  return StoragePolicy::getXML(this);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
Type OptionT<Type>::defaultValue() const
{
  return m_defvalue;
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
std::string OptionT<Type>::getDefaultValueAsString() const
{
  return StoragePolicy::getDefaultValueAsString(this);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
std::string OptionT<Type>::getValueAsString() const
{
  return StoragePolicy::getValueAsString(this);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
void OptionT<Type>::onProcessed()
{
  Option::onProcessed();
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
void OptionT<Type>::onComplete(const OptionList& opts)
{
  Option::onComplete(opts);
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
inline std::string OptionT<Type>::helpDescription() const
{
  std::ostringstream os;
  os << StoragePolicy::description(this);
  return os.str();
}

//////////////////////////////////////////////////////////////////////////////

/**
 * "true" => true and "yes" => true.
 */
template <>
inline void OptionT<bool>::setValue(const std::string& value)
{
  std::string aslower = value;
  Common::StringOps::toLower(aslower);
  (*m_value) = (( aslower == "true" ) || ( aslower == "yes" ));
  this->validate();
  setLinkedOptions(value);
}

//////////////////////////////////////////////////////////////////////////////

/// prints "true" and "false", not "1" and "0".
template <>
inline std::string OptionT<bool>::helpDescription() const
{
  std::ostringstream os;
  os << description() << " (default = " << (m_defvalue ? "true" : "false") << ")";
  return os.str();
}

//////////////////////////////////////////////////////////////////////////////

/**
 * special case, because:
 * instead of:          option is:
 *     --foo true   ->  --foo
 *     --foo false  ->  --nofoo
 */
template <>
inline std::string OptionT<bool>::getDefaultValueAsString() const
{
  return m_defvalue ? "true" : "false";
}

//////////////////////////////////////////////////////////////////////////////

/**
 * special case, because:
 * instead of:          option is:
 *     --foo true   ->  --foo
 *     --foo false  ->  --nofoo
 */
template <>
inline std::string OptionT<bool>::getValueAsString() const
{
  return *m_value ? "true" : "false";
}

//////////////////////////////////////////////////////////////////////////////

template < class Type >
CFuint OptionT<Type>::setCommandLineValues(
  const std::string& arg,
  CFuint position,
  int argc,
  char** argv)
{
  return StoragePolicy::setCommandLineValues(arg,position,argc,argv,this);
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_OptionT_hh
