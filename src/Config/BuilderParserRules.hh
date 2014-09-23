// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Config_BuilderParserRules_hh
#define COOLFluiD_Config_BuilderParserRules_hh

//////////////////////////////////////////////////////////////////////////////

#include <string>
#include <map>

#include "Config/BuilderParserException.hh"
#include "Config/BuilderParserMandType.hh"
#include "Config/Config.hh"

typedef std::map<std::string, COOLFluiD::Config::BuilderParserMandType> BuilderParserFrameAttrs;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
 
 namespace Config {
  
//////////////////////////////////////////////////////////////////////////////
  
  /// @brief Manages protocol rules.
  
  /// @section overview Overview
  ///
  /// The purpose of this class is to provide a generic, protocol-independant,
  /// way to handle XML protocol rules. A frame parser can use it to check 
  /// whether a frame respects the protocol defined. A frame buider can use it 
  /// to build a frame that respects the difined protocol. Note that this class 
  /// is not capable of parsing or building a frame, but provides methods to 
  /// easily check the frame correctness. @n
  ///
  /// @section protocoldescr Protocol description
  /// 
  /// For this class, a protocol is typically a set of XML frames that can be 
  /// exchanged between two communicating entities (processes, objects, 
  /// threads,...). Each frame has a type (an unsigned integer to which a 
  /// string has been assigned), may content attributes and/or data is XML 
  /// format and is wrapped in a global tag, the root tag, defined by the 
  /// protocol and that allows to check that the frame belongs to this 
  /// protocol. Typically, the root tag is the protocol name. @n
  ///
  /// A frame as described above would be formatted like this : 
  /// \code
  /// <ProtocoleName>
  ///  <FrameType attr1="value1" attr2="value2">
  ///   <!-- (data) -->
  ///  </FrameType>
  /// </ProtocoleName>
  /// \endcode
  /// 
  /// A valid protocol:
  /// @li defines a <b>root type</b>, used to identify the protocol;
  /// @li defines an <b>error type</b>, used to return an invalid type;
  /// @li defines an error tag that is different from the root tag;
  /// @li defines a non-empty name for its root type;
  /// @li may not define a name for its error tag;
  ///
  /// @warning A frame type is considered as valid if and only if it as has 
  /// been bound to a non-empty name string. @n
  /// 
  /// @section mandatoriness Mandatoriness
  ///
  /// A protocol can define the mandatoriress of a frame attribute or frame 
  /// data. Available policies for this feature are defined by 
  /// @c #BuilderParserMandType. 
  /// \note @li @c #MAND_UNDEFINED will not be accepted if used to set a 
  /// mandatoriness policy since this enumeration item is used to signal an 
  /// error.
  /// @li only @c #MAND_OPTIONAL and @c #MAND_MANDATORY policies will be 
  /// accepted for attributes. Creating an attribute with @c #MAND_FORBIDDEN is
  /// the same as if it does not exist.
  /// @li an attribute with @c #MAND_MANDATORY policy must be present with a 
  /// non-empty value
  /// @li attributes having @c #MAND_OPTIONAL policy with an empty value are 
  /// ignored by both builder and parser.
  ///
  /// @section copies Instance copies
  ///
  /// Protocol rules are internally handled by maps. Because a protocol can
  /// define a big amount of rules, the memory used by an object of this class
  /// can easily become big. That is, it is recommended to avoid as much as
  /// possible copying an instance of this class.@n
  /// 
  /// @section newprotocol Create a new protocol
  /// 
  /// To create a new protocol, it is recommmended to subclass this class. 
  /// Defining a new rule, be done like this (assuming that @c MY_TYPE is a 
  /// enumeration item):
  /// @code
  /// this->setTypeName(MY_TYPE, "myType");
  /// this->addAttribute(MY_TYPE, "myAttribute", MAND_MANDATORY);
  /// this->setDataMandatoriness(MY_TYPE, MAND_OPTIONAL);
  /// @endcode
  /// It is recommended to manage frame types through a enumeration.
  /// @warning Please remember that the type name has to be set before 
  /// attributes or data mandatoriness. Not doing so will result on an error.
  /// 
  /// @author Quentin Gasper
    
  class Config_API BuilderParserRules
  {
   public:
   
    /// @brief Constructor.
    
    /// @param errorType Error type ID
    /// @param rootType Root type ID
    /// @param rootName Name for the root tag.
    /// @throw #BuilderParserException If @c errorType and @c rootType are 
    /// equal or if @c rootName is empty
    BuilderParserRules(unsigned int errorType, unsigned int rootType, 
                       const std::string & rootName);
   
    /// @brief Destructor.
    ~BuilderParserRules();
   
    /// @brief Check whether a type is valid. 
    
    /// A type is valid if it has a name associated. Thus @c #errorType is
    /// not considered as a valid type by this method.
    /// @param type The type to check.
    /// @return Returns @c true if the type is valid, otherwise returns @c false.    
    bool isValid(unsigned int type) const;

    /// @brief Gives the root type.
    
    /// @return Returns the root type of this protocol.
    unsigned int getRootType() const;
    
    /// @brief Gives the error type.
    
    /// @return Returns the error type of this protocol.
    unsigned int getErrorType() const;
    
    /// @brief Gives the name associated to the root type
    
    /// This is an overloaded method, provided for convinience.
    /// This is equivalent to calling @c #getTypeName with the root type as 
    /// parameter.
    /// @return Returns the name associated to the root type.
    std::string getRootName() const;
    
    /// @brief Gives the name associated to a given frame type.
    
    /// @param type Frame type. 
    /// @return Returns @c true if @c type was valid and type name was written
    /// into @c name; otherwise, returns @c false.
    std::string getTypeName(unsigned int type) const;
    
    /// @brief Gives the frame type to which the given name is associated.
    
    /// @param name Type name
    /// @return Returns the frame type associated, or the error type define by 
    /// the protocol if the name is not associated to any frame type.
    unsigned int getFrameType(const std::string & name) const;
    
    /// @brief Gives attributes associated to a given frame type.
    
    /// @param attrs Object to which attributes will be written. May be empty
    /// after a successful call if the frame type has no attribute.
    /// @param type Frame type.
    /// @param ok Pointer to a bool variable. If not @c CFNULL, the pointed 
    /// variable will be set to @c true on success, or to @c false on error. 
    /// @return Returns attributes associated to the node, or an empty
    /// object if @c type is not valid or the type has no attributes.
    BuilderParserFrameAttrs getAttributes(unsigned int type, 
                                               bool * ok = CFNULL) const;
    
    /// @brief Associates a name to a frame type. 
    
    /// If the frame already has name, it is replaced. Data mandatoriness is
    /// set to @c #MAND_FORBIDDEN if no name was associated before.
    /// @param type Frame type.
    /// @param name Frame type name
    /// @return Returns @c true if the name was successfuly associated, or 
    /// @c false if the name was empty or if the frame type is ther error type.
    bool setTypeName(unsigned int type, const std::string & name);
    
    /// @brief Adds an attribute to a frame type
    
    /// If an attribute with that name already exists for that frame type, the
    /// name is overwritten but data mandatoriness and attribute remain the same.
    /// @param type Frame type. 
    /// @param name Attribute name.
    /// @param mandType Mandatoriness policy. Default value is @c #MAND_MANDATORY.
    /// @return Returns @c false if @c mandType is @c #MAND_UNDEFINED or
    /// @c #MAND_FORBIDDEN, or if the frame type is not valid. Otherwise, 
    /// returns @c true.
    bool addAttribute(unsigned int type, const std::string & name, 
                      BuilderParserMandType mandType = MAND_MANDATORY);
    
    /// @brief Gives the mandatoriness policy of a given attribute in a given
    /// frame type.
    
    /// @param type Frame type. 
    /// @param name Attribute name
    /// @return Returns the mandatoriness policy, or @c #MAND_UNDEFINED if 
    /// @c type was not valid or if @c name attribute was not found for this 
    /// type.
    BuilderParserMandType getAttrMandatoriness(unsigned int type, 
                                               const std::string & name) const;
    
    /// @brief Sets data madatoriness policy for a frame type.
    
    /// @param type Frame type. 
    /// @param mandType New mandatoriness policy.
    /// @return Returns @c false if @c mandType is equal to @c #MAND_UNDEFINED
    /// or if the frame type is not valid. Otherwise, returns @c true.
    bool setDataMandatoriness(unsigned int type, BuilderParserMandType mandType);
    
    /// @brief Gives data mandatoriness policy for a given frame type.
    
    /// @param type Frame type. 
    /// @return Returns the mandatoriness policy, or @c #MAND_UNDEFINED if 
    /// @c type was not valid.
    BuilderParserMandType getDataMandatoriness(unsigned int type) const;
    
   private: // data
   
    /// @brief Error type ID.
    unsigned int m_errorType;
   
    /// @brief Root type ID.
    unsigned int m_rootType;
   
    /// @brief Frame names.
    
    /// Key is the type ID. Value is the associated name.
    std::map<unsigned int, std::string> m_frameNames;
    
    /// @brief Frame attributes.
    
    /// Key is the type ID. Value is the associated attributes and their 
    /// mandatoriness policy.
    std::map<unsigned int, BuilderParserFrameAttrs> m_frameAttributes;
    
    /// @brief Frame data mandatoriness policies.
    
    /// Key is the type ID. Value is the associated mandatoriness policy.
    std::map<unsigned int, BuilderParserMandType> m_frameDataMand;
    
    /// @brief Check parameters validity.
    
    /// @return Returns @c true if @c type is a valid type and @c str is
    /// a non-empty string.
    bool check(unsigned int type, const std::string & str) const;

  }; // class BuilderParserRules
 
//////////////////////////////////////////////////////////////////////////////

 } // namespace Config
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Config_BuilderParserRules_h