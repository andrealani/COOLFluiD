#ifndef COOLFluiD_client_GraphicalOption_h
#define COOLFluiD_client_GraphicalOption_h

/////////////////////////////////////////////////////////////////////////////

#include "ClientServer/client/OptionTypes.hh"

class QFormLayout;
class QHBoxLayout;
class QLabel;
class QLineEdit;
class QWidget;

namespace COOLFluiD
{
  namespace client
  {
    
    /// @brief Displays an option graphically.
    
    /// The value component is adapted to the type of the option.
    
    /////////////////////////////////////////////////////////////////////////////
    
    class GraphicalOption
    {
      
    public:
      
      /// @brief Constructor. 
      
      /// @param type Option type. Must be one of those defined by TOptionTypes 
      /// enum.
      GraphicalOption(TOptionTypes type);
      
      /// @brief Destructor.
      
      /// Frees all allocated memory.
      ~GraphicalOption();
      
      /// @brief Gives the option name.
      
      /// @return Returns the option name. 
      QString getName() const;
      
      /// @brief Sets option name.
      
      /// @param name Option name.
      void setName(const QString & name);
      
      /// @brief Gives the option value.
      
      /// @return Returns the option value.
      /// @note If the type of the option is @c OptionTypes::TYPE_FILES,
      /// the variant object returned is a @c QStringList.
      QVariant getValue() const;
      
      /// @brief Gives the option value in string format.
      
      /// This is an overloaded method, provided for convinience. This method 
      /// calls @c getValue() and converts the returned variant object 
      /// to a @c QString.
      /// @return Returns the value in a string.
      /// @note If the type of the option is @c OptionTypes::TYPE_FILES,
      /// this method returnes a string with files separated by a white space.
      QString getValueString() const;
      
      TOptionTypes getType() const;
      
      /// @brief Adds this option to the provided m_layout.
      
      /// @param m_layout Layout to which the m_options has to be added.
      void addToLayout(QFormLayout * layout);
      
      /// @brief Sets a new value to the option.
      
      /// @param value New value. Must be in a format compatible with the option 
      /// type. Compatible formats list for each type is available in QVariant 
      /// class documentation.
      /// @throw InvalidValueException If the value could not be converted to
      /// the option type. 
      /// @note If the type of the option is @c OptionTypes::TYPE_FILES or 
      /// @c OptionTypes::TYPE_LIBRAIRIES, the value is a @c QStringList. 
      /// Neither files existence nor paths rightness are checked. The list 
      /// may be empty.
      void setValue(const QVariant & value);
      
      /// @brief Enables or disables the value component.
      
      /// If the component is enabled, its value is modifiable.
      /// @param enabled If @c true, the component is enabled. 
      /// Otherwise it is disabled.
      void setEnabled(bool enabled);
      
      /// @brief Indicates wether the value component is enabled or not.
      
      /// @return Returns @c true if the component is enabled.
      bool isEnabled() const;
      
      /// @brief Sets a tooltip.
      
      /// @param toolTip Tool tip to set.
      void setToolTip(const QString & toolTip);
      
      /// @brief Indicates wether the value has been modified.
      
      /// The value is considered to have been modified if it is different from
      /// the last value assigned using @c setValue.
      /// @return Returns @c true if the value has been modified.
      bool isModified() const;
      
      /// @brief Gives the option original value.
      
      /// The original value is the last one set using @c setValue.
      /// @return Returns the option original value.
      /// @note If the type of the option is @c OptionTypes::TYPE_FILES,
      /// the variant object returned is a @c QStringList.
      QVariant getOrginalValue() const;
      
      /// @brief Gives the option value in string format.
      
      /// This is an overloaded method, provided for convinience. This method 
      /// calls @c getOrginalValue() and converts the returned variant object 
      /// to a @c QString.
      /// @return Returns the value in a string.
      /// @note If the type of the option is @c OptionTypes::TYPE_FILES or 
      /// @c OptionTypes::TYPE_LIBRAIRIES, this method returnes a string with 
      /// files separated by a white space.
      QString getOrginalValueString() const;
      
      /// @brief Commits changes.
      
      /// Calling this method will set the current option value, given by
      /// @c #getValue, as original value.
      void commit();
      
    private:
      
      /// @brief Label for the option name.
      QLabel * m_name;
      
      /// @brief Line edit for the option value.
      QWidget * m_valueWidget;
      
      /// @brief Type of the option, according to the type ids defined by  
      /// OptionTypes class.
      TOptionTypes m_type;
      
      /// @brief Indicates wether the value component is enabled (allows 
      /// modification) or not.
      bool m_enabled;
      
      QVariant m_originalValue;
      
      
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_GraphicalOption_h
