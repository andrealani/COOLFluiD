#ifndef COOLFLuiD_client_OptionTypes_h
#define COOLFLuiD_client_OptionTypes_h

/////////////////////////////////////////////////////////////////////////////

#include <QHash>

namespace COOLFluiD
{
  namespace client
  {
    
    /// @brief Defines available option types
    /// @todo Find a better name !
    /// @author Quentin Gasper.
    
    enum TOptionTypes
    {
      /// @brief Type id used to indicate an invalid type.
      NO_TYPE,
      
      /// @brief Type id used to indicate a @e boolean type.
      TYPE_BOOL,
      
      /// @brief Type id used to indicate a @e signed integer type.
      TYPE_INT,
      
      /// @brief Type id used to indicate a @e integer type.
      TYPE_UNSIGNED_INT,
      
      /// @brief Type id used to indicate a @e double type.
      TYPE_DOUBLE,
      
      /// @brief Type id used to indicate a @e string type.
      TYPE_STRING,
      
      /// @brief Type id used to indicate a @e files type.
      TYPE_FILES,
      
      /// @brief Type id used to indicate a @e libraries type.
      TYPE_LIBRARIES,
      
      /// @brief Type id used to indicate a <i>selectable item list</i> type.
      TYPE_HOST_LIST
    };
    
    /// @brief Defines and manages the option types.
    
    /// @author Quentin Gasper.
    
    /////////////////////////////////////////////////////////////////////////////
    
    class OptionTypes
    {
      
    public:
      
      /// @brief Checks if a type id is valid.
      
      /// A type id is valid if it exists and is it has a type name associated. 
      /// Thus @c TOptionTypes::NO_TYPE will not be considered as valid by this 
      /// function.
      /// @param id The type id to check.
      /// @return Returns @c true if the type id is valid, 
      /// otherwise returns @c false.
      static bool isValid(TOptionTypes id);
      
      /// @brief Gives the type id of a given type name.
      
      /// @param type The type name.
      /// @return Returns the type id corresponding to the given type name, or
      /// @c TOptionTypes::NO_TYPE if the type name is unknown.
      static TOptionTypes getTypeId(const QString & type);
      
      /// @brief Gives the type name for a given type id.
      
      /// @param type The type id.
      /// @return Returns the type name for the provided type id, or an empty 
      /// string if the type id does not exist or if it is @c TOptionTypes::NO_TYPE.
      static QString getTypeString(TOptionTypes type);
      
      /// @brief Gives a types list.
      
      /// This list contains all types that have a name associated to their id. 
      /// This list is not sorted.
      /// @returns Returns the types list.
      static QStringList getTypesList();
      
    private:
      
      /// @brief Hash map with all types.
      
      /// The key is the type id defined by one the public constant interger
      /// attributes of this class. The value is the type name for this id. All 
      /// types ids have a name except @c TOptionTypes::NO_TYPE.
      static QHash<TOptionTypes, QString> types;
      
      /// @brief Builds the types hash map.
      
      /// This function builds the hash map at most once during runtime. If it 
      /// is called a second time, it returns without doing anything. 
      static void buildTypes();
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFLuiD_client_OptionTypes_h
