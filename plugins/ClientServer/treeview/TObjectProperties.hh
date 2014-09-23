#ifndef COOLFluiD_treeview_TObjectProperties_h
#define COOLFluiD_treeview_TObjectProperties_h

#include "ClientServer/treeview/TreeViewAPI.hh"

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace treeview
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Handles object properties.
    
    struct TreeView_API TObjectProperties
    {
      public :
      
      /// @brief Object type name.
      QString m_type;
      
      /// @brief Object abstract type name.
      QString absType;
      
      /// @brief If @c true, the object is basic, otherwise it is advanced.
      bool basic;
      
      /// @brief If @c true, the object is static, otherwise it is dynamic.
      bool dynamic;
      
      /// @brief Constructor.
      
      /// Provided for convinience.
      /// @param type Object type name.
      /// @param absType Object abstract type name.
      /// @param basic If @c true, the object is basic, otherwise it is advanced.
      /// @param dynamic If @c true, the object is dynamic, otherwise it is not.
      TObjectProperties(const QString & m_type = QString(), 
                        const QString & absType = QString(), 
                        bool basic = false, 
                        bool dynamic = false)
      {
        this->m_type = m_type;
        this->absType = absType;
        this->basic = basic;
        this->dynamic = dynamic;
      }
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace treeview
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_treeview_TObjectProperties_h