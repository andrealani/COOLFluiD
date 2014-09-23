#ifndef COOLFluiD_client_MenuActionInfo_h
#define COOLFluiD_client_MenuActionInfo_h

/////////////////////////////////////////////////////////////////////////////

class QIcon;
class QKeySequence;
class QMenu;
class QString;
class QToolBar;
class QAction;
class QString;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Manages an action configuration.
    
    /// The goal of this structure is to make action creation process more clear.
    /// Instead of giving many parameters (sometimes with a default value) to
    /// function that will build an action, @c MenuActionInfo allows to only set
    /// useful parameters. Other ones will have a default value.
    
    /// @author Quentin Gasper.
    
    struct MenuActionInfo : public QObject
    {
      public :
      
      /// @brief Action text.
      
      /// The default value is an empty string.
      QString m_text;
      
      /// @brief Slot to call when the action is triggered.
      
      /// It is recommended to use Qt's @c SLOT() macro. The default value is 
      /// @c NULL.
      const char * m_slot;
      
      /// @brief Indicates wether the action is enabled.
      
      /// The default value is @c true.
      bool m_enabled;
      
      /// @brief Indicates whether the action is a checkable info.
      
      /// The default value is @c false.
      bool m_checkable;
      
      /// @brief The menu the action will be added to.
      
      /// The default value is @c NULL.
      QMenu * m_menu;
      
      /// @brief The toolbar the action will be added to.
      
      /// The default value is @c NULL.
      QToolBar * m_toolbar;
      
      /// @brief Action icon.
      
      /// The default value is an empty string.
      QIcon m_icon;
      
      /// @brief Shortcut to trigger the action.
      
      /// The default value is an empty key sequence built using default 
      /// @c QKeySequence constructor..
      QKeySequence m_shortcut;
      
      /// @brief Constructor.
      
      /// Initializes all parameters to their default value.
      MenuActionInfo();
      
      /// @brief Initializes all parameters to their default value.
      
      /// @note Zones pointed by @c #slot, @c #menu and @c #toolbar are @b never
      /// destroyed.
      void initDefaults();
      
      /// @brief Builds an action from this config.
      
      /// @c NULL pointers are ignored. If a toolbar is set but no icon, icon
      /// is replace by action text. 
      /// @param parent Parent object that contains @c slot. This object will be
      /// thie action parent.
      /// @return Returns the built action or a @c NULL pointer if the parent
      /// is @c NULL.
      QAction * buildAction(QObject * parent);
      
    }; // struct MenuActionInfo
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_MenuActionInfo_h
