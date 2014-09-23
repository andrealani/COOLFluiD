#ifndef COOLFluiD_client_TypeAndNameDialog_h
#define COOLFluiD_client_TypeAndNameDialog_h

#include <QDialog>
#include <QObject>

class QComboBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QMainWindow;
class QDialogButtonBox;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    /// @brief Dialog used to add a node.
    
    /// This class inherits from @c QDialog and is used to show a 
    /// dialog allowing the user to create a new node. The dialog is modal, 
    /// which means that once it is visible, the calling code execution is 
    /// stopped until the dialog is invisible again. The user is invited to 
    /// type the name of the new node and select the concrete type of this 
    /// node. If the dialog has a parent window, it is centered on this parent.
    /// Otherwise, it is centered on the screen.@n @n
    
    /// After calling the constructor, the dialog is not visible. 
    /// @c show() method has to be called to show it. This is a 
    /// blocking method: it will not return until is invisible again. This 
    /// method returns either the name entered by the user (if he clicked on 
    /// "OK" to validate his entry) or an empty string (if he clicked on 
    /// "Cancel" or closed the dialog to cancel his entry).@n @n
    
    /// If the user validates his entry, the @c concreteType 
    /// parameter is used to store the selected concrete type. The method 
    /// guarantees that the selected type will be one of the provided list. If 
    /// the user cancels his entry, the parameter is not modified. If the user 
    /// clicks on "OK" without typing any name, it is considered as a 
    /// cancellation. @n @n
    
    /// A typical use of this class is (assuming that @c this is a
    /// @c QMainWindow object and @c concreteTypes is a
    /// @c QStringList with some concrete types) : @n
    /// @code
    /// TypeAndNameDialog dialog(this);
    /// QString type;                // used to store the chosen concrete type
    /// QString name = dialog.show(concreteTypes, type);
    ///
    /// if(name != "")
    /// {
    ///  // some treatements
    /// }
    /// @endcode
    
    /// @author Quentin Gasper.
    
    class TypeAndNameDialog : QDialog
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor.
      
      /// @param fieldLabelText Text to set on the buddy label of the text field.
      /// @param dropListText Text to set on the buddy label of the drop down list.
      /// @param parent Dialog parent. May be null.
      TypeAndNameDialog(const QString & fieldLabelText,
                        const QString & dropListText, 
                        QWidget * parent = NULL);
      
      /// @brief Destructor.
      
      /// Frees all allocated before the object is deleted. The parent is not
      /// destroyed.
      ~TypeAndNameDialog();
      
      /// @brief Shows the dialog.
      
      /// This is a blocking method. It will not return until the dialog is
      /// invisible. 
      /// @param types List of the available concrete types.
      /// @param concreteType Reference to a @c QString where the 
      /// selected type will be stored if and only if the user clicked on "OK" 
      /// and the name is not empty, otherwise the value is unchanged.
      /// @return If the user clicked on "OK", returns the name entered in the
      /// line edit component (may be empty if nothing was entered). Otherwise, 
      /// returns an empty string by calling the default @c QString 
      /// constructor. If the provided list is empty, an empty string is 
      /// returned.
      QString show(const QStringList & types, QString & concreteType);
      
      /// @brief Sets a new name value
      
      /// @param newName New name. May be empty.
      void setName(const QString & newName);
      
      private slots:
      
      /// @brief Slot called when "OK" button is clicked.    
      void btOkClicked();
      
      /// @brief Slot called when "Cancel" button is clicked.    
      void btCancelClicked();
      
    private:
      
      /// @brief Drop-down list that allows the user to select a concrete type.
      QComboBox * m_cbTypes;
      
      /// @brief Line edit that allows the user to enter the new object name.    
      QLineEdit * m_editName;
      
      /// @brief Button box containing "OK" and "Cancel" m_buttons.
      QDialogButtonBox * m_buttons;
      
      /// @brief The parent window.
      
      /// Can be null.
      QMainWindow * m_parent;
      
      /// @brief Label for the line edit
      QLabel * m_labName;
      
      /// @brief Label for the drop-down list
      QLabel * m_labConcreteType;
      
      /// @brief Layout on which the components will be placed.
      QFormLayout * m_layout;
      
      /// @brief Indicates whether the user clicked on "OK" button or not.
      
      /// If the user clicked on "OK" button, the attribute value is 
      /// @c true, otherwise (if the user closed the window or 
      /// clicked on "Cancel" button) it is @c false.
      bool m_okClicked;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_TypeAndNameDialog_h
