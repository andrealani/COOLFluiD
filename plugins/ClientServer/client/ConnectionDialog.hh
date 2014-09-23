#ifndef COOLFluiD_client_ConnectionDialog_h
#define COOLFluiD_client_ConnectionDialog_h

/////////////////////////////////////////////////////////////////////////////

#include <QDialog>

#include "ClientServer/treeview/TSshInformation.hh"

class QCheckBox;
class QDialogButtonBox;
class QFormLayout;
class QHBoxLayout;
class QLabel;
class QLineEdit;
class QMainWindow;
class QSpinBox;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    //   struct TSshInformation;
    
    /// @brief Dialog used to gather information to connect to a server.
    
    /// This class inherits from @c QDialog and is used to show a 
    /// dialog gathering the needed information to connect to a server. <b>It 
    /// does not realize this connection</b>. The dialog is modal, wich means 
    /// that once it is visible, the calling code execution is stopped until 
    /// the dialog is invisible again. If the dialog has a parent window, it is
    /// centered on this parent. Otherwise, it is centered on the screen.@n @n
    
    /// This dialog has two modes : basic and advanced. The difference is that 
    /// in advanced mode, the user is able to choose the port number, but not 
    /// in basic mode. In both modes, the user is invited to enter the hostname 
    /// to connect to (default value: @e localhost ) and, if he chose to 
    /// start a new server m_instance, the username used to authenticate on the 
    /// remote machine (default value: the username of the process owner).@n @n
    
    /// After calling the constructor, the dialog is invisible. The show method 
    /// has to be called to show it. This is a blocking method: it will not 
    /// return until is invisible again. This method returns @c true 
    /// if user clicked on "OK" to validate his entry or @c false if 
    /// he clicked on "Cancel" or closed the dialog to cancel his entry.@n @n
    
    /// If the user validates his entry, gathered information are stored in 
    /// the TSshInformation structure parameter. If the user did not choose to 
    /// launch a new server m_instance, @c username parameter is this
    /// structure is not modified. The method guarantees that all other 
    /// attributes will be correctly set. If the user cancels his entry, the
    /// structure is not modified. If the user clicks on "OK" without typing 
    /// any name, it is considered as a cancellation.@n @n
    
    /// A typical use of this class is (assuming that @c this is a @c QMainWindow 
    /// object) : @n 
    /// @code
    /// ConnectionDialog dialog(this);  
    /// TSshInformation sshInfos;        // used to store gathered information
    /// 
    /// if(dialog.show(sshInfos)         // show advanced connection dialog
    /// {                                // if user clicked on "OK"
    ///  // some treatements
    /// }
    /// @endcode
    
    /// @author Quentin Gasper.
    
    class ConnectionDialog : public QDialog
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor.
      
      /// @param parent Parent window.
      ConnectionDialog(QMainWindow * parent);
      
      /// @brief Desctructor.
      ~ConnectionDialog();
      
      /// @brief Shows the dialog. 
      
      /// This is a blocking method. It will not return while the dialog 
      /// is visible.
      /// @param hidePort If @c true, user will not be able to select the port 
      /// number.
      /// @param sshInfos Reference to TSshInformation structure where grabbed
      /// information will be written if and only if the user clicked on "OK" 
      /// and the name is not empty, otherwise the structure is unchanged.
      /// @return If the user clicked on "OK", returns "true". Otherwise, 
      /// returns @c false.
      bool show(bool hidePort, COOLFluiD::treeview::TSshInformation & sshInfos);
      
      void setSshInfos(const COOLFluiD::treeview::TSshInformation & sshInfos);
      
      public slots:
      
      /// @brief Slot called when "OK" button is clicked.
      
      /// Sets @c #m_okClicked to @c true and then sets the dialog to an 
      /// invisible state.
      void btOkClicked();
      
      /// @brief Slot called when "Cancel" button is clicked. 
      
      /// Sets @c #m_okClicked to @c false and then 
      /// sets the dialog to an invisible state.
      void btCancelClicked();
      
      /// @brief Slot called when @c #m_chkLaunchServer has been checked or 
      /// unchecked. 
      
      /// If it is checked, username line edit will be enabled for 
      /// modification, otherwise it will be disabled.
      /// @param state New state of @c #m_chkLaunchServer (based on 
      /// @c Qt::CheckState enum). If the value is @c Qt::Checked, the username 
      /// line edit is set to enabled.
      void chkLaunchServerChecked(int state);
      
    private:
      
      /// @brief Label for the hostame line edit.
      QLabel * m_labHostname;
      
      /// @brief Label for the username line edit.
      QLabel * m_labUsername;
      
      /// @brief Label for the port number spin box.
      QLabel * m_labPortNumber;
      
      /// @brief Line edit for the hostame.
      QLineEdit * m_editHostname;
      
      /// @brief Line edit for the username.
      QLineEdit * m_editUsername;
      
      /// @brief Spin box for the port number.
      QSpinBox * m_spinPortNumber;
      
      /// @brief Main m_layout.
      QFormLayout * m_layout;
      
      /// @brief Button box containing "OK" and "Cancel" m_buttons.
      QDialogButtonBox * m_buttons;
      
      /// @brief Layout for hostname and port number components.
      QHBoxLayout * m_infosLayout;
      
      /// @brief Ckeck box to check of the user wants to launch a new 
      /// server m_instance.
      QCheckBox * m_chkLaunchServer;
      
      /// @brief Indicates whether the user clicked on "OK" button. 
      
      /// If the user clicked on "OK" button, the attribute value is @c true, (if 
      /// the user closed the window or clicked on "Cancel" button) it is @c false.
      bool m_okClicked;
      
      QSize m_minSize;
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_ConnectionDialog_h
