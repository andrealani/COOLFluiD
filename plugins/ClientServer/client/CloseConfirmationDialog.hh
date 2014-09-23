#ifndef COOLFluiD_client_CloseConfirmationDialog_h
#define COOLFluiD_client_CloseConfirmationDialog_h

/////////////////////////////////////////////////////////////////////////////

#include <QDialog>
#include <QHash>

class QDialogButtonBox;
class QGridLayout;
class QLabel;
class QMainWindow;
class QTextEdit;
class QVBoxLayout;

namespace COOLFluiD
{
  namespace client
  {
    class CloseConfirmationPanel;
    struct CloseConfirmationInfos;
    
    enum CloseConfirmationType
    {
      CLOSE_SAVE_FILE,
      CLOSE_COMMIT,
      CLOSE_SHUT_DOWN
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class is a dialog that asks important questions to the user
    /// before he closes the application.
    
    /// The aim of this class is to preserve user's work from loss of data. When
    /// user closes the application, the class is called if some elements request
    /// his attention. He is requested to solve all the "problems" before the
    /// application closes. The dialog guarantees that it will not be closeable
    /// before user answers to all questions.@n @n
    
    /// After calling the constructor, the dialog is invisible. The @c #show 
    /// method has to be called to show it. This is a blocking method: it will  
    /// not return until is invisible again. This method returns @c true 
    /// if user clicked on "OK" to validate his entry or @c false if 
    /// he clicked on "Cancel" or closed the dialog to cancel his entry.@n @n
    
    /// Calling @c #show method directly after the constructor will do nothing
    /// and directly return @c true. This is because no confirmation panel has
    /// been set. Each confirmation panel asks one question and checks user's 
    /// entry, independantly of other m_panels. To work, this dialog needs at least
    /// one panel and at most one panel of each available type. Avaiblable types
    /// are defined @c #CloseConfirmationType enum. Panels can be set by calling 
    /// @c #addConfirmation method before callig @c #show. 
    
    /// If the user validates his entry, gathered information are stored in 
    /// the @c #CloseConfirmationInfos structure parameter. If the user cancels 
    /// his entry, the structure is not modified. @n @n
    
    /// A typical use of this class is (assuming that @c this is a @c QMainWindow 
    /// object) : @n 
    /// @code
    /// CloseConfirmationDialog dialog(this);  
    /// CloseConfirmationInfos infos;    // used to store gathered information
    /// 
    /// dialog.addConfirmation(CLOSE_COMMIT);
    /// 
    /// if(dialog.show(infos)            // show advanced connection dialog
    /// {                                // if user clicked on "OK"
    ///  // some treatements
    /// }
    /// @endcode
    
    /// @author Quentin Gasper.
    
    class CloseConfirmationDialog : public QDialog
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor
      
      /// @param parent Parent window. May be @c NULL.
      CloseConfirmationDialog(QMainWindow * parent = 0);
      
      /// @brief Destructor
      
      /// Frees all allocated memory. Parent is not destroyed.
      ~CloseConfirmationDialog();
      
      /// @brief Adds a confirmation panel to the dialog
      
      /// If a panel of this type already exists in this dialog, nothing is done.
      /// @param type Panel type. Available value are defined in 
      /// @c #CloseConfirmationType enum.
      /// @param becauseCommit If @c type is @c #CLOSE_SAVE_FILE, indicates 
      /// whether the presence of this panel is due to a commit. See 
      /// @c SaveConfirmationPanel class documentation for more details. This 
      /// parameter if ignored for other types.
      void addConfirmation(CloseConfirmationType type, bool becauseCommit = false);
      
      /// @brief Shows the dialog.
      
      /// Before the dialog is showed, @link CloseConfirmationPanel::setData 
      /// @c setData @endlink methods of set m_panels are called to set default 
      /// data (if any).
      /// @param infos Structure where user's entries will be stored if he clicks
      /// on "OK". If he does not, structure is not modified. 
      /// @return Returns @c true if user clicked on "OK", @c false if he did not.
      /// If no comfirmation panel is set, the method returns directly @c true 
      /// (dialog is not shown).
      bool show(CloseConfirmationInfos & infos);
      
      private slots:
      
      /// @brief Slot called when user clicks on "OK"
      void btOkClicked();
      
      /// @brief Slot called when user clicks on "Cancel"
      void btCancelClicked();
      
      /// @brief Slot called when a panel requests to show some help.
      
      /// @param m_helpText Help text to show
      void showHelp(const QString & helpText);
      
      /// @brief Slot called when a panel has been resized.
      
      /// This slot resizes the dialog to fit the new panel size.
      void panelResized();
      
    private: // data
      
      /// @brief Label at top of the window displaying a text to explain the
      /// situation.
      QLabel * m_labInstructions;
      
      /// @brief "OK" and "Cancel" m_buttons.
      QDialogButtonBox * m_buttons;
      
      /// @brief Dialog m_layout.
      QVBoxLayout * m_layout;
      
      /// @brief Text area that displays help about a question.
      QTextEdit * m_editHelp;
      
      /// @brief Dialog where help is displayed.
      QDialog * m_helpDialog;
      
      /// @brief List of available m_panels.
      
      /// The key is the type of the panel, the value is either a pointer to 
      /// the object (if the panel has been set) or @c NULL (if the panel has not 
      /// been set).
      QHash<CloseConfirmationType, CloseConfirmationPanel *> m_panels;
      
      /// @brief Indicates whether user clicked on "OK".
      
      /// If @c true, user clicked on "OK".
      bool m_okClicked;
      
      /// @brief Indicates whether help has never been shown
      
      /// If @c true, help was never shown. In this case, help dialog will be 
      /// correctly placed on next show.
      bool m_helpNeverShown;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_CloseConfirmationDialog_h