#ifndef COOLFluiD_client_RemoteFSBrowser_h
#define COOLFluiD_client_RemoteFSBrowser_h

/////////////////////////////////////////////////////////////////////////////

#include <QObject>
#include <QDialog>
#include <QDialogButtonBox>
#include <QIcon>
#include <QModelIndex>

#include "ClientServer/network/NetworkFrameType.hh"

class QCompleter;
class QDialogButtonBox;
class QDomDocument;
class QEvent;
class QHBoxLayout;
class QIcon;
class QLabel;
class QLineEdit;
class QListView;
class QMainWindow;
class QVBoxLayout;
class QSortFilterProxyModel;
class QStandardItemModel;

namespace COOLFluiD
{
  
  namespace treeview
  {
    class TreeModel;
  }
  
  namespace client
  {
    class ClientKernel;
    class FilesListItem;
    
    /// @brief This enum defines a validation m_protocol between @c RemoteFSBrowser
    /// and its subclasses. 
    
    /// It allows subclasses to inform the base class that the currently selected 
    /// item(s) is (are) conformed to the validation policy defined by these
    /// subclasses. See @c RemoteFSBrowser class documentation for more details.
    enum ValidationPolicy
    {
      /// @brief The associated selection has been validated by the subclass 
      /// policy.
      /// This validation can be used to tell the base class that everything is 
      /// "Ok", the dialog can be m_hidden.
      POLICY_VALID,
      
      /// @brief The associated selection has been invalidated by the subclass 
      /// policy.   
      /// This validation can be used to tell the base class that the selection 
      /// does not respect the validation policy, the dialog can not be m_hidden.
      POLICY_NOT_VALID,
      
      /// @brief The associated selection has been validated by the subclass 
      /// policy.
      /// This validation can be used to tell the base class that the selection is 
      /// a directory and respects the validation policy. Using this validation 
      /// means for the base class that the subclass accepts the selection but 
      /// not as a final selection. Thus, the base class will enter the selected 
      /// directory instead of hiding the dialog.
      POLICY_ENTER_DIRECTORY
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Dialog used to select a file to open.
    
    /// This class inherits from @c QDialog and is used to show a
    /// dialog allowing the user to select a file to open. The dialog is modal,
    /// wich means that once it is visible, the calling code execution is 
    /// stopped until the dialog is invisible again. If the dialog has a parent
    /// window, it is centered on this parent. Otherwise, it is centered on the
    /// screen.@n @n
    
    /// This class allows user to browse server files system. Double-clicking 
    /// on a directory will send a request to the server to open a directory 
    /// and return its contents. Thus the class needs a @c ClientNetworkComm 
    /// object with an open m_socket to communicate with the server. The 
    /// constructor sends a request to the server to open the default
    /// directory.@n @n
    
    /// After calling the constructor, the dialog is invisible. The @c show
    /// method has to be called to show it. This is a blocking method: it will 
    /// not return until the dialog is invisible again. This method returns 
    /// either  the name (with path) of the selected file (if he clicked on 
    /// "OK" to validate his selection) or an empty string (if he clicked on 
    /// "Cancel" or closed the dialog to cancel his selection).@n @n
    
    /// A typical use of this class is (assuming that @c this is a
    /// @c QMainWindow object): @n @n
    /// \code
    /// RemoteFSBrowser dialog(this);
    /// QString name = dialog.show();
    /// 
    /// if(name != "")
    /// {  
    /// // some treatements
    /// }
    /// \endcode
    
    /// @author Quentin Gasper.
    
    class RemoteFSBrowser : public QDialog
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor.
      
      /// @param index Simulation index.
      /// @param parent Parent window of the dialog. May be @c NULL.
      /// @throw std::invalid_argument If a connection to the server does not exist.
      RemoteFSBrowser(const QModelIndex & index, QMainWindow * parent = NULL);
      
      /// @brief Destructor.
      
      /// Frees all allocated memory. Parent window and @c communication object 
      /// are not destroyed.
      ~RemoteFSBrowser();
      
      /// @brief Shows the dialog. 
      
      /// This is a blocking method. It will not return while the dialog is
      /// visible. 
      /// @param startingDir Directory path where the browsing starts. If empty, 
      /// two cases are possible : 
      /// @li if the dialog has been already shown before, the current directory 
      /// is used as starting directory
      /// @li if the dialog is shown for the first time, the browsing starts in 
      /// the default directory defined by the server. 
      ///
      /// The path may be absolute or relative. If the string is null, it is 
      /// considered as empty.
      /// @return Returns the absolute path of the selected file to open or an
      /// empty string if user clicked on "Cancel" or closed the dialog.
      /// @see showMultipleSelect
      QString show(const QString & startingDir = "");
      
      /// @brief Shows the dialog. 
      
      /// This is a blocking method. It will not return while the dialog is
      /// visible. The difference with @c show() is that this method allows user 
      /// to select more than one file.
      /// @param startingDir Directory path where the browsing starts. If empty, 
      /// two cases are possible: 
      /// @li if the dialog has been already shown before, the current directory 
      /// is used as starting directory
      /// @li if the dialog is shown for the first time, the browsing starts in 
      /// the default directory defined by the server. 
      ///
      /// The path may be absolute or relative. If the string is null, it is 
      /// considered as empty.
      /// @return Returns the absolute paths of the selected files to open or an
      /// empty list if user clicked on "Cancel" or closed the dialog.
      /// @see show
      QStringList showMultipleSelect(const QString & startingDir = "");
      
      /// @brief Sets no extensions
      
      /// Old extensions are deleted.
      /// @param extensions New extensions
      void setExtensions(const QStringList & extensions);
      
      /// @brief Sets whether files need to be shown or not
      
      /// @param includeFiles If @c true, sub-directories and files will be 
      /// requested. If @c false, only sub-directories will be requested.
      /// @note This method may do nothing if it was disabled by the subclass 
      /// from which it is called. See @c #allowModifyBools for more information
      /// about this disabling feature. Read also the subclass documentation to 
      /// see if it disables this method.
      void setIncludeFiles(bool includeFiles);
      
      /// @brief Sets whether files without extension need to be shown or not
      
      /// @param includeNoExtension If @c true, files without any extension will 
      /// be requested. If @c false, they will not.
      /// @note This method may do nothing if it was disabled by the subclass 
      /// from which it is called. See @c #allowModifyBools for more information
      /// about this disabling feature. Read also the subclass documentation to 
      /// see if it disables this method.
      void setIncludeNoExtension(bool includeNoExtension);
      
      /// @brief Gives the current extensions
      
      /// @return Returns current extensions
      QStringList getExtensions() const;
      
      /// @brief Indicates whether files are included.
      
      /// @return Returns @c true if files are included; otherwise, returns 
      /// @c false.
      bool getIncludeFiles() const;
      
      /// @brief Indicates whether files without extension are included.
      
      /// @return Returns @c true if files without extension are included; 
      /// otherwise, returns @c false.
      bool getIncludeNoExtension() const;
      
      /// @brief Gives the current path
      
      /// @return Returns the current path, or an empty string if the dialog
      /// does not have any current path.
      QString getCurrentPath() const;
      
    protected:
      
      /// @brief Client kernel.
      ClientKernel * clientKernel;
      
      QModelIndex index;
      
      /// @brief Allows subclasses to disable some methods.
      
      /// Certain subclasses may need to prevent their calling code from executing 
      /// some methods susceptible to modify their purpose but may also need to 
      /// keep a public inheritance. For example, a subclass having as goal to 
      /// select a directory will call @c setIncludeFiles(false) so that only 
      /// directories will be shown. But it may also want to prohibit its calling 
      /// code to call @c setIncludeFiles(true) because doing so would not 
      /// respect the purpose of the subclass. This attribute allows to avoid 
      /// such calls. @n @n
      /// Concerned methods:
      /// @li @c #setIncludeFiles
      /// @li @c #setIncludeNoExtension
      ///
      /// If this attribute is set to @c false, calling one of these methods will 
      /// do nothing (and will not raise any error). @n
      /// The default value is @c true. An usual use of this attribute is to call
      /// these methods in the subclass constructor and then set the attribute to
      /// @c false. @n @n
      /// The use of this attribute is not mandatory and subclasses can 
      /// temporarily set back the value to @c true if needed. @n
      /// It is highly recommended to mention in the subclass documentation if
      /// these methods are disabled at a certain moment.
      bool allowModifyBools;
      
      /// @brief Allows subclasses to disable @c show method.
      
      /// If this attribute is set to @c false and @c show is called, an
      /// error message will be displayed, then @c show returns. Default value
      /// is @c true.
      bool allowSingleSelect;
      
      /// @brief Allows subclasses to disable @c showMultipleSelect method.
      
      /// If this attribute is set to @c false and @c showMultipleSelect is 
      /// called, an error message will be displayed, then @c showMultipleSelect 
      /// returns. Default value is @c true.
      bool allowMultipleSelect;
      
      /// @brief Method called when user presses a key or a combination of keys
      
      /// This methods overrides base class method and calls it before any other 
      /// treatement. Depending on what is pressed, the method has 3 working modes:
      /// @li If @e Enter key is pressed and one of the m_buttons is focus, it is 
      /// like if user clicked on this button;
      /// @li If either no modifier key (such as ctrl, shift, alt, etc...) is 
      /// pressed or shift key and another key are pressed, the combination is 
      /// appended to the filter and @c #editFilter is focused in. To 
      /// avoid confusion, list m_view is focused out.
      /// @li In all other cases, nothing is done. The method calls base class 
      /// method.
      /// @param event Event that occured.
      virtual void keyPressEvent(QKeyEvent * event);
      
      /// @brief Method called when user presses "Tab" key.
      
      /// This methods overrides base class method. If the path text edit has 
      /// the focus, the completer is displayed, the first item is selected and 
      /// appended to the path. If the path text edit does not have the focus,
      /// the base method is called.
      virtual bool focusNextPrevChild(bool next); 
      
      /// @brief Allows a subclass to add a button to the dialog button group.
      
      /// @note It is up to the subclass to destroy the created m_buttons in its
      /// destructor.
      /// @param text Text for the new button
      /// @param role Role of the new button. 
      /// @return Returns a pointer to the newly created button, or a @c NULL 
      /// pointer if the text is empty/null or the role is invalid.
      QPushButton * addButton(const QString & text, 
                              QDialogButtonBox::ButtonRole role);
      
      /// @brief Allows a subclass to check if a sub-directory or a file with
      /// a given name already exists in the current directory.
      
      /// @param name The name to search
      /// @return Returns @c true if a sub-directory or a file with the given 
      /// name exists in the current directory; otherwise, returns 
      /// @c false.
      bool itemExists(const QString & name) const;
      
      /// @brief Checks if a given name names a directory.
      
      /// The check is done locally, no request is done to the server. This means
      /// that the given item needs to in be the current directory.
      /// @param name Directory name to check.
      /// @return Returns @c true if a directory of that name was found. 
      /// Otherwise (including if a file of that name was found), returns 
      /// @c false.
      bool isDirectory(const QString & name) const;
      
      /// @brief Checks if a given name names a files.
      
      /// The check is done locally, no request is done to the server. This means
      /// that the given item needs to in be the current directory.
      /// @param name File name to check.
      /// @return Returns @c true if a file of that name was found. 
      /// Otherwise (including if a directory of that name was found), returns 
      /// @c false.
      bool isFile(const QString & name) const;
      
      
      /// @brief Checks if the selected m_items respect the policy.
      
      /// This method is called each time user clicks on "Ok" if @c show method
      /// was called and an item is selected. Subclasses that let @c show enabled 
      /// should override this method to define their own policy. For example, a 
      /// subclass used to select a file to open should reimplement this method 
      /// to return @c POLICY_VALID if the selected item is a file or 
      /// @c POLICY_ENTER_DIRECTORY if the selected item is a directory. 
      /// Default implementation always returns @c POLICY_VALID.
      /// @param name Path of the select directory or file.
      /// @param isDir @c true if @c name is a path to a directory, @c false if 
      /// it is a path to a file.
      /// @return Returns the validation.    
      virtual ValidationPolicy isAcceptable(const QString & name, bool isDir);
      
      /// @brief Checks if the selected item respect the policy.
      
      /// This method is called each time user clicks on "Ok" if 
      /// @c showMultipleSelect method was called and an item is selected. 
      /// Subclasses that let @c showMultipleSelect enabled should override this 
      /// method to define their own policy. For example, a subclass used to 
      /// select files should reimplement this method to return @c POLICY_VALID 
      /// if the selected m_items are all files or @c POLICY_NOT_VALID if the 
      /// selected m_items contain at least one directory. Default implementation 
      /// always returns @c POLICY_VALID.
      /// @param name Path of the select directory or file.
      /// @param isDir @c true if @c name is a path to a directory, @c false if 
      /// it is a path to a file.
      /// @return Returns the validation.
      virtual ValidationPolicy isAcceptable(const QStringList & names);
      
      /// @brief Allows subclasses to reinitialize their environment before the 
      /// dialog is showed again.
      
      /// This method is called each time the dialog is showed. Subclasses may 
      /// override this method to reset their internal data to default values, 
      /// if needed. The default implementation does nothing.
      virtual void reinitValues();
      
      /// @brief Displays a message box with an error message
      
      /// @param message Message to display
      void showError(const QString & message);
      
      /// @brief Cat two path strings.
      
      /// Path separator is added, is needed. The created path is guaranteed to be 
      /// compatible with operating system on which the server is running.
      /// @param part1 First part
      /// @param part2 Second part
      /// @warning @c part1 parameter is modified: @c part2 is directly appended 
      /// to it.
      void assemblePath(QString & part1, const QString & part2) const;
      
      /// @brief Gives the selected file or directory.
      
      /// This method is called at the end of @c show method. Subclasses that 
      /// reimplemented <code>isAcceptable(const QString &, bool)</code> should
      /// reimplement this method.
      /// @return Return the selected file directory, or an empty string if no
      /// item is selected.
      virtual QString getSelectedFile() const;
      
      /// @brief Gives the selected files and/or directories.
      
      /// This method is called at the end of @c showMultipleSelect method. 
      /// Subclasses that reimplemented <code>isAcceptable(const QStringList &)</code> 
      /// should reimplement this method.
      /// @return Return the selected file directory, or an empty string if no
      /// item is selected.
      virtual QStringList getSelectedFileList() const;
      
      /// @brief Allows a subclass to set text on the satus label.
      
      /// @param text Text for the status label
      void setStatus(const QString & text);
      
      private slots:
      
      /// @brief Slot called when "OK" button is clicked. 
      
      /// Sets @c #m_okClicked to @c true and then sets 
      /// the dialog to an invisible state.
      void btOkClicked();
      
      /// @brief Slot called when "Cancel" button is clicked. 
      
      /// Sets @c #m_okClicked to @c false and then sets the dialog to an 
      /// invisible state.
      void btCancelClicked();
      
      /// @brief Slot called each time the text in the line edit is modified. 
      
      /// The filter is modified according to this text and the list m_view is
      /// updated.
      /// @param text New text in the line edit.
      void filterUpdated(const QString & text);
      
      /// @brief Slot called when path line edit has been modified
      
      /// @param text New text in the line edit.
      void pathUpdated(const QString & text);
      
      /// @brief Slot called when the server sends the contents of a directory 
      /// to the client. 
      
      /// The model is updated and the line edit is cleared.
      /// @param path Absolute path of the directoy of which contents belong to.
      /// @param dirs Directories list. Each element is a directory.
      /// @param files Files list. Each element is a file.
      void dirContents(const QString & path, const QStringList & dirs, 
                       const QStringList & files);
      
      /// @brief Slot called when an error comes from the network level.
      
      /// @param error
      /// @param fromServer @c true if the error message comes from the server, 
      /// otherwise @c false. This parameter is never used.
      void error(const QString & error, bool fromServer);
      
      /// @brief Slot called when the user double-click on an item in the 
      /// list m_view. 
      
      /// If the item is a directory, a request to open this directory is sent 
      /// to the server. If the item is a file, the dialog reacts as if the 
      /// user clicked on "OK" button with this item selected.
      /// @param index Clicked item index in the model filter.
      void doubleClick(const QModelIndex & index);
      
      /// @brief Slot called if the completer has been activated
      
      /// @param text Text that activated the completer
      void completerActivated(const QString & text);
      
      /// @brief Slot called if an ack comes from the server
      
      /// Only "create directory" ACK is taken in account. Others are ignored.
      /// @param type Type of the ACK
      void ack(COOLFluiD::network::NetworkFrameType type);
      
    private:
      
      /// @brief Label for the filter.
      QLabel * m_labFilter;
      
      /// @brief Label for the files list.
      QLabel * m_labFilesList;
      
      /// @brief Text edit for the current directory path
      QLineEdit * m_editPath;
      
      /// @brief Layout for path widgets (label and text edit)
      QHBoxLayout * m_pathLayout;
      
      /// @brief Bottom m_layout (filter and m_buttons)
      QHBoxLayout * m_bottomLayout;
      
      /// @brief Line edit for the filter.
      QLineEdit * m_editFilter;
      
      /// @brief Button box for "OK" and "Cancel" m_buttons.
      QDialogButtonBox * m_buttons;
      
      /// @brief List m_view used to display files list.
      QListView * m_listView;
      
      /// @brief Model for the list m_view
      QStandardItemModel * m_viewModel;
      
      /// @brief Model for the completer
      QStandardItemModel * m_completerModel;
      
      /// @brief Automatic completion for path
      QCompleter * m_pathCompleter;
      
      /// @brief Filter for the model.
      QSortFilterProxyModel * m_filterModel;
      
      /// @brief Dialog m_layout.
      QVBoxLayout * m_layout;
      
      /// @brief Indicates whether the user clicked on "OK" button or not. 
      
      /// If the user clicked on "OK" button, the attribute value is @c true, 
      /// otherwise (if the user closed the window or clicked on "Cancel" button) 
      /// it is @c false.
      bool m_okClicked;
      
      /// @brief Path of the current directory.
      QString m_currentPath;
      
      /// @brief File selected by the user.
      QString m_currentFile;
      
      /// @brief Files selected by the user is @c showMultipleSelect() was called.
      QStringList m_currentFilesList;
      
      /// @brief Parent window
      QMainWindow * m_parentWindow;
      
      /// @brief Indicates whether user is able to select multiple item at the 
      /// same time.
      bool m_multipleSelectAllowed;
      
      /// @brief Indicates whether files should be asked to the server.
      
      /// If @c true, files will be asked; otherwise, they will not.
      bool m_includeFiles;
      
      /// @brief List of extensions.
      
      /// Only files having one of these extensions will be returned by the server.
      /// Note that this attribute is ignored if @c #includeFiles is @c false.
      QStringList m_extensions;
      
      /// @brief Indicates wether files without extension should be asked to 
      /// the server.
      
      /// Note that this attribute is ignored if @c #includeFiles is @c false.
      bool m_includeNoExtension;
      
      /// @brief Indicates wether the completer is being updated.
      
      /// This attribute is @c true if the completer is waiting for the server to
      /// send a directory content.
      bool m_updatingCompleter;
      
      /// @brief List of m_items in the list
      QList<FilesListItem *> m_items;
      
      /// @brief List of m_items in the completer
      QList<FilesListItem *> m_itemsCompleter;
      
      /// @brief Path separator
      
      /// This represents the path separator on the operating system on which 
      /// the server is running.
      QString m_pathSep;
      
      /// @brief Satus label
      QLabel * m_labStatus;
      
      /// @brief Applies new data to a standard item model.
      
      /// Old model data are deleted.
      /// @param model Model to update
      /// @param path The path. It will be prepended to all item. May be empty.
      /// @param dirs Directory list
      /// @param files File list
      /// @param modelItems List where new m_items will be stored. This list
      /// is cleared before usage.
      void updateModel(QStandardItemModel * model, const QString & path,
                       const QStringList & dirs, const QStringList & files, 
                       QList<FilesListItem *> & modelItems);
      
      /// @brief Sends a request to the server to open a directory
      
      /// This methods change the mouse cursor to a "wait cursor" (usualy a
      /// hourglass).
      /// @param path Directory path to open
      void openDir(const QString & path);
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_RemoteFSBrowser_h
