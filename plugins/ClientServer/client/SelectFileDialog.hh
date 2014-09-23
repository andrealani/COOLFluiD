#ifndef COOLFluiD_SelectFileDialog_h
#define COOLFluiD_SelectFileDialog_h

#include <QFileDialog>
#include <QHash>

namespace COOLFluiD
{
  namespace client
  {
    /// @brief Open/Save file dialog.
    
    /// This class subclasses @c QFileDialog. It allows to select a file to open 
    /// or to save. The main feature it adds to the base class is the file
    /// extension management. It will not allow user to select a file that does 
    /// not have the expected extension. @n @n
    
    /// The dialog is modal, wich means that once it is visible, the calling code 
    /// execution is stopped until the dialog is invisible again. If the dialog 
    /// has a parent window, it is centered on this parent. Otherwise, it is 
    /// centered on the screen.@n @n
    
    /// The class possesses two modes (inherited from the base class) : 
    /// @c QFileDialog::AcceptOpen (for a "Open file" dialog) and 
    /// @c QFileDialog::AcceptSave (for a "Save file" dialog). In @e open mode,
    /// user can not select a file that does not exist. In @e save mode, user
    /// can select a file that already exist, but he has to confirm that he is
    /// aware of overwriting this file. @n @n
    
    /// The class possesses its own extensions management, allowing the calling 
    /// code to specify the available file types and extensions, as the base class
    /// does, but it also garantees that only a file name with the extension 
    /// specified by the selected type will be returned. @n @n
    
    /// After calling the constructor, the dialog is invisible. The @c show
    /// method has to be called to show it. This is a blocking method: it will 
    /// not return until the dialog is invisible again. This method returns 
    /// either the name (with path) of the selected file (if user clicked on 
    /// "OK" to validate his selection) or an null string (if user clicked on 
    /// "Cancel" or closed the dialog to cancel his selection), built with default
    /// constructeur.@n @n
    
    /// The extension management system consists of two parts : the file type 
    /// (i.e. : "application", "text", "source",...) and the file extension. 
    /// The resulting string will be: @code <type> files (*.<extension>) @endcode 
    /// One type may have more than one extension, the resulting string will be 
    /// then: @code 
    /// <type> files (*.<extension_1> *.<extension_2> ... *.<extension_n>) 
    /// @endcode
    /// The file types and extensions must be set before calling @c show. New ones
    /// can be set between two consecutive calls to @c show on the same object.
    
    /// A typical use of this class to select a file where to save data is 
    /// (assuming that @c this is a @c QWidget object): @n @n
    /// @code
    /// SelectFileDialog dialog(this);
    /// QString filename;
    /// 
    /// dialog.addFileType("CFcase", "CFcase");
    /// dialog.addFileType("XML", "xml");
    /// dialog.addFileType("Web", "html htm css js");
    /// 
    /// filename = dialog.show(QFileDialog::AcceptSave);
    /// 
    /// if(!name.isNull())
    /// {  
    ///  // open the file
    /// }
    /// @endcode
    /// In the above code, if user selects "CFcase" @c show method will always 
    /// return a file name with "CFcase" extension. This is to garantee that 
    /// file extension and file type are always compatible.
    
    /// @warning <ul> 
    /// <li> Do not call base class method @c exec() to show the dialog or 
    /// the file extension management will not work properly. 
    /// <li> Never use the following base class methods: <ul>
    /// <li> <code>setAcceptMode(AcceptMode mode)</code>
    /// <li> <code>setConfirmOverwrite(bool enabled)</code>
    /// <li> <code>setNameFilter(const QString & filter)</code>@n @n
    /// The modifications they do on the object will be overwritten when 
    /// calling @c show method.
    /// </ul></ul>@n @n
    
    /// @author Quentin Gasper.
    class SelectFileDialog : public QFileDialog
    {
      Q_OBJECT
      
    public:
      /// @brief Constructor.
      
      /// @param parent Parent widget.
      SelectFileDialog(QWidget * parent = NULL);
      
      /// @brief Destructor.
      
      /// Frees all allocated memory. Parent is not destroyed.
      ~SelectFileDialog();
      
      /// @brief Shows the dialog.
      
      /// @param mode Dialog mode. Either @c QFileDialog::AcceptOpen (for a 
      /// "Open file" dialog) and @c QFileDialog::AcceptSave (for a "Save file" 
      /// dialog)
      /// @return Returns the file name selected, or a null string (built with
      /// @c QString ) if user clicked on "Cancel" or closed the dialog
      QString show(AcceptMode mode);
      
      /// @brief Adds the type and its associated file extensions.
      
      /// Some rules to know:
      /// @li if the type already exists, its value is replaced by the new one.
      /// @li all '*' and '?' characters are removed from @c extensions string.
      /// @li if @c extensions string contains an extension twice, it is set only 
      /// once
      /// @li if either trimmed @c type or trimmed @c extentions are empty, nothing is 
      /// done
      /// @li if both trimmed @c type and trimmed @c extentions are empty, 
      /// @c type is set to "All" and @c extensions remains empty. 
      /// An "All files" is then created where all files with or without 
      /// extension are accepted. This allows to have at most one "All files" type.
      /// 
      /// Examples:
      /// @code
      /// dialog.addFileType("Web", "ht*ml ht?m css js css");
      /// // the resulting string is : "Web files (*.html *.htm *.css *.js)"
      /// 
      /// dialog.addFileType("", "");
      /// // the resulting string is : "All files (*)"
      /// 
      /// dialog.addFileType(" ", "");     // does nothing
      /// dialog.addFileType("Web", "\n"); // does nothing
      /// dialog.addFileType("", "css");   // does nothing
      /// dialog.addFileType("Web", "");   // does nothing
      /// @endcode
      /// @param type File type.
      /// @param extensions Extensions. May contains several extension, separeted
      /// by a white space.
      void addFileType(const QString & type, const QString & extensions);
      
      /// @brief Removes a file type.
      
      /// If the file type does not exist, nothing is done.
      /// @param type File type name to remove.
      void removeFileType(const QString & type);
      
      /// @brief Appends one or more file extensions to a type.
      
      /// Some rules to know:
      /// @li if the type does not exist, it is created. 
      /// @li if a extension is already defined for this type, it is not 
      /// redefined. 
      /// @li all '*' and '?' characters are removed from @c extensions string. 
      /// @li if either trimmed @c type or trimmed @c extentions are empty, 
      /// nothing is done. 
      ///
      /// For example :
      /// @code
      /// dialog.addFileType("Web", "html htm");
      /// dialog.appendExtensionToFileType("Web", "css html js css");
      /// // the resulting string is : "Web files (*.html *.htm *.css *.js)"
      /// @endcode
      /// @param type Type name to which extension will be appended to.
      /// @param extensions Extensions. May contains several extension, separeted
      /// by a white space.
      void appendExtensionToFileType(const QString & type, const QString & extensions);
      
      /// @brief Removes all previously defined file types.
      void clearFileTypes();
      
      public slots:
      /// @brief Overrides @c QFileDialog::accept() slot for extension management
      /// system.
      virtual void accept();
      
    private:
      /// @brief Hash map for file types and extensions.
      
      /// The key is the file type and the value contains associated extension(s).
      QHash<QString, QString> m_extensions;
      
      /// @brief Current selected file name and path.
      QString m_filename;
      
      /// @brief Checks if @c filename has an extension associated to the select
      /// file type when selecting a file for save purpose.
      
      /// Three cases are possible: 
      /// @li the filename has a correct extension
      /// @li the filename has an extension that is not defined by any type or 
      /// has no extension. The current extension (if any) is removed and the
      /// first extension defined by the selected type is appended to the filename.
      /// @li the filename has an extension that is defined by another file type 
      /// than the selected one. User is resquested to choose between keep the 
      /// current file type, switch to the one that defines the extension or 
      /// cancel. <ul>
      /// <li> If user chooses to keep the file type, the extension is replaced with
      /// the first one defined by this type.
      /// <li> If user chooses to switch to the other file type, the selected file
      /// changed. The extension remains the same.
      /// <li> If user chooses to cancel, nothing is done.</ul>
      /// @param extensions Regular used to check if extension is defined by any
      /// file type.
      /// @return Returns @c false if user chooses to cancel in the third case. 
      /// Otherwise, always returns @c true.
      bool fixExtensionSave(const QRegExp & extensions);
      
      /// @brief Checks if @c filename has an extension associated to the select
      /// file type when selecting a file for open purpose.
      
      /// If the filename has an extension that is not defined by any type or 
      /// has no extension, the current extension (if any) is removed and the
      /// first extension defined by the selected type is appended to the 
      /// filename.
      void fixExtensionOpen();
      
      
      /// @brief Removes from a string every word that appears more than once to
      /// only one occurence of it.
      
      /// @param string String
      /// @return Returns the resulting string
      QString unique(const QString & string);
    };
  }
}

#endif // COOLFluiD_SelectFileDialog_h
