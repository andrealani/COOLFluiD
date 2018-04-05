#ifndef COOLFluiD_RadiativeTransfer_UTILITIES_STRING_UTILS_H
#define COOLFluiD_RadiativeTransfer_UTILITIES_STRING_UTILS_H

#include <string>
#include <vector>

namespace String {

/**
 * Tokenizes a string into a vector of sub strings that are separated by the 
 * characters in teh delim string.
 */
void tokenize(
    const std::string &str, std::vector<std::string> &tokens, 
    const std::string &delim = " ", const bool multi_delim = true);


/**
 * Removes all characters from the string belonging to to_erase.
 */
std::string& eraseAll(
    std::string &str, const std::string &to_erase = " \n\t");

/**
 * Returns the string with the whitespace at the beginning and end of the string
 * removed.
 */
std::string trim(
    const std::string& str, const std::string &to_erase = " \t\f\v\n\r");
    
/**
 * Returns the string with all characters belonging to the set {' ', '\t', '\f',
 * '\v', '\r', '\n'} removed.
 */
std::string removeWhiteSpace(const std::string& str);

/**
 * Returns the string in upper case.
 */
std::string toUpperCase(const std::string& str);

/**
 * Returns the string in lower case.
 */
std::string toLowerCase(const std::string& str);

/**
 * Returns true if the string represents a number in ascii where number_base is
 * the base of the number space the string represents.
 */
bool isNumeric(const std::string& input, int number_base = 10);

/**
 * Returns true if all the strings in the vector represent numbers in ascii.
 * @see isNumeric()
 */
bool isNumeric(
    const std::vector<std::string>& array, int number_base = 10);


} // namespace String

#endif // UTILITIES_STRING_UTILS_H
