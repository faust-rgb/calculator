#ifndef PARSER_PARSER_UTILS_H
#define PARSER_PARSER_UTILS_H

#include <string>
#include <string_view>
#include <vector>

namespace parser_utils {

/**
 * @struct BalancedScanResult
 * @brief Results of a balanced scan (brackets, parens, quotes).
 */
struct BalancedScanResult {
    bool balanced = true;
    std::size_t first_mismatch_pos = std::string_view::npos;
    int paren_depth = 0;
    int bracket_depth = 0;
    int brace_depth = 0;
};

/**
 * @brief Scans a string to find top-level delimiter positions while respecting brackets and quotes.
 */
BalancedScanResult scan_balanced(std::string_view text);

/**
 * @brief Checks if a string is wrapped by balanced delimiters and can be stripped.
 */
bool is_wrapped_by(std::string_view text, char open, char close);

/**
 * @brief Finds the position of a top-level character (not inside brackets/quotes).
 */
std::size_t find_top_level(std::string_view text, char target);

/**
 * @brief Splits a string by a delimiter, only at the top nesting level.
 */
std::vector<std::string> split_top_level(std::string_view text, char delimiter);

/**
 * @brief Checks if a string contains script-specific syntax (like curly braces).
 */
bool contains_script_syntax(std::string_view text);

} // namespace parser_utils

#endif
