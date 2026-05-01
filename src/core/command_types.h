#ifndef CALCULATOR_COMMAND_TYPES_H
#define CALCULATOR_COMMAND_TYPES_H

#include <string>
#include <string_view>
#include <utility>

enum class CommandSyntax {
    kCall,
    kMeta
};

struct CommandKey {
    CommandSyntax syntax = CommandSyntax::kCall;
    std::string name;

    bool operator<(const CommandKey& other) const {
        if (syntax != other.syntax) {
            return static_cast<int>(syntax) < static_cast<int>(other.syntax);
        }
        return name < other.name;
    }

    bool operator==(const CommandKey& other) const {
        return syntax == other.syntax && name == other.name;
    }
};

inline CommandKey call_command_key(std::string_view name) {
    return {CommandSyntax::kCall, std::string(name)};
}

inline CommandKey meta_command_key(std::string_view name) {
    return {CommandSyntax::kMeta, std::string(name)};
}

inline std::string command_key_display(const CommandKey& key) {
    return key.syntax == CommandSyntax::kMeta ? ":" + key.name : key.name;
}

#endif
