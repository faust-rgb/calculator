import os
import re

with open('src/symbolic/calculator_symbolic_commands.cpp', 'r') as f:
    content = f.read()

# We want to grab the whole handle_symbolic_command function
start_idx = content.find('bool handle_symbolic_command(')
end_idx = content.find('std::string SymbolicModule::execute_args')

if start_idx == -1 or end_idx == -1:
    print("Could not find bounds")
    exit(1)

handle_body = content[start_idx:end_idx]

# Let's extract the vector field commands and integral commands
with open('extracted_body.cpp', 'w') as f:
    f.write(handle_body)

print("Extracted to extracted_body.cpp")
