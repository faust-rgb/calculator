# Commands Reference

## General

- `help`
- `:help`
- `:help commands`
- `:help functions`
- `:help matrix`
- `:help examples`
- `:help exact`
- `:help variables`
- `:help persistence`
- `:help programmer`
- `exit`
- `quit`

## Exact Mode

- `:exact on`
  Enable exact fraction mode
- `:exact off`
  Disable exact fraction mode
- `:exact`
  Show current exact mode status

## Symbolic Constants Mode

- `:symbolic on`
  Preserve `pi` and `e` in scalar display results
- `:symbolic off`
  Return to normal numeric display
- `:symbolic`
  Show current symbolic constants mode status

## Programmer Formatting

- `:hexprefix on`
  Show hex results with a `0x` prefix
- `:hexprefix off`
  Hide the `0x` prefix
- `:hexprefix`
  Show current hex prefix mode status
- `:hexcase upper`
  Use uppercase hex digits such as `0xFF`
- `:hexcase lower`
  Use lowercase hex digits such as `0xff`
- `:hexcase`
  Show current hex letter-case mode

## Variables

- `name = expression`
  Assign a variable
- `:vars`
  List all stored variables
- `:clear name`
  Remove a single variable
- `:clear`
  Remove all variables

## History

- up arrow
  Recall previous command in interactive mode
- down arrow
  Move toward newer history entries
- `:history`
  Print the current session history

## Persistence

- `:save file`
  Save variables to a file
- `:load file`
  Load variables from a file
- `:run file`
  Execute a script file

Notes:

- current persistence supports scalar variables
- matrix variables can be used in-session, but are not yet saved or restored

## Scripting

- `fn name(args) { ... }`
  Define a script function
- `if (cond) { ... } else { ... }`
  Conditional execution
- `while (cond) { ... }`
  Loop while condition is true
- `for (init; cond; step) { ... }`
  C-style loop
- `return expr;`
  Return from a script function
- `break;`
  Exit the current loop
- `continue;`
  Skip to the next loop iteration
- string literals such as `"hello"`
  Script string value
- `print(a, b, c);`
  Emit one or more formatted values inside a script

## Autocomplete

- `Tab`
  Autocomplete commands, functions, variables, and custom functions
- double `Tab`
  Show the current candidate list when multiple completions match

Examples:

- type `:he` then press `Tab` -> `:help`
- type `sq` then press `Tab` -> `sqrt(`
- type `:help ma` then press `Tab` -> `:help matrix`
- type `:help ` then double `Tab` -> show help topics
- type `g` inside `diff(g` then press `Tab` -> complete a matching custom function

## Prefixed Integer Literals

These are expression features rather than commands, but they are commonly used
like shell-style numeric shortcuts:

- `0b1010`
- `0o77`
- `0xFF`

They can be mixed directly into expressions:

- `0b1010 + 0xF`
- `and(0xF, 0b1010)`
