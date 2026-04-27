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

## 2.0 Precision Runtime

- `:v2 expression`
  Force evaluation with the 2.0 exact/high-precision runtime
- `:precision digits`
  Set 2.0 decimal precision
- `:v2precision digits`
  Compatibility alias for `:precision`
- `:v2vars`
  List only 2.0 runtime variables
- `:v2clear`
  Clear only 2.0 runtime variables

Notes:

- ordinary `:vars`, `:clear name`, and `:clear` include both 1.0 and 2.0
  variables
- migrated 2.0 forms route to 2.0 without the prefix, including complex
  literals using `i`, large integer arithmetic, nested matrix literals, and
  2.0 CAS commands such as `simplify(...)`

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
- `:run file.calc`
  Execute a `.calc` script file

Notes:

- current persistence writes `STATE_V5`
- state files preserve scalars, strings, matrices, custom functions, script
  functions, 2.0 precision, and 2.0 variables
- older `STATE_V1` through `STATE_V4` files remain loadable

## Scripting

Use `bin/calculator file.calc` or `:run file.calc` to execute scripts. The
dedicated syntax guide is `test/script/SYNTAX_GUIDE.md`.

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
- `rol(1, 3)`
- `popcount(0xF0)`
