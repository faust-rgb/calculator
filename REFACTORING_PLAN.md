# Calculator Architecture Refactoring Plan

## Executive Summary

This document outlines a comprehensive refactoring plan to address performance bottlenecks and architectural issues in the calculator's parsing and execution pipeline. The plan is organized into four major phases with clear dependencies and priorities.

---

## Phase 1: Lazy Tokenization & Unified Lexer (Priority: HIGH)

### Problem Statement
The current `CommandParser::peek_token()` triggers full tokenization of the entire input on first call, causing unnecessary memory allocation and processing for long inputs or scripts.

### Solution: On-Demand Token Generation

#### 1.1 Create LazyTokenStream Class

**File:** `src/core/lazy_token_stream.h`

```cpp
class LazyTokenStream {
public:
    explicit LazyTokenStream(std::string_view source);

    // Peek without consuming (generates token if needed)
    const Token& peek();
    const Token& peek(size_t offset);  // Lookahead

    // Consume and return current token
    Token advance();

    // Check if at end
    bool is_at_end() const;

    // Reset to specific position (for backtracking)
    void reset_to(size_t position);

    // Get current position for potential backtrack
    size_t checkpoint() const;

private:
    std::string_view source_;
    std::vector<Token> cache_;      // Cache generated tokens
    size_t cache_pos_ = 0;          // Current position in cache
    size_t source_pos_ = 0;         // Current position in source
    bool end_reached_ = false;

    // Generate next token from source
    Token generate_next();
};
```

**Benefits:**
- Tokens generated only when needed
- Memory usage proportional to lookahead depth, not input length
- Supports arbitrary lookahead without full tokenization

#### 1.2 Refactor BaseParser for Lazy Tokenization

Modify `BaseParser` to use `LazyTokenStream` instead of manual position tracking:

```cpp
class BaseParser {
protected:
    explicit BaseParser(std::string_view source, bool lazy_mode = true);

    LazyTokenStream tokens_;  // Replace manual pos_ tracking

    // Updated methods
    const Token& peek_token() { return tokens_.peek(); }
    Token advance_token() { return tokens_.advance(); }
    size_t save_checkpoint() { return tokens_.checkpoint(); }
    void restore_checkpoint(size_t cp) { tokens_.reset_to(cp); }
};
```

#### 1.3 Update CommandParser

Replace the current eager tokenization:

```cpp
// Before: peek_token() generates ALL tokens
const CommandToken& CommandParser::peek_token() {
    if (!tokens_scanned_) {
        // Generates ALL tokens immediately
        while (true) {
            CommandToken tok = next_token();
            tokens_.push_back(tok);
            if (tok.type == CommandTokenType::kEnd) break;
        }
        tokens_scanned_ = true;
    }
    return tokens_[token_pos_];
}

// After: LazyTokenStream handles on-demand generation
const Token& CommandParser::peek_token() {
    return tokens_.peek();
}
```

---

## Phase 2: Deep AST Integration (Priority: HIGH)

### Problem Statement
`CommandParser` stops at recognizing command types and returns raw expression strings. Expression parsing is deferred to evaluation time, causing repeated parsing in loops.

### Solution: Complete AST Construction at Parse Time

#### 2.1 Extend CommandASTNode with Compiled Expressions

```cpp
struct ExpressionInfo {
    std::string_view text;                    // Original text (for display)
    std::unique_ptr<ExpressionAST> ast;       // Compiled AST
    ExpressionHint hint;                      // Type hint
    ExpressionFeature features;               // Features
    bool is_compiled = false;

    // Slot bindings for fast variable lookup
    std::vector<std::pair<std::string, int>> variable_slots;
};

struct AssignmentInfo {
    std::string_view variable;
    ExpressionInfo expression;  // Was: std::string_view
};

struct FunctionCallInfo {
    std::string_view name;
    std::vector<ExpressionInfo> arguments;  // Was: vector<string_view>
};
```

#### 2.2 Integrate Expression Compiler into CommandParser

```cpp
class CommandParser : public BaseParser {
public:
    CommandASTNode parse() {
        // Parse command structure
        CommandASTNode node = parse_command();

        // Compile all expressions within
        compile_expressions(node);

        return node;
    }

private:
    void compile_expressions(CommandASTNode& node) {
        switch (node.kind) {
            case CommandKind::kAssignment: {
                auto* info = node.as_assignment();
                compile_expression(info->expression);
                break;
            }
            case CommandKind::kFunctionCall: {
                auto* info = node.as_function_call();
                for (auto& arg : info->arguments) {
                    compile_expression(arg);
                }
                break;
            }
            case CommandKind::kExpression: {
                auto* expr = node.as_expression();
                compile_expression(*expr);
                break;
            }
            // ... other cases
        }
    }

    void compile_expression(ExpressionInfo& info) {
        info.ast = compile_expression_ast(std::string(info.text));
        if (info.ast) {
            info.is_compiled = true;
            info.hint = analyze_expression_hint(std::string(info.text));
            info.features = analyze_expression_features(std::string(info.text));
        }
    }
};
```

#### 2.3 Pre-compile Script Statements

Modify `ScriptParser` to compile expressions during parsing:

```cpp
// In script_parser.cpp
std::unique_ptr<IfStatement> ScriptParser::parse_if_statement() {
    auto stmt = std::make_unique<IfStatement>();

    // Parse condition
    stmt->condition = parse_expression_text();

    // PRE-COMPILE the condition
    stmt->cache = std::make_shared<ExpressionCache>();
    stmt->cache->expanded = stmt->condition;
    stmt->cache->compiled_ast = compile_expression_ast(stmt->condition);
    stmt->cache->is_compiled = (stmt->cache->compiled_ast != nullptr);
    stmt->cache->hint = analyze_expression_hint(stmt->condition);
    stmt->cache->features = analyze_expression_features(stmt->condition);

    // Parse branches
    stmt->then_branch = parse_statement();
    // ...

    return stmt;
}
```

---

## Phase 3: Smart Lookahead & Backtracking (Priority: MEDIUM)

### Problem Statement
Current backtracking in `parse_function_call` uses `token_pos_ = 0` to reset the entire parser state, which is inefficient and error-prone.

### Solution: Structured Parsing with Checkpoints

#### 3.1 Implement Checkpoint-based Backtracking

```cpp
class CommandParser : public BaseParser {
private:
    CommandASTNode parse_definition_or_assignment(const Token& id_tok) {
        // Save checkpoint before attempting to parse as function definition
        size_t checkpoint = save_checkpoint();

        advance_token(); // consume identifier

        if (peek_token().kind == TokenKind::kLParen) {
            // Try to parse as function definition
            if (auto def = try_parse_function_definition(id_tok)) {
                return *def;
            }

            // Backtrack to checkpoint and try function call
            restore_checkpoint(checkpoint);
            return parse_function_call(id_tok);
        }

        // ... handle assignment or expression
    }

    std::optional<CommandASTNode> try_parse_function_definition(const Token& id_tok) {
        // Parse parameter list
        advance_token(); // consume '('

        std::vector<std::string_view> params;
        while (peek_token().kind == TokenKind::kIdentifier) {
            params.push_back(advance_token().text);
            if (!match_token(TokenKind::kComma)) break;
        }

        if (!match_token(TokenKind::kRParen)) {
            return std::nullopt;  // Not a valid parameter list
        }

        if (!match_token(TokenKind::kEqual)) {
            return std::nullopt;  // Not a function definition
        }

        // Success - parse body
        std::string_view body = parse_remaining_as_expression();
        return CommandASTNode::make_function_definition(id_tok.text, params, body);
    }
};
```

#### 3.2 Replace Full Reset with Local Backtrack

```cpp
// Before: Brute-force reset
CommandASTNode CommandParser::parse_function_call(const Token& id_tok) {
    // ... parse arguments ...

    if (peek_token().type != CommandTokenType::kEnd) {
        token_pos_ = 0;  // FULL RESET - inefficient!
        return parse_expression();
    }
    // ...
}

// After: Local backtrack
CommandASTNode CommandParser::parse_function_call(const Token& id_tok) {
    size_t start_checkpoint = save_checkpoint();

    // ... parse arguments ...

    if (peek_token().kind != TokenKind::kEnd) {
        // Backtrack to start and parse as full expression
        restore_checkpoint(start_checkpoint);
        return parse_expression();
    }
    // ...
}
```

---

## Phase 4: Unified Scope Model (Priority: HIGH)

### Problem Statement
`script_runtime.cpp` contains duplicated code paths for `FlatScopeStack` and `std::map<std::string, StoredValue>` scope management, with `if (impl->use_flat_scopes)` checks throughout.

### Solution: Complete Migration to FlatScopeStack

#### 4.1 Remove Legacy Scope Path

1. **Delete `local_scopes` member from `Calculator::Impl`:**
```cpp
struct Calculator::Impl {
    // REMOVE: std::vector<std::map<std::string, StoredValue>> local_scopes;

    FlatScopeStack flat_scopes;  // Only scope storage
    // REMOVE: bool use_flat_scopes = true;  // No longer needed
};
```

2. **Simplify `visible_variables`:**
```cpp
VariableResolver visible_variables(const Calculator::Impl* impl) {
    return VariableResolver(&impl->variables, &impl->flat_scopes);
}
```

3. **Simplify `assign_visible_variable`:**
```cpp
void assign_visible_variable(Calculator::Impl* impl,
                             const std::string& name,
                             const StoredValue& value) {
    // Try flat_scopes first
    if (VariableSlot* existing = impl->flat_scopes.find(name)) {
        existing->value = value;
        return;
    }

    // Check global variables
    if (impl->variables.find(name) != impl->variables.end()) {
        impl->variables[name] = value;
        return;
    }

    // Create in current scope
    if (impl->flat_scopes.scope_depth() > 0) {
        impl->flat_scopes.set(name, value);
    } else {
        impl->variables[name] = value;
    }
}
```

#### 4.2 Simplify Script Execution

Remove all `if (impl->use_flat_scopes)` branches:

```cpp
// Before: Duplicated paths
ScriptSignal execute_script_statement(...) {
    case script::Statement::Kind::kFor: {
        if (impl->use_flat_scopes) {
            // FlatScopeStack path (~50 lines)
        } else {
            // std::map path (~50 lines, nearly identical)
        }
    }
}

// After: Single path
ScriptSignal execute_script_statement(...) {
    case script::Statement::Kind::kFor: {
        impl->flat_scopes.push_scope();
        // Single implementation (~50 lines)
        impl->flat_scopes.pop_scope();
    }
}
```

**Estimated code reduction:** ~400 lines in `script_runtime.cpp`

---

## Phase 5: Enhanced Syntax Validation (Priority: MEDIUM)

### Problem Statement
`has_obvious_syntax_error` only checks bracket balance and trailing operators, missing many common errors.

### Solution: Comprehensive Syntax Validation

#### 5.1 Create SyntaxValidator Class

```cpp
class SyntaxValidator {
public:
    struct Error {
        std::string message;
        size_t position;
        Severity severity;  // Error, Warning
    };

    std::vector<Error> validate(std::string_view expression);

private:
    bool check_bracket_balance(std::string_view expr, std::vector<Error>& errors);
    bool check_operator_sequences(std::string_view expr, std::vector<Error>& errors);
    bool check_operand_context(std::string_view expr, std::vector<Error>& errors);
    bool check_function_syntax(std::string_view expr, std::vector<Error>& errors);
};
```

#### 5.2 Implement Validation Rules

```cpp
std::vector<Error> SyntaxValidator::validate(std::string_view expr) {
    std::vector<Error> errors;

    // 1. Bracket balance (existing)
    check_bracket_balance(expr, errors);

    // 2. Invalid operator sequences: ++, **, //, ^^
    check_operator_sequences(expr, errors);

    // 3. Missing operands: binary op at start/end
    check_operand_context(expr, errors);

    // 4. Function call syntax: identifier followed by (
    check_function_syntax(expr, errors);

    // 5. String literal termination
    check_string_termination(expr, errors);

    // 6. Invalid character sequences
    check_invalid_chars(expr, errors);

    return errors;
}
```

---

## Phase 6: Unified Parser Factory (Priority: MEDIUM)

### Problem Statement
`evaluate_expression_value` tries multiple parsers sequentially (exact, matrix, symbolic, scalar), causing wasted work when the wrong path is attempted.

### Solution: Context-Aware Single Dispatch

#### 6.1 Create UnifiedParserFactory

```cpp
class UnifiedParserFactory {
public:
    struct ParseContext {
        bool exact_mode = false;
        bool allow_matrix = true;
        bool allow_complex = true;
        bool need_symbolic = false;
        const std::set<std::string>* known_functions = nullptr;
    };

    // Analyze once and return appropriate parser
    Parser select_parser(const std::string& expr, const ParseContext& ctx);

private:
    // Token-based analysis (single pass)
    Parser analyze_and_select(const std::vector<Token>& tokens, const ParseContext& ctx);
};

enum class Parser {
    kStringLiteral,
    kIdentifier,
    kExact,
    kMatrix,
    kComplex,
    kScalar,
    kSymbolic,
};
```

#### 6.2 Single-Pass Expression Analysis

```cpp
Parser UnifiedParserFactory::select_parser(const std::string& expr, const ParseContext& ctx) {
    // Single tokenization pass
    ExpressionLexer lexer(expr);
    std::vector<Token> tokens;
    lexer.tokenize(&tokens);

    // Analyze tokens to determine parser
    return analyze_and_select(tokens, ctx);
}

Parser UnifiedParserFactory::analyze_and_select(const std::vector<Token>& tokens,
                                                 const ParseContext& ctx) {
    bool has_i = false;
    bool has_bracket = false;
    bool has_matrix_func = false;
    bool has_rat_call = false;
    int paren_depth = 0;

    for (const auto& tok : tokens) {
        // Single pass analysis
        if (tok.kind == TokenKind::kIdentifier) {
            if (tok.text == "i") has_i = true;
            if (tok.text == "rat") has_rat_call = true;
            if (is_matrix_function(tok.text)) has_matrix_func = true;
        }
        if (tok.kind == TokenKind::kLBracket) has_bracket = true;
        // ...
    }

    // Dispatch decision
    if (ctx.exact_mode && !has_i && !has_bracket) return Parser::kExact;
    if (has_bracket || has_matrix_func) return Parser::kMatrix;
    if (has_i) return Parser::kComplex;
    if (has_rat_call) return Parser::kScalar;  // Special rat handling

    return Parser::kScalar;
}
```

#### 6.3 Refactor evaluate_expression_value

```cpp
StoredValue evaluate_expression_value(Calculator* calculator,
                                      Calculator::Impl* impl,
                                      const std::string& expression,
                                      bool exact_mode,
                                      std::shared_ptr<ExpressionCache>* cache) {
    // Get or create cache
    auto expr_cache = get_or_create_cache(cache, calculator, expression);

    // Single dispatch decision
    UnifiedParserFactory::ParseContext ctx;
    ctx.exact_mode = exact_mode;
    ctx.need_symbolic = impl->symbolic_constants_mode;

    Parser parser = UnifiedParserFactory().select_parser(expr_cache->expanded, ctx);

    // Direct dispatch - no fallback attempts
    switch (parser) {
        case Parser::kStringLiteral:
            return parse_string_literal(expr_cache->expanded);
        case Parser::kIdentifier:
            return lookup_variable(expr_cache->expanded, impl);
        case Parser::kExact:
            return evaluate_exact(expr_cache->expanded, impl);
        case Parser::kMatrix:
            return evaluate_matrix(expr_cache->expanded, impl);
        case Parser::kComplex:
            return evaluate_complex(expr_cache->expanded, impl);
        case Parser::kScalar:
            return evaluate_scalar(expr_cache->expanded, impl, expr_cache);
        case Parser::kSymbolic:
            return evaluate_symbolic(expr_cache->expanded, impl);
    }
}
```

---

## Phase 7: On-Demand Symbolic Computation (Priority: LOW)

### Problem Statement
`attach_symbolic_text` is called after every scalar calculation, even when symbolic results aren't needed.

### Solution: Lazy Symbolic Attachment

#### 7.1 Defer Symbolic Processing

```cpp
struct StoredValue {
    double decimal = 0.0;
    // ... other fields ...

    // Symbolic text is now computed on demand
    mutable bool has_symbolic_text = false;
    mutable std::string symbolic_text;
    mutable bool symbolic_computed = false;

    // Lazy getter
    std::string get_symbolic_text(const std::string& expr) const {
        if (!symbolic_computed) {
            symbolic_text = compute_symbolic_text(expr);
            symbolic_computed = true;
            has_symbolic_text = !symbolic_text.empty();
        }
        return symbolic_text;
    }
};
```

#### 7.2 Conditional Symbolic Mode

```cpp
// Only compute symbolic text when:
// 1. symbolic_constants_mode is enabled, AND
// 2. Result is being displayed (not intermediate), AND
// 3. Expression contains symbolic constants

StoredValue evaluate_scalar(const std::string& expr, ...) {
    StoredValue result;
    result.decimal = compute_value(expr);

    // Don't compute symbolic text here
    // It will be computed on demand when formatting output

    return result;
}

std::string format_stored_value(const StoredValue& value,
                                bool symbolic_mode,
                                const std::string& expr = "") {
    if (symbolic_mode && !value.symbolic_computed && !expr.empty()) {
        value.get_symbolic_text(expr);
    }

    if (value.has_symbolic_text) {
        return value.symbolic_text;
    }
    // ... format numeric value
}
```

---

## Phase 8: Command Registry for Module Self-Registration (Priority: MEDIUM)

### Problem Statement
`try_process_function_command` rebuilds `CoreServices` for each command, and module registration is scattered.

### Solution: Centralized Command Registry

#### 8.1 Create CommandRegistry

```cpp
class CommandRegistry {
public:
    using Handler = std::function<bool(const std::string&, std::string*, bool)>;

    // Register a command handler
    void register_command(const std::string& name,
                          Handler handler,
                          const std::string& help_text = "");

    // Try to process a command
    bool try_process(const std::string& input,
                     std::string* output,
                     bool exact_mode);

    // Get all registered commands (for completion)
    std::vector<std::string> get_commands() const;

    // Get help for a command
    std::string get_help(const std::string& name) const;

private:
    std::map<std::string, Handler> handlers_;
    std::map<std::string, std::string> help_texts_;
};
```

#### 8.2 Module Self-Registration

```cpp
// In each module's initialization
class MatrixModule : public CalculatorModule {
public:
    void register_commands(CommandRegistry& registry) override {
        registry.register_command("mat", [this](const std::string& args, std::string* out, bool exact) {
            return this->handle_mat(args, out, exact);
        }, "Create a matrix: mat([1,2,3]) or mat([[1,2],[3,4]])");

        registry.register_command("transpose", [this](...) { ... });
        registry.register_command("inverse", [this](...) { ... });
        // ... more commands
    }
};
```

#### 8.3 Integrate with Calculator

```cpp
struct Calculator::Impl {
    CommandRegistry command_registry;  // Centralized registry
    // ... other members
};

// In calculator initialization
void Calculator::initialize() {
    impl_->command_registry = CommandRegistry();

    for (auto& module : impl_->registered_modules) {
        module->register_commands(impl_->command_registry);
    }
}

// Simplified command processing
bool Calculator::try_process_function_command(const std::string& input,
                                               std::string* output,
                                               bool exact_mode) {
    return impl_->command_registry.try_process(input, output, exact_mode);
}
```

---

## Phase 9: Script Function Return Type Enhancement (Priority: LOW)

### Problem Statement
`invoke_script_function_decimal` forces script functions to return `double`, preventing them from returning matrices, complex numbers, or strings.

### Solution: Return StoredValue from Script Functions

#### 9.1 Create Generic Script Function Invocation

```cpp
StoredValue invoke_script_function(Calculator* calculator,
                                   Calculator::Impl* impl,
                                   const std::string& name,
                                   const std::vector<StoredValue>& arguments) {
    auto it = impl->script_functions.find(name);
    if (it == impl->script_functions.end()) {
        throw std::runtime_error("unknown function: " + name);
    }

    const ScriptFunction& function = it->second;

    // Create scope with arguments
    impl->flat_scopes.push_scope();
    for (size_t i = 0; i < arguments.size(); ++i) {
        impl->flat_scopes.set(function.parameter_names[i], arguments[i]);
    }

    ScopeGuard guard([&]() { impl->flat_scopes.pop_scope(); });

    // Execute function body
    std::string ignored_output;
    ScriptSignal signal = execute_script_block(
        calculator, impl, *function.body, false, &ignored_output, false);

    if (signal.kind != ScriptSignal::Kind::kReturn || !signal.has_value) {
        throw std::runtime_error("script function must return a value");
    }

    return signal.value;  // Return full StoredValue
}
```

#### 9.2 Update Expression Evaluator to Handle Typed Returns

```cpp
// In evaluate_compiled_ast
double evaluate_compiled_ast(...) {
    // ...
    if (node->is_script_function) {
        // Get full result
        StoredValue result = invoke_script_function(...);

        // Extract scalar or throw
        if (result.is_matrix) {
            throw std::runtime_error("cannot use matrix result in scalar expression");
        }
        if (result.is_complex) {
            throw std::runtime_error("cannot use complex result in scalar expression");
        }
        return result.decimal;
    }
    // ...
}
```

---

## Implementation Order & Dependencies

```
Phase 1 (Lazy Tokenization)
    ↓
Phase 3 (Smart Lookahead) ← depends on Phase 1
    ↓
Phase 2 (Deep AST) ← depends on Phase 3
    ↓
Phase 4 (Unified Scope) ← independent, can parallelize
    ↓
Phase 5 (Syntax Validation) ← independent
    ↓
Phase 6 (Unified Parser) ← depends on Phase 2
    ↓
Phase 7 (Lazy Symbolic) ← depends on Phase 6
    ↓
Phase 8 (Command Registry) ← independent
    ↓
Phase 9 (Typed Returns) ← depends on Phase 4
```

**Recommended execution order:**
1. **Phase 4** (Unified Scope) - Largest immediate impact, removes ~400 lines of duplicated code
2. **Phase 1** (Lazy Tokenization) - Foundation for parser improvements
3. **Phase 3** (Smart Lookahead) - Depends on Phase 1
4. **Phase 2** (Deep AST) - Major performance win for loops
5. **Phase 6** (Unified Parser) - Clean up evaluation dispatch
6. **Phase 8** (Command Registry) - Architectural cleanup
7. **Phase 5** (Syntax Validation) - User experience improvement
8. **Phase 7** (Lazy Symbolic) - Minor performance optimization
9. **Phase 9** (Typed Returns) - Feature enhancement

---

## Expected Outcomes

### Performance Improvements
- **Loop execution:** 2-5x faster (pre-compiled expressions, no string parsing in loops)
- **Memory usage:** 30-50% reduction (lazy tokenization, unified scope)
- **Parse time:** 20-40% faster (single-pass analysis, no fallback attempts)

### Code Quality
- **Lines of code:** ~500-700 lines removed (duplicated scope paths, fallback logic)
- **Complexity:** Reduced cyclomatic complexity in `evaluate_expression_value` and `execute_script_statement`
- **Maintainability:** Single scope model, unified parser dispatch

### User Experience
- **Error messages:** More precise syntax error reporting with position information
- **Script capabilities:** Functions can return matrices and complex numbers
- **Responsiveness:** Faster startup for large scripts (lazy parsing)

---

## Testing Strategy

### Unit Tests
1. `LazyTokenStream` - verify on-demand generation, lookahead, backtracking
2. `FlatScopeStack` - verify scope push/pop, variable lookup, slot binding
3. `UnifiedParserFactory` - verify correct parser selection for each expression type
4. `SyntaxValidator` - verify all error detection rules

### Integration Tests
1. Script execution benchmarks (before/after timing)
2. Memory usage profiling (valgrind massif)
3. Expression evaluation correctness (existing test suite)

### Regression Tests
- All existing calculator tests must pass
- Script function tests with various return types
- Error message quality tests

---

## Risk Assessment

### High Risk
- **Phase 4 (Unified Scope):** Touches all script execution code. Mitigation: comprehensive test coverage, incremental rollout with feature flag.

### Medium Risk
- **Phase 2 (Deep AST):** Changes parsing pipeline. Mitigation: fallback to string parsing if AST compilation fails.
- **Phase 6 (Unified Parser):** Changes evaluation dispatch. Mitigation: maintain existing dispatch as fallback initially.

### Low Risk
- **Phase 1 (Lazy Tokenization):** Internal change, same external behavior.
- **Phase 5 (Syntax Validation):** Additive, doesn't change parsing.
- **Phase 7 (Lazy Symbolic):** Additive, doesn't change computation.
- **Phase 8 (Command Registry):** Refactoring with same behavior.

---

## Migration Path

### Stage 1: Preparation (Week 1)
- Add `LazyTokenStream` class
- Add `CommandRegistry` class
- Add `SyntaxValidator` class
- Add feature flag for unified scope

### Stage 2: Core Changes (Weeks 2-3)
- Implement Phase 4 (Unified Scope) behind feature flag
- Implement Phase 1 (Lazy Tokenization)
- Implement Phase 3 (Smart Lookahead)

### Stage 3: Integration (Weeks 4-5)
- Implement Phase 2 (Deep AST)
- Implement Phase 6 (Unified Parser)
- Implement Phase 8 (Command Registry)

### Stage 4: Polish (Week 6)
- Implement Phase 5 (Syntax Validation)
- Implement Phase 7 (Lazy Symbolic)
- Implement Phase 9 (Typed Returns)
- Remove feature flags and legacy code

### Stage 5: Validation (Week 7)
- Performance benchmarking
- Memory profiling
- Documentation update
- Release
