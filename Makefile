CXX ?= g++
BASE_CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
OPT_CXXFLAGS ?= -O2 -static
CXXFLAGS ?= $(BASE_CXXFLAGS) $(OPT_CXXFLAGS)
LDFLAGS ?=

BIN_DIR := bin
BUILD_DIR := build
APP := $(BIN_DIR)/calculator
TEST_APP := $(BIN_DIR)/calculator_tests
BENCH_APP := $(BIN_DIR)/calculator_benchmark
SRC_DIR := src
TEST_DIR := test
SRC_DIRS := $(SRC_DIR)/app $(SRC_DIR)/core $(SRC_DIR)/math $(SRC_DIR)/matrix $(SRC_DIR)/analysis $(SRC_DIR)/algebra $(SRC_DIR)/symbolic $(SRC_DIR)/script $(SRC_DIR)/numeric $(SRC_DIR)/expression $(SRC_DIR)/runtime $(SRC_DIR)/compat $(SRC_DIR)/cas
INCLUDES := $(addprefix -I,$(SRC_DIRS))
CPPFLAGS += $(INCLUDES) -MMD -MP

MAIN_SRC := $(SRC_DIR)/app/main.cpp
BENCH_SRC := $(SRC_DIR)/app/benchmark.cpp
COMMON_SRCS := $(SRC_DIR)/core/calculator_lifecycle.cpp \
	$(SRC_DIR)/core/calculator_help.cpp \
	$(SRC_DIR)/core/core_helpers.cpp \
	$(SRC_DIR)/core/precise_decimal_parser.cpp \
	$(SRC_DIR)/core/decimal_parser.cpp \
	$(SRC_DIR)/core/exact_and_symbolic_render.cpp \
	$(SRC_DIR)/core/script_runtime.cpp \
	$(SRC_DIR)/core/calculator_commands.cpp \
	$(SRC_DIR)/core/state_persistence.cpp \
	$(SRC_DIR)/math/mymath.cpp \
	$(SRC_DIR)/math/mymath_special_functions.cpp \
	$(SRC_DIR)/matrix/matrix.cpp \
	$(SRC_DIR)/matrix/matrix_expression.cpp \
	$(SRC_DIR)/matrix/matrix_linear_algebra.cpp \
	$(SRC_DIR)/matrix/number_matrix.cpp \
	$(SRC_DIR)/numeric/bigint.cpp \
	$(SRC_DIR)/numeric/rational.cpp \
	$(SRC_DIR)/numeric/decimal.cpp \
	$(SRC_DIR)/numeric/functions.cpp \
	$(SRC_DIR)/numeric/complex.cpp \
	$(SRC_DIR)/numeric/number.cpp \
	$(SRC_DIR)/numeric/conversion.cpp \
	$(SRC_DIR)/expression/expr.cpp \
	$(SRC_DIR)/expression/parser.cpp \
	$(SRC_DIR)/expression/printer.cpp \
	$(SRC_DIR)/expression/evaluator.cpp \
	$(SRC_DIR)/cas/simplify.cpp \
	$(SRC_DIR)/cas/differentiate.cpp \
	$(SRC_DIR)/cas/integrate.cpp \
	$(SRC_DIR)/runtime/value.cpp \
	$(SRC_DIR)/runtime/matrix_value.cpp \
	$(SRC_DIR)/runtime/environment.cpp \
	$(SRC_DIR)/runtime/function_registry.cpp \
	$(SRC_DIR)/runtime/line_executor.cpp \
	$(SRC_DIR)/compat/legacy_commands.cpp \
	$(SRC_DIR)/analysis/function_analysis.cpp \
	$(SRC_DIR)/analysis/multivariable_integrator.cpp \
	$(SRC_DIR)/analysis/ode_solver.cpp \
	$(SRC_DIR)/symbolic/node_parser.cpp \
	$(SRC_DIR)/symbolic/algebra_helpers.cpp \
	$(SRC_DIR)/symbolic/polynomial_helpers.cpp \
	$(SRC_DIR)/symbolic/simplify.cpp \
	$(SRC_DIR)/symbolic/transforms.cpp \
	$(SRC_DIR)/symbolic/symbolic_expression_calculus.cpp \
	$(SRC_DIR)/symbolic/symbolic_expression_transforms.cpp \
	$(SRC_DIR)/algebra/polynomial.cpp \
	$(SRC_DIR)/script/script_parser.cpp
COMMON_HDRS := $(SRC_DIR)/core/calculator.h \
	$(SRC_DIR)/core/calculator_internal_types.h \
	$(SRC_DIR)/math/mymath.h \
	$(SRC_DIR)/math/mymath_internal.h \
	$(SRC_DIR)/matrix/matrix.h \
	$(SRC_DIR)/matrix/matrix_internal.h \
	$(SRC_DIR)/matrix/number_matrix.h \
	$(SRC_DIR)/numeric/bigint.h \
	$(SRC_DIR)/numeric/rational.h \
	$(SRC_DIR)/numeric/decimal.h \
	$(SRC_DIR)/numeric/functions.h \
	$(SRC_DIR)/numeric/complex.h \
	$(SRC_DIR)/numeric/number.h \
	$(SRC_DIR)/numeric/precision_context.h \
	$(SRC_DIR)/numeric/conversion.h \
	$(SRC_DIR)/expression/expr.h \
	$(SRC_DIR)/expression/parser.h \
	$(SRC_DIR)/expression/printer.h \
	$(SRC_DIR)/expression/evaluator.h \
	$(SRC_DIR)/cas/simplify.h \
	$(SRC_DIR)/cas/differentiate.h \
	$(SRC_DIR)/cas/integrate.h \
	$(SRC_DIR)/runtime/value.h \
	$(SRC_DIR)/runtime/matrix_value.h \
	$(SRC_DIR)/runtime/environment.h \
	$(SRC_DIR)/runtime/function_registry.h \
	$(SRC_DIR)/runtime/line_executor.h \
	$(SRC_DIR)/analysis/function_analysis.h \
	$(SRC_DIR)/analysis/multivariable_integrator.h \
	$(SRC_DIR)/analysis/ode_solver.h \
	$(SRC_DIR)/symbolic/symbolic_expression.h \
	$(SRC_DIR)/symbolic/symbolic_expression_internal.h \
	$(SRC_DIR)/algebra/polynomial.h \
	$(SRC_DIR)/script/script_parser.h \
	$(SRC_DIR)/script/script_ast.h
COMMON_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(COMMON_SRCS))
MAIN_OBJ := $(BUILD_DIR)/$(MAIN_SRC:.cpp=.o)
BENCH_OBJ := $(BUILD_DIR)/$(BENCH_SRC:.cpp=.o)
TEST_OBJ := $(BUILD_DIR)/$(TEST_DIR)/tests.o
DEPS := $(MAIN_OBJ:.o=.d) $(BENCH_OBJ:.o=.d) $(TEST_OBJ:.o=.d) $(COMMON_OBJS:.o=.d)

.PHONY: all test script-test benchmark check debug asan ubsan clean

all: $(APP)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(APP): $(MAIN_OBJ) $(COMMON_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $(APP)

$(BENCH_APP): $(BENCH_OBJ) $(COMMON_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $(BENCH_APP)

$(TEST_APP): $(TEST_OBJ) $(COMMON_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $(TEST_APP)

test: $(TEST_APP)
	$(TEST_APP)

script-test: $(APP)
	test/script/run_symbolic_cli_validation.sh

benchmark: $(BENCH_APP)
	$(BENCH_APP)

check: test script-test

debug:
	$(MAKE) OPT_CXXFLAGS="-O0 -g"

asan:
	$(MAKE) OPT_CXXFLAGS="-O1 -g -fsanitize=address,undefined" LDFLAGS="-fsanitize=address,undefined"

ubsan:
	$(MAKE) OPT_CXXFLAGS="-O1 -g -fsanitize=undefined" LDFLAGS="-fsanitize=undefined"

clean:
	rm -rf $(BUILD_DIR) $(APP) $(TEST_APP) $(BENCH_APP)

-include $(DEPS)
