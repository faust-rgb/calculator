CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -O2 -static

BIN_DIR := bin
APP := $(BIN_DIR)/calculator
TEST_APP := $(BIN_DIR)/calculator_tests
SRC_DIR := src
TEST_DIR := test
SRC_DIRS := $(SRC_DIR)/app $(SRC_DIR)/core $(SRC_DIR)/math $(SRC_DIR)/matrix $(SRC_DIR)/analysis $(SRC_DIR)/algebra $(SRC_DIR)/symbolic $(SRC_DIR)/script
INCLUDES := $(addprefix -I,$(SRC_DIRS))

MAIN_SRC := $(SRC_DIR)/app/main.cpp
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
	$(SRC_DIR)/matrix/matrix_linear_algebra.cpp \
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
	$(SRC_DIR)/analysis/function_analysis.h \
	$(SRC_DIR)/analysis/multivariable_integrator.h \
	$(SRC_DIR)/analysis/ode_solver.h \
	$(SRC_DIR)/symbolic/symbolic_expression.h \
	$(SRC_DIR)/symbolic/symbolic_expression_internal.h \
	$(SRC_DIR)/algebra/polynomial.h \
	$(SRC_DIR)/script/script_parser.h \
	$(SRC_DIR)/script/script_ast.h
COMMON_PARTS :=

.PHONY: all test clean

all: $(APP)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(APP): $(MAIN_SRC) $(COMMON_SRCS) $(COMMON_HDRS) $(COMMON_PARTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MAIN_SRC) $(COMMON_SRCS) -o $(APP)

$(TEST_APP): $(TEST_DIR)/tests.cpp $(COMMON_SRCS) $(COMMON_HDRS) $(COMMON_PARTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEST_DIR)/tests.cpp $(COMMON_SRCS) -o $(TEST_APP)

test: $(TEST_APP)
	$(TEST_APP)

clean:
	rm -f $(APP) $(TEST_APP)
