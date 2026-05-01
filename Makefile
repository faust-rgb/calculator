CXX ?= g++
BASE_CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
OPT_CXXFLAGS ?= -O0 -g -static
CXXFLAGS ?= $(BASE_CXXFLAGS) $(OPT_CXXFLAGS)
LDFLAGS ?=

BIN_DIR := bin
BUILD_DIR := build
APP := $(BIN_DIR)/calculator
TEST_APP := $(BIN_DIR)/calculator_tests
SRC_DIR := src
TEST_DIR := test
TEST_SUITE_DIR := $(TEST_DIR)/suites
SRC_DIRS := $(SRC_DIR)/app $(SRC_DIR)/core $(SRC_DIR)/math $(SRC_DIR)/matrix $(SRC_DIR)/analysis $(SRC_DIR)/polynomial $(SRC_DIR)/symbolic $(SRC_DIR)/script $(SRC_DIR)/statistics $(SRC_DIR)/dsp $(SRC_DIR)/plot $(SRC_DIR)/types $(SRC_DIR)/precise $(SRC_DIR)/parser $(SRC_DIR)/math/helpers
INCLUDES := -I$(SRC_DIR) $(addprefix -I,$(SRC_DIRS)) -I$(TEST_DIR)
CPPFLAGS += $(INCLUDES) -MMD -MP

MAIN_SRC := $(SRC_DIR)/app/main.cpp
COMMON_SRCS := $(SRC_DIR)/core/calculator_lifecycle.cpp \
	$(SRC_DIR)/script/script_signal.cpp \
	$(SRC_DIR)/core/variable_resolver.cpp \
	$(SRC_DIR)/core/utils.cpp \
	$(SRC_DIR)/core/calculator_basic_commands.cpp \
	$(SRC_DIR)/core/calculator_help.cpp \
	$(SRC_DIR)/core/system_module.cpp \
	$(SRC_DIR)/core/expression_compiler.cpp \
	$(SRC_DIR)/core/calculator_service_factory.cpp \
	$(SRC_DIR)/core/module_registration.cpp \
		$(SRC_DIR)/core/command_parser.cpp \
		$(SRC_DIR)/precise/rational.cpp \
		$(SRC_DIR)/precise/precise_decimal.cpp \
		$(SRC_DIR)/precise/precise_parser.cpp \
		$(SRC_DIR)/precise/precise_module.cpp \
	$(SRC_DIR)/parser/decimal_parser.cpp \
	$(SRC_DIR)/parser/exact_parser.cpp \
	$(SRC_DIR)/parser/symbolic_render_parser.cpp \
	$(SRC_DIR)/math/helpers/integer_helpers.cpp \
	$(SRC_DIR)/math/helpers/combinatorics.cpp \
	$(SRC_DIR)/math/helpers/bitwise_helpers.cpp \
	$(SRC_DIR)/math/helpers/unit_conversions.cpp \
	$(SRC_DIR)/math/helpers/base_conversions.cpp \
	$(SRC_DIR)/script/script_runtime.cpp \
	$(SRC_DIR)/core/calculator_commands.cpp \
	$(SRC_DIR)/analysis/calculator_simplex.cpp \
	$(SRC_DIR)/core/calculator_state_persistence.cpp \
	$(SRC_DIR)/math/mymath.cpp \
	$(SRC_DIR)/math/mymath_special_functions.cpp \
	$(SRC_DIR)/math/standard_math_module.cpp \
	$(SRC_DIR)/math/integer_math_module.cpp \
	$(SRC_DIR)/matrix/matrix.cpp \
	$(SRC_DIR)/matrix/matrix_module.cpp \
	$(SRC_DIR)/matrix/matrix_ops.cpp \
	$(SRC_DIR)/matrix/matrix_transform.cpp \
	$(SRC_DIR)/matrix/matrix_utility.cpp \
	$(SRC_DIR)/matrix/matrix_solvers.cpp \
	$(SRC_DIR)/matrix/matrix_stats.cpp \
	$(SRC_DIR)/matrix/matrix_interpolation.cpp \
	$(SRC_DIR)/matrix/matrix_expression.cpp \
	$(SRC_DIR)/matrix/matrix_linear_algebra.cpp \
	$(SRC_DIR)/matrix/calculator_matrix_commands.cpp \
	$(SRC_DIR)/analysis/function_analysis.cpp \
	$(SRC_DIR)/analysis/multivariable_integrator.cpp \
	$(SRC_DIR)/analysis/multidim_integration.cpp \
	$(SRC_DIR)/analysis/ode_solver.cpp \
	$(SRC_DIR)/analysis/calculator_ode.cpp \
	$(SRC_DIR)/analysis/calculator_integration.cpp \
	$(SRC_DIR)/analysis/vector_field_theorems.cpp \
	$(SRC_DIR)/analysis/calculator_analysis_cmds.cpp \
	$(SRC_DIR)/analysis/calculator_rootfinding.cpp \
	$(SRC_DIR)/analysis/calculator_optimization.cpp \
	$(SRC_DIR)/analysis/calculator_series.cpp \
	$(SRC_DIR)/analysis/optimization_helpers.cpp \
	$(SRC_DIR)/symbolic/node_parser.cpp \
	$(SRC_DIR)/symbolic/algebra_helpers.cpp \
	$(SRC_DIR)/symbolic/polynomial_helpers.cpp \
	$(SRC_DIR)/symbolic/simplify.cpp \
	$(SRC_DIR)/symbolic/transforms.cpp \
	$(SRC_DIR)/symbolic/symbolic_expression_calculus.cpp \
	$(SRC_DIR)/symbolic/symbolic_expression_transforms.cpp \
	$(SRC_DIR)/symbolic/symbolic_polynomial.cpp \
	$(SRC_DIR)/symbolic/integration_engine.cpp \
	$(SRC_DIR)/symbolic/calculator_transforms.cpp \
	$(SRC_DIR)/symbolic/calculator_symbolic_commands.cpp \
	$(SRC_DIR)/polynomial/calculator_polynomial.cpp \
	$(SRC_DIR)/polynomial/polynomial.cpp \
	$(SRC_DIR)/script/script_parser.cpp \
	$(SRC_DIR)/statistics/statistics.cpp \
	$(SRC_DIR)/statistics/probability.cpp \
	$(SRC_DIR)/statistics/calculator_statistics.cpp \
	$(SRC_DIR)/dsp/fft.cpp \
	$(SRC_DIR)/dsp/convolution.cpp \
	$(SRC_DIR)/dsp/window_functions.cpp \
	$(SRC_DIR)/dsp/filter_design.cpp \
	$(SRC_DIR)/dsp/time_frequency.cpp \
	$(SRC_DIR)/dsp/residue.cpp \
		$(SRC_DIR)/dsp/dsp_module.cpp \
	$(SRC_DIR)/dsp/calculator_signal_commands.cpp \
	$(SRC_DIR)/plot/plot_renderer.cpp \
	$(SRC_DIR)/plot/calculator_plot.cpp \
	$(SRC_DIR)/plot/plot_styles.cpp \
	$(SRC_DIR)/plot/svg_renderer.cpp
COMMON_HDRS := $(SRC_DIR)/core/calculator.h \
	$(SRC_DIR)/core/calculator_internal_types.h \
	$(SRC_DIR)/core/system_module.h \
		$(SRC_DIR)/analysis/calculator_simplex.h \
	$(SRC_DIR)/math/mymath.h \
	$(SRC_DIR)/math/mymath_complex.h \
	$(SRC_DIR)/math/mymath_internal.h \
	$(SRC_DIR)/math/standard_math_module.h \
	$(SRC_DIR)/matrix/matrix.h \
	$(SRC_DIR)/matrix/matrix_internal.h \
	$(SRC_DIR)/matrix/calculator_matrix_commands.h \
	$(SRC_DIR)/analysis/function_analysis.h \
	$(SRC_DIR)/analysis/multivariable_integrator.h \
	$(SRC_DIR)/analysis/multidim_integration.h \
	$(SRC_DIR)/analysis/ode_solver.h \
	$(SRC_DIR)/analysis/calculator_ode.h \
	$(SRC_DIR)/analysis/calculator_integration.h \
	$(SRC_DIR)/analysis/calculator_analysis_cmds.h \
	$(SRC_DIR)/analysis/calculator_rootfinding.h \
	$(SRC_DIR)/analysis/calculator_optimization.h \
	$(SRC_DIR)/analysis/calculator_series.h \
	$(SRC_DIR)/analysis/optimization_helpers.h \
	$(SRC_DIR)/symbolic/symbolic_expression.h \
	$(SRC_DIR)/symbolic/symbolic_expression_internal.h \
	$(SRC_DIR)/symbolic/calculator_transforms.h \
	$(SRC_DIR)/symbolic/calculator_symbolic_commands.h \
	$(SRC_DIR)/polynomial/calculator_polynomial.h \
	$(SRC_DIR)/polynomial/polynomial.h \
	$(SRC_DIR)/script/script_parser.h \
	$(SRC_DIR)/script/script_ast.h \
	$(SRC_DIR)/dsp/dsp.h \
	$(SRC_DIR)/dsp/calculator_signal_commands.h \
	$(SRC_DIR)/plot/plot_renderer.h \
	$(SRC_DIR)/plot/calculator_plot.h \
	$(SRC_DIR)/plot/plot_styles.h \
	$(SRC_DIR)/plot/svg_renderer.h
COMMON_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(COMMON_SRCS))
MAIN_OBJ := $(BUILD_DIR)/$(MAIN_SRC:.cpp=.o)
TEST_SRCS := $(TEST_DIR)/main.cpp $(wildcard $(TEST_SUITE_DIR)/*.cpp)
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))
DEPS := $(MAIN_OBJ:.o=.d) $(TEST_OBJS:.o=.d) $(COMMON_OBJS:.o=.d)

.PHONY: all test script-test check debug asan ubsan clean

all: $(APP)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BUILD_DIR)/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(APP): $(MAIN_OBJ) $(COMMON_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $(APP)

$(TEST_APP): $(TEST_OBJS) $(COMMON_OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $(TEST_APP)

test: $(TEST_APP)
	$(TEST_APP)

script-test: $(APP)
	test/script/run_symbolic_cli_validation.sh

check: test script-test

debug:
	$(MAKE) OPT_CXXFLAGS="-O0 -g"

asan:
	$(MAKE) OPT_CXXFLAGS="-O1 -g -fsanitize=address,undefined" LDFLAGS="-fsanitize=address,undefined"

ubsan:
	$(MAKE) OPT_CXXFLAGS="-O1 -g -fsanitize=undefined" LDFLAGS="-fsanitize=undefined"

clean:
	rm -rf $(BUILD_DIR) $(APP) $(TEST_APP)

-include $(DEPS)
