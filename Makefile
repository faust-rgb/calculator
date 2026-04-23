CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -O2

APP := calculator
TEST_APP := calculator_tests
SRC_DIR := src
TEST_DIR := test
SRC_DIRS := $(SRC_DIR)/app $(SRC_DIR)/core $(SRC_DIR)/math $(SRC_DIR)/matrix $(SRC_DIR)/analysis $(SRC_DIR)/algebra $(SRC_DIR)/symbolic $(SRC_DIR)/script
INCLUDES := $(addprefix -I,$(SRC_DIRS))

MAIN_SRC := $(SRC_DIR)/app/main.cpp
COMMON_SRCS := $(SRC_DIR)/core/calculator.cpp \
	$(SRC_DIR)/core/calculator_help.cpp \
	$(SRC_DIR)/math/mymath.cpp \
	$(SRC_DIR)/matrix/matrix.cpp \
	$(SRC_DIR)/analysis/function_analysis.cpp \
	$(SRC_DIR)/analysis/multivariable_integrator.cpp \
	$(SRC_DIR)/analysis/ode_solver.cpp \
	$(SRC_DIR)/symbolic/symbolic_expression.cpp \
	$(SRC_DIR)/algebra/polynomial.cpp \
	$(SRC_DIR)/script/script_parser.cpp
COMMON_HDRS := $(SRC_DIR)/core/calculator.h \
	$(SRC_DIR)/math/mymath.h \
	$(SRC_DIR)/matrix/matrix.h \
	$(SRC_DIR)/analysis/function_analysis.h \
	$(SRC_DIR)/analysis/multivariable_integrator.h \
	$(SRC_DIR)/analysis/ode_solver.h \
	$(SRC_DIR)/symbolic/symbolic_expression.h \
	$(SRC_DIR)/algebra/polynomial.h \
	$(SRC_DIR)/script/script_parser.h \
	$(SRC_DIR)/script/script_ast.h

.PHONY: all test clean

all: $(APP)

$(APP): $(MAIN_SRC) $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(MAIN_SRC) $(COMMON_SRCS) -o $(APP)

$(TEST_APP): $(TEST_DIR)/tests.cpp $(COMMON_SRCS) $(COMMON_HDRS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEST_DIR)/tests.cpp $(COMMON_SRCS) -o $(TEST_APP)

test: $(TEST_APP)
	./$(TEST_APP)

clean:
	rm -f $(APP) $(TEST_APP)
