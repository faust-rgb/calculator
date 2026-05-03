CXX ?= g++
BASE_CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic
OPT_CXXFLAGS ?= -O2 -static
CXXFLAGS ?= $(BASE_CXXFLAGS) $(OPT_CXXFLAGS)
LDFLAGS ?=
BIN_DIR := bin
BUILD_DIR := build
APP := $(BIN_DIR)/calculator
TEST_APP := $(BIN_DIR)/calculator_tests
SRC_DIR := src
TEST_DIR := test
TEST_SUITE_DIR := $(TEST_DIR)/suites

# Collect all source directories using wildcard
SRC_DIRS := $(shell find $(SRC_DIR) -type d)
INCLUDES := -I$(SRC_DIR) $(addprefix -I,$(SRC_DIRS)) -I$(TEST_DIR)
CPPFLAGS += $(INCLUDES) -MMD -MP

# Collect all source files using wildcard (exclude main.cpp)
MAIN_SRC := $(SRC_DIR)/app/main.cpp
COMMON_SRCS := $(shell find $(SRC_DIR) -name "*.cpp" -not -path "$(SRC_DIR)/app/main.cpp")
COMMON_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(COMMON_SRCS))
MAIN_OBJ := $(BUILD_DIR)/$(MAIN_SRC:.cpp=.o)

# Test sources
TEST_SRCS := $(TEST_DIR)/main.cpp $(wildcard $(TEST_SUITE_DIR)/*.cpp)
TEST_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SRCS))

# Dependencies
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