GCC_OPTIONS := -O3 -I./MQSMSim -std=c++20 -Wall -Wextra -pedantic
DEP_OPTIONS := -MMD
#DEP_OPTIONS := -MMD -MP
COMPILER := g++

BIN_PATH := ./bin
OBJ_PATH := ./obj
SRC_PATH := ./src

.PHONY: test clean

examples: $(SRC_PATH)/examples.cpp $(BIN_PATH)
	$(COMPILER) -o$(BIN_PATH)/examples $(SRC_PATH)/examples.cpp $(GCC_OPTIONS)


run_examples: $(BIN_PATH)/examples
	@echo "Running examples (RUN_TIMEOUT=$(RUN_TIMEOUT)s)..."
	@if [ "$(RUN_TIMEOUT)" -gt 0 ]; then \
		timeout "$(RUN_TIMEOUT)s" $(BIN_PATH)/examples; \
		status=$$?; \
		if [ $$status -eq 124 ]; then \
			echo "Examples stopped after $(RUN_TIMEOUT)s timeout."; \
		elif [ $$status -ne 0 ]; then \
			exit $$status; \
		fi; \
	else \
		$(BIN_PATH)/examples; \
	fi

$(BIN_PATH):
	mkdir -p $(BIN_PATH)

$(OBJ_PATH):
	mkdir -p $(OBJ_PATH)




testall: $(SRC_PATH)/testall.cpp $(BIN_PATH)
	$(COMPILER) -o$(BIN_PATH)/testall $(SRC_PATH)/testall.cpp $(GCC_OPTIONS)
	$(BIN_PATH)/testall



clean:
	rm -rf $(BIN_PATH) $(OBJ_PATH) 