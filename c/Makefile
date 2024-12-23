# Compilator 
CC = gcc 

# Files 
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Source files 
SRCS = $(SRC_DIR)/functions.c $(SRC_DIR)/mainForVectors.c

# Object files 
OBJS = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRCS))

# Executable name 
TARGET = $(BIN_DIR)/programme 

# Create folders if necessary 
MKDIR_P = mkdir -p 

# Default rule
all: $(TARGET)

# Executable rule  
$(TARGET): $(OBJS)
	$(MKDIR_P) $(BIN_DIR)
	$(CC) $(OBJS) -o $(TARGET) -lm

# Rule to create objects files in obj folder
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(MKDIR_P) $(OBJ_DIR)
	$(CC) -c $< -o $@

# Clean generate files
clean: 
	rm -rf $(OBJ_DIR) $(BIN_DIR)

# Clean only objects files 
clean-obj: 
	rm -rf $(OBJ_DIR)

# Rule to run gnuplot for plotting
gnuplot:
	gnuplot -p -e "plot 'data.txt' with lines"


.PHONY: all clean clean-obj gnuplot 