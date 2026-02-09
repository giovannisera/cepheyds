# Compiler
FC = gfortran

# Flags 
FFLAGS = -O3 -Wall -Wextra 

# Executable file
TARGET = analysis_cepheyds

# Source file
SRC = 0_main.f90

## RULES
all: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET) *.mod *.o
	@echo "The cleaning has been completed. Files .mod and .o removed."

# Help
help:
	@echo "Available commands"
	@echo "  make        - Compiles the programm $(TARGET)"
	@echo "  make clean  - Removes the binary files, the modules and the objects"
