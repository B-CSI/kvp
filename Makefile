# Compiler
CC = gcc

# Flags
FLAGS = -O3 -march=native -shared -fPIC -fopenmp
FLAGS_WHL = -O3 -shared -fPIC -fopenmp

# Files
FILES = kvp/src/*.c

# Lib
LIB = kvp/lib/libkvp-
UNIX = unix.so
MAC = maco.so


# Linux for same machine
all : kvp

# MacOS for same machine
mac : kvp_mac

# Linux for wheel
linuxwhl : kvp_linux_wheel

# MacOS for wheel
macoswhl : kvp_macos_wheel

kvp : $(FILES)
	$(CC) -o $(LIB)$(UNIX) $(FILES) $(FLAGS)

kvp_mac : $(FILES)
	$(CC) -o $(LIB)$(MAC) $(FILES) $(FLAGS)

kvp_linux_wheel : $(FILES)
	$(CC) -o $(LIB)$(UNIX) $(FILES) $(FLAGS_WHL)

kvp_macos_wheel : $(FILES)
	$(CC) -o $(LIB)$(MAC) $(FILES) $(FLAGS_WHL)
