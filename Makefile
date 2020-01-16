all: MinimizeGrid.cpp ReconstituteGrid.cpp MinimizeGrid ReconstituteGrid

clean:
	/bin/rm -f MinimizeGrid ReconstituteGrid

.cpp:
	g++ -O2 $< -o $@
