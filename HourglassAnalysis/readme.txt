To compile:
javac HourglassAnalysis.java

To run:
java HourglassAnalysis

To change input specification:
Go to  main function: 
	change the input file name, source file name and target file name
	recompile and run
A sample input file, source file and target file is provided
To change the coverageThreshold(tau), change the variable pathCoverageTau value (default set at 90%)

Input file format:
Edge list, one edge per line in the following format: from_node <spaces> to_node

Source file format:
One source node identifier per line

Target file format:
One target node identifier per line

Output:
The code print some network statistics in the console.
For each node it prints generality, complexity and path centrality.
It also prints the total number of paths, one sample core and core size.
The above is printed for both original and flattened network.
It also prints the hScore.
Note that the code creates a temp file named "flat.txt" which is the flattened network.
