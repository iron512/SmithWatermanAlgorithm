# Smith Waterman Algorithm
This repository contains all the material regarding the final project for the course of Algorithms for bioinformatics, for the year 2020/2021.
The task is an implementation of the Smith Waterman Algorithm (provided here in python)

Ivan Martini #2075971

### Setup
0) **Prerequisites**: Have `python` installed. Nothing more should be required on any operating system. I had the chance to test the code only on Archlinux and Ubuntu 20.04, but everything went fine.

1) Clone the repository
```
git clone https://github.com/iron512/SmithWatermanAlgorithm.git
```

2) Run the code
```
python3 smith_waterman.py TGTTGTAG TGTTGAT
```

*Note: the script has an helper. using the command -h it will show all the parameters that can be added to tune the algorithm*

3) Examples

Run the code changing the default scoring system (+3,-3,-2)
```
python3 smith_waterman.py TGTTACGG GGTTGACTA -sm +1 -sx -1 -sg -2
```

Run the code showing in the end the tracing matrix
```
python3 smith_waterman.py TGTTACGG GGTTGACTA -m std
```

Run the code showing in the end the tracing matrix, highlighting also the path in yellow. suitable for a terminal display, but not for storing on a txt file (use -m std instead)
```
python3 smith_waterman.py TGTTACGG GGTTGACTA -m color
```