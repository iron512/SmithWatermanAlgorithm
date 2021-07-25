#!/usr/bin/python

import argparse
import numpy

def cmatch(a, b):
	if a == b:
		return +3
	else:
		return -3

def pprint(s1,s2,matrix):
	for i in range(0, len(s1)+1):
		for j in range(0, len(s2)+1):
			if i == 0 and j == 0:
				print(0, end="\t")
			elif i == 0:
				print(s2[j-1], end="\t")
			elif j == 0:
				print(s1[i-1], end="\t")
			else:
				print(matrix[i][j], end="\t")		
		print()


def main():
	#SETUP THE PARSER
	parser = argparse.ArgumentParser(description = "A basic implementation of the smith-waterman algorithm for the local alignment.\nThis code has been developed as part of the final exam of the course of Algorithms for Bioinformatics 2020/2021, University of Trento",formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	parser.add_argument("s1", type=str, help="the first sequence, coded as a String")
	parser.add_argument("s2", type=str, help="the second sequence, coded as a String")

	parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit\n\n')

	parser.add_argument("-sm","--scoring_match",type = int, help="The match value for the custom storing system (default = 3)\n", default= 3)
	parser.add_argument("-sx","--scoring_mismatch",type = int, help="The mismatch value for the custom storing system (default = -3)\n", default = -3)
	parser.add_argument("-sg","--scoring_gap",type = int, help="The gap value for the custom storing system (default = -2)\n\n", default = -2)

	#parser.add_argument("-m","--matrix", help="a file CSV containing the PAM matrix wanted")

	#PARSE THE USER'S INPUT
	args = parser.parse_args()

	seq1 = [char for char in args.s1.upper()]
	seq2 = [char for char in args.s2.upper()]

	match = args.scoring_match
	mismatch = args.scoring_mismatch
	gap = args.scoring_gap 

	#TGTTACGG GGTTGACTA
	matrix = numpy.zeros((len(seq1) + 1, len(seq2) + 1), dtype=object)
	direction = numpy.zeros((len(seq1) + 1, len(seq2) + 1), dtype=object)
	current_max = -1
	current_max_pos = []

	for i in range(1, len(seq1)+1):
		for j in range(1, len(seq2)+1):
			direction[i][j] = []
			
			diag = matrix[i-1][j-1] + cmatch(seq1[i-1], seq2[j-1])  
			up = matrix[i-1][j] -2
			left = matrix[i][j-1] -2 
			new = 0

			matrix[i][j] = max(diag, up, left, new)
			if matrix[i][j] > current_max:
				current_max = matrix[i][j]
				current_max_pos = [i,j]

			if (matrix[i][j] == diag):
				direction[i][j].append('d')
			if (matrix[i][j] == up):
				direction[i][j].append('u')
			if (matrix[i][j] == left):
				direction[i][j].append('l')
			if (matrix[i][j] == new):
				direction[i][j].append('n')

	x = current_max_pos[0]
	y = current_max_pos[1]

	alignA = []
	alignB = []

	while matrix[x][y] != 0:
		if direction[x][y][0] == 'd':		
			alignA.append(seq1[x-1])
			alignB.append(seq2[y-1])
			x-=1
			y-=1

		elif direction[x][y][0] == 'u':	
			alignA.append('-')
			alignB.append(seq1[x-1])
			x-=1

		elif direction[x][y][0] == 'l':	
			alignA.append(seq2[y-1])
			alignB.append('-')
			y-=1

	pprint(seq1,seq2,matrix)

	alignA.reverse()
	alignB.reverse()
	print("".join(alignA))
	print("".join(alignB))
	print()



#MAIN
if __name__ == '__main__':
	main()