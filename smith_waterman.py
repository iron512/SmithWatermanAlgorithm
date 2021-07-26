#!/usr/bin/python

#TGTTACGG GGTTGACTA

import argparse

#auxiliary function for applying the colors
class bcolors:
	GREEN = '\033[92m'
	YELLOW = '\033[93m'
	CYAN = '\033[96m'
	RED = '\033[91m'
	RESET = '\033[0m'

def grn(text):
	return bcolors.GREEN + str(text) + bcolors.RESET

def ylw(text):
	return bcolors.YELLOW + str(text) + bcolors.RESET

def red(text):
	return bcolors.RED + str(text) + bcolors.RESET

#the custom match function. Can be updated to provide differents rules other than match and mistmatch.
def cmatch(a, b, positive, negative):
	if a == b:
		return positive
	else:
		return negative

#the custom print function. The first row/columns are substituted by the actual letters of the sequences.
#online versions also prompt the first rows/columns (filled with zeros) but I feel that this way of displaying the 
#result is more intuitive and understandable. 
def pprint(s1, s2, matrix):
	for i in range(0, len(s1)+1):
		for j in range(0, len(s2)+1):
			if i == 0 and j == 0:
				print("\\", end="\t")
			elif i == 0:
				print(s2[j-1], end="\t")
			elif j == 0:
				print(s1[i-1], end="\t")
			else:
				print(matrix[i][j], end="\t")		
		print()

#A little modification of the function presented above. The headers are colored in green, but also the actual sequence is highlighted in yellow.
#Since it was quite hard to display the arrows in the terminal, I adopted this solution.
#It works using the terminal colors. Not suited for storage in txt files (-m std fits it better.)
def pprint_color(s1, s2, matrix, color):
	for i in range(0, len(s1)+1):
		for j in range(0, len(s2)+1):
			if i == 0 and j == 0:
				print(grn("\\"), end="\t")
			elif i == 0:
				print(grn(str(s2[j-1])), end="\t")
			elif j == 0:
				print(grn(str(s1[i-1])), end="\t")
			else:
				if (i,j) not in color:
					print(matrix[i][j], end="\t")
				else:
					print(ylw(str(matrix[i][j])), end="\t")	
		print()

def main():
	#setup the parser
	parser = argparse.ArgumentParser(description = "A basic implementation of the smith-waterman algorithm for the local alignment.\nThis code has been developed as part of the final exam of the course of Algorithms for Bioinformatics 2020/2021, University of Trento",formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	parser.add_argument("s1", type=str, help="the first sequence, coded as a String")
	parser.add_argument("s2", type=str, help="the second sequence, coded as a String")

	parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit\n\n')

	parser.add_argument("-sm","--scoring_match",type = int, help="The match value for the custom storing system (default = 3)\n", default= 3)
	parser.add_argument("-sx","--scoring_mismatch",type = int, help="The mismatch value for the custom storing system (default = -3)\n", default = -3)
	parser.add_argument("-sg","--scoring_gap",type = int, help="The gap value for the custom storing system (default = -2)\n\n", default = -2)

	parser.add_argument("-m","--matrix", help="Show the traceback matrix.\n", default='none', choices=['none', 'std', 'color'])

	#parse the user's input
	args = parser.parse_args()

	#transform the input in an uppercase array of chars
	seq1 = [char for char in args.s1.upper()]
	seq2 = [char for char in args.s2.upper()]

	match = args.scoring_match
	mismatch = args.scoring_mismatch
	gap = args.scoring_gap 

	output = args.matrix

	#initialize through numpy the required matrices
	matrix = []
	direction = []
	for x in range(0,len(seq1)+1):
		matrix.append([])
		direction.append([])
		for y in range(0,len(seq2)+1):
			matrix[x].append(0)
			direction[x].append([])


	#store the maximal values in just one pass
	current_max = -1
	current_max_pos = []

	for i in range(1, len(seq1)+1):
		for j in range(1, len(seq2)+1):			
			#generate the candidates
			diag = matrix[i-1][j-1] + cmatch(seq1[i-1], seq2[j-1], match, mismatch)  
			up = matrix[i-1][j] + gap
			left = matrix[i][j-1] + gap 
			new = 0

			#store the highest candidate(s)
			matrix[i][j] = max(diag, up, left, new)
			
			if matrix[i][j] > current_max:
				current_max = matrix[i][j]
				current_max_pos = [[i,j]]
			#Store all the highest values. Just ONE path each is explored
			elif matrix[i][j] == current_max:
				current_max_pos.append([i,j])
			
			#keep track of the candidates, needed for backtracking the algorithm.
			if (matrix[i][j] == diag):
				direction[i][j].append('d')
			if (matrix[i][j] == up):
				direction[i][j].append('u')
			if (matrix[i][j] == left):
				direction[i][j].append('l')
			if (matrix[i][j] == new):
				direction[i][j].append('n')

	print("SEQ1:","".join(seq1))
	print("SEQ2:","".join(seq2))
	print("SCORE:",current_max)
	print()

	#the pair will always be one (The last, arbitrary choice), but the algorithm is ready to explore multiple paths.
	for pair in current_max_pos[-1:]:
		x = pair[0]
		y = pair[1]

		#store the backtracked values and their position (to color the path)
		alignA = []
		matches = []
		alignB = []
		color = []

		#start from the end, go until a cell found has been inizialized by a "new" (or the boundaries are exceeded).
		while x>0 and y>0 and (not direction[x][y][0] == 'n'):
			#parse the match/mismatch
			color.append((x,y))
			if direction[x][y][0] == 'd':		
				alignA.append(seq1[x-1])
				alignB.append(seq2[y-1])
				if (seq1[x-1] == seq2[y-1]):
					matches.append("|")
				else:
					matches.append("X")
				x-=1
				y-=1

			#parse the gaps
			elif direction[x][y][0] == 'u':	
				alignA.append('-')
				alignB.append(seq1[x-1])
				matches.append(" ")
				x-=1

			elif direction[x][y][0] == 'l':	
				alignA.append(seq2[y-1])
				alignB.append('-')
				matches.append(" ")
				y-=1
		color.append((x,y))

		#reverse the backtracked values, in order to produce the correct string and print them as a string
		alignA.reverse()
		matches.reverse()
		alignB.reverse()

		print("".join(alignA))
		print("".join(matches))
		print("".join(alignB))
		print()

		if output == "color":
			pprint_color(seq1,seq2,matrix,color)

	if output == "std":
		pprint(seq1,seq2,matrix)

#MAIN
if __name__ == '__main__':
	main()