import os
import sys
from sys import argv

def main():
	script, filename = argv

	file = open(filename, 'r')
	counter = 0
	correct = 0

	for line in file:
		if "0000" in line:
			counter += 1
		if(counter < 513):
			if "d0 22" in line:
				correct += 1

	file.close()
	print(counter)
	print(correct)

main()

