import numpy as np
import nibabel as nib

def read_bruker_header_file(filename):
	
	f = open(filename, 'r')
	
	line = f.readline()
	if (line[0:7] == "##TITLE"):
		header = {}

		while line:
			if (line[0:3] == "##$"):
				attribute = line[3:].split('=')
				header[attribute[0]] = []
			line = f.readline()

	f.close()

	return header
