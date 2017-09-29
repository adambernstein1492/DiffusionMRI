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

				if (attribute[1][0] != '('):
					header[attribute[0]] = attribute[1][:-1]

				if (attribute[1][0] == '('):
					attribute[1] = attribute[1].replace("(", "")
					attribute[1] = attribute[1].replace(")", "")
					attribute[1] = attribute[1].replace(" ", "")

					field_size = attribute[1][:-1].split(",")
					for i in range(len(field_size)):
						field_size[i] = int(field_size[i])

					field_size = tuple(field_size)
					header[attribute[0]] = np.zeros(field_size)

			line = f.readline()

	f.close()

	return header
