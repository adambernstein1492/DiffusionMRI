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

				# If parameter value is just a number
				if (attribute[1][0] != '('):
					header[attribute[0]] = attribute[1][:-1]

				# If parameter value is a string, vector, or array
				if (attribute[1][0] == '('):
					attribute[1] = attribute[1].replace("(", "")
					attribute[1] = attribute[1].replace(")", "")
					attribute[1] = attribute[1].replace(" ", "")

					field_size = attribute[1][:-1].split(",")
					for i in range(len(field_size)):
						field_size[i] = int(field_size[i])

					field_size = tuple(field_size)
					header[attribute[0]] = []


					# Save spot in file
					prev_spot = f.tell()
					next_char = f.read(1)
					info = ''
					
					# Read file until next attribute is reached
					while (next_char != '#') and (next_char != '$'): 
						info += next_char

						next_char = f.read(1)
						if (next_char != '#') and (next_char != '$'):
							prev_spot = f.tell()
					
					# Split Info Field based on type of information
					# Handle Strings
					if info[0] == '<':
						info = info.replace('\n','').split('> <')

						for i in range(len(info)):
							info[i] = info[i].replace('<', '').replace('>', '')

					# Handle Cells
					elif info[0] == '(':
						info = info.replace('\n', '').split(') (')

						for i in range(len(info)):
							info[i] = info[i].replace('(', '').replace(')', '')

					elif (info[0] != '<') and (info[0] != '('):
						info = info.replace('\n', '').split()
						info = np.reshape(info, field_size)

						# Get rid of duplicate fields
						if len(info) > 10:
							is_repeated = True

							for i in range(len(info)):
								if (info[i] != info[0]):
									is_repeated = False

							if is_repeated:
								info = info[0]

					header[attribute[0]] = info
					# Rewind to beginning of next attribute and/or comment
					f.seek(prev_spot)

			line = f.readline()

	f.close()

	return header
