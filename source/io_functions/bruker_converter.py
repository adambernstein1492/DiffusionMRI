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
							info[i] = info[i].replace('(', '').replace(')', '').replace('<','').replace('>','')
							info[i] = info[i].split(', ')

					# Handle Arrays and Vectors
					elif (info[0] != '<') and (info[0] != '('):
						info = info.replace('\n', '').split()
						info = np.reshape(info, field_size)

						# Get rid of duplicate fields
						if len(info) > 10:
							is_repeated = True

							for i in range(len(info)):
								if (np.any(info[i] != info[0])):
									is_repeated = False

							if is_repeated:
								info = info[0]

					header[attribute[0]] = info
					# Rewind to beginning of next attribute and/or comment
					f.seek(prev_spot)

			line = f.readline()

	f.close()

	return header

def get_affine(visu_pars):

	affine_orientation = visu_pars['VisuCoreOrientation']
	affine_orientation = np.reshape(affine_orientation, (3,3))
	for i in range(3):
		for j in range(3):
			affine_orientation[i,j] = float(affine_orientation[i,j])

	resolution = get_resolution(visu_pars)
	offset = get_offset(visu_pars)

	affine = np.eye(4, dtype=np.float32)
	affine[0:3,0:3] = affine_orientation
	affine[0:3,3] = offset

	affine = np.linalg.inv(affine)

	affine[0:3,0:3] = affine[0:3,0:3].dot(np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]]))

	for i in range(3):
		max = 0
		index = 0
		for j in range(3):
			if np.abs(affine[j,i]) > max:
				max = np.abs(affine[j,i])
				index = j

		if i != 2:
			if affine[index,i] > 0:
				affine[:,i] *= -1
		if i == 2:
			if affine[index,i] < 0:
				affine[:,i] *= -1


	affine[0:3,0:3] = affine[0:3,0:3].dot(np.diag(resolution))

	return affine

def get_offset(visu_pars):

	offset = np.zeros(3)
	offset_values = np.zeros(visu_pars['VisuCorePosition'].shape)

	for i in range(visu_pars['VisuCorePosition'].shape[0]):

		offset_values[i][0] = float(visu_pars['VisuCorePosition'][i][0])
		offset_values[i][1] = float(visu_pars['VisuCorePosition'][i][1])
		offset_values[i][2] = float(visu_pars['VisuCorePosition'][i][2])

	x = np.amin(offset_values,axis=0)

	offset[0] = x[0]
	offset[1] = x[1]
	offset[2] = x[2]

	return offset

def get_resolution(visu_pars):

	resolution = np.zeros(3)

	# 2D case
	if visu_pars['VisuCoreDim'] == '2':
		resolution[0] = float(visu_pars['VisuCoreExtent'][0]) / float(visu_pars['VisuCoreSize'][0])
		resolution[1] = float(visu_pars['VisuCoreExtent'][1]) / float(visu_pars['VisuCoreSize'][1])

		depth = np.zeros(3)
		for i in range(len(depth)):
			depth[i] = (float(visu_pars['VisuCorePosition'][-1][i]) - float(visu_pars['VisuCorePosition'][0][i]))

		depth = np.linalg.norm(depth)

		for i in range(len(visu_pars['VisuFGOrderDesc'])):
			if visu_pars['VisuFGOrderDesc'][i][1] == 'FG_SLICE':
				num_slices = float(visu_pars['VisuFGOrderDesc'][i][0])

		resolution[2] = depth / (num_slices - 1)

	# 3D case
	if visu_pars['VisuCoreDim'] == '3':
		resolution[0] = float(visu_pars['VisuCoreExtent'][0]) / float(visu_pars['VisuCoreSize'][0])
		resolution[1] = float(visu_pars['VisuCoreExtent'][1]) / float(visu_pars['VisuCoreSize'][1])
		resolution[2] = float(visu_pars['VisuCoreExtent'][2]) / float(visu_pars['VisuCoreSize'][2])

	return resolution

def get_matrix(visu_pars):

	matrix = np.zeros(3)

	# 2D case
	if visu_pars['VisuCoreDim'] == '2':
		if visu_pars['VisuFGOrderDescDim'] == '2':
			matrix = np.zeros(4)

		matrix[0] = int(visu_pars['VisuCoreSize'][0])
		matrix[1] = int(visu_pars['VisuCoreSize'][0])

		for i in range(len(visu_pars['VisuFGOrderDesc'])):
			if visu_pars['VisuFGOrderDesc'][i][1] == 'FG_SLICE':
				matrix[2] = int(visu_pars['VisuFGOrderDesc'][i][0])

			if visu_pars['VisuFGOrderDesc'][i][1] == 'FG_MOVIE':
				matrix[3] = int(visu_pars['VisuFGOrderDesc'][i][0])

	# 3D case
	if visu_pars['VisuCoreDim'] == '3':
		matrix[0] = int(visu_pars['VisuCoreSize'][0])
		matrix[1] = int(visu_pars['VisuCoreSize'][1])
		matrix[2] = int(visu_pars['VisuCoreSize'][2])

	return tuple([int(i) for i in matrix])

def convert_to_nifti(recon_data, visu_pars, out_path):

	# Get Word Size
	word_size = visu_pars['VisuCoreWordType']

	if word_size == '_16BIT_SGN_INT':
		word_size = np.int16

	# Get Matrix Size
	matrix = get_matrix(visu_pars)

	# Get affine transformation
	affine = get_affine(visu_pars)

	img_data = np.fromfile(recon_data, dtype=word_size)
	img_data = np.reshape(img_data, matrix, order='F')

	# Convert to Nifti and Save
	img = nib.Nifti2Image(img_data, affine)
	nib.save(img, out_path)
