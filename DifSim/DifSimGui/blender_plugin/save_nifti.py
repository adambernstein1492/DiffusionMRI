import sys
import numpy
from nipy.core.image.image import Image, AffineTransform, CoordinateSystem
from nipy import load_image, save_image

def save_nifti(signal_file, output_file):
    S = (numpy.genfromtxt(signal_file)[:,0])
    ndir = S.shape[-1]
    # No longer need to stack b0 value, included by default
    #S = numpy.hstack(([1],S))
    S = S.reshape([1,1,1,ndir])
    datdims = len(S.shape)
    aff = numpy.zeros((datdims+1,datdims+1))
    aff[0,0] = 1
    aff[1,1] = 1
    aff[2,2] = 1
    aff[0,datdims] = S.shape[0]
    aff[1,datdims] = S.shape[1]
    aff[2,datdims] = S.shape[2]
    aff[3,3] = 1
    aff[3,datdims] = S.shape[3]
    if (datdims > 3):
      aff[4,datdims] = 1
    
    at = AffineTransform(function_domain=CoordinateSystem(coord_names=('i', 'j', 'k', 'l'), name='', coord_dtype=numpy.float64),function_range=CoordinateSystem(coord_names=('x', 'y', 'z', 't'), name='', coord_dtype=numpy.float64),affine=aff)
    nim = Image(S, at)
    save_image(nim,output_file)
    
if __name__ == "__main__":
    signal_file, output_file = sys.argv[1:]
    save_nifti(signal_file, output_file)
