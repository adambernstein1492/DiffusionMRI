# lrf utils

import os
import sys
import gtk
import time as tm

import matplotlib
matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from pylab import *

import numpy as np
import numpy.polynomial.polynomial as P
import scipy as sp

from scipy.sparse import issparse
from scipy.sparse import coo_matrix
from itertools import *
import code

# Interaction

def keyboard(banner=None):
    ''' Function that mimics the matlab keyboard command '''
    # use exception trick to pick up the current frame
    try:
        raise None
    except:
        frame = sys.exc_info()[2].tb_frame.f_back
    print "# Use quit() to exit :) Happy debugging!"
    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner, local=namespace)
    except SystemExit:
        return 

# Graphics

# open the disk file 'fname' in for 'mode'
def open_file(fname,mode):
    if (mode=='w'):
        f = open(fname,mode)
        return f
    elif (mode=='r'):
        if (os.path.exists(fname)):
            f = open(fname,mode)
        else:
            exit_str = '**** %s not found!' %(fname) + '...aborting'
            sys.exit(exit_str)
        return f

# delete a file from disk
def delete_file(fname):
    delete_str = 'rm -f ' + fname
    os.system(delete_str)

def clearticklabels(ax):
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

def clearticks(ax):
    ax.set_xticks(())
    ax.set_yticks(())

def add_quit_button():
    manager = get_current_fig_manager()
    toolbar = manager.toolbar
    quit_button_loc = 8; # where to insert this in the mpl toolbar
    quit_button = gtk.Button("Quit")
    quit_button.connect('clicked',gtk.main_quit)
    quit_button.show()
    toolitem = gtk.ToolItem()
    toolitem.show()
    toolitem.add(quit_button)
    toolbar.insert(toolitem, quit_button_loc); 

# Structural

def vectorize(x):
    x=x.reshape(np.prod(np.shape(x)),1)
    return x

def unvectorize(x):
    x = np.ndarray.flatten(x)
    return x

def zeropad(x,npx,npy):
    return np.pad(x,(npx,npy),'constant',constant_values=(0,0))

def spacingsN(origin, size, delta):
    return [linspace(ox,ox+s*d,s) for ox,s,d in zip(origin,size,delta)]

def tosparse(x):
    idx = x.nonzero()
    row = idx[0]
    col = idx[1]
    dat = x[idx]
    z = coo_matrix((dat,(row,col)),shape=x.shape)
    return z

def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a

def create_false_list(n):
    l=[]
    for k in range(n): l.append(False)
    return l

# Combinatorics

# works like matlab's version
def nchoosek(n,k,m=0):
    return np.array(list(combinations(range(m,n),k)))

# Math

def nextpow2(n):
    nn = np.rint(exp(log(2)*ceil(log(n+1)/log(2))))
    return nn

# sparse matrix power                                                                                
def sparseMatrixPower(x,n):
    xpow = x
    for i in range(n-1):
        xpow = xpow.dot(x)
    return xpow

# find maximum of each column of array arr
def max_over_cols(arr):
    x = np.array(map(max,zip(*arr)))
    return x.reshape((len(x),1))

# find maximum of each row of array arr
def max_over_rows(arr):
    x = np.array(map(max,zip(*arr.T)))
    return x.reshape((len(x),1))

def evenQ(n):
    return (mod(n,2)==0)

def oddQ(n):
    return (mod(n,2)!=0)

def rms(data,axis=None):
    return np.sqrt(np.mean(data**2,axis))

# I/O

def printRange(x,xname,tab=False):
    if issparse(x)==True:
        if tab==True:
            print '\trange(%s) = (%.2g,%.2g)' % (xname, x.data.min(), x.data.max())
        else:
            print 'range(%s) = (%.2g,%.2g)' % (xname, x.data.min(), x.data.max())
    else:
        if tab==True:
            print '\trange(%s) = (%.2g,%.2g)' % (xname, x.min(), x.max())
        else:
            print 'range(%s) = (%.2g,%.2g)' % (xname, x.min(), x.max())

def printShape(x,xname,tab=False):
    if tab==True:
        print '\tshape(%s) = ' % (xname),x.shape
    else:
        print 'shape(%s) = ' % (xname),x.shape

def logqyn(qstr,qdef,tab=False):
    rawstr = qstr + '? (y/n): ['+ str(qdef)+']: '
    if tab==True:
        rawstr = '\t'+rawstr
    q = raw_input(rawstr)
    if not q: q=qdef

    if q=='y':
        return True
    elif q=='n':
        return False
    else:
        print "logqyn: unknown data type %s !" % q
        print "(returning False)"
        return False

def pinput(pstr,pdef,dtype='float',lims=None,tab=False):
    if lims!=None:
        rawstr = 'Enter ' + pstr + ' ('+ str(lims)+') ' + '['+ str(pdef)+']: '
    else:
        rawstr = 'Enter ' + pstr + ' ['+ str(pdef)+']: '

    if tab==True:
        rawstr = '\t'+rawstr

    p = raw_input(rawstr)
    if not p: p=pdef
    if dtype=='float':
        p = float(p)
    elif dtype=='int':
        p=int(float(p))
    elif dtype=='char':
        p=p

    if (lims!=None) :
        if dtype=='char':
            while len([y for y in lims if p in y])==0:
                print 'input',p,'not a valid option',lims,'- try again...'
                p = raw_input(rawstr)
                if not p: p=pdef
            
        else:
            while any([p<lims[0] , p>lims[1]]):
                print 'input',p,'not within limits',lims,'- try again...'
                p = raw_input(rawstr)
                if not p: p=pdef
                if dtype=='float':
                    p = float(p)
                elif dtype=='int':
                    p = int(p)

    if dtype=='float':
        return float(p)
    elif dtype=='int':
        return int(float(p))
    elif dtype=='char':
        return p
    else:
        print "pinput: unknown data type!"
        print "(returning float)"
        return float(p)

