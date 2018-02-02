#!/usr/bin/env python2.7
"""
**** Plot results of Bloch simulation bloch.c ****
"""
print __doc__

#-------- Libraries --------------

import sys
sys.path.insert(0, '/Users/lrf/bin')
from lrf_utils import *

#-------- Functions --------------

def plot_sphere(fig,n=100,rad=1):
    u = np.linspace(0, 2*np.pi, n) 
    v = np.linspace(0, np.pi, n) 
    x = rad*np.outer(np.cos(u), np.sin(v)) 
    y = rad*np.outer(np.sin(u), np.sin(v)) 
    z = rad*np.outer(np.ones(np.size(u)), np.cos(v)) 
    ax = fig.gca(projection='3d') 
    ax.plot_surface(x, y, z, rstride=2, cstride=2, 
                    linewidth=.1, color='y', alpha=0.5) 
    minbound = -1.
    maxbound =  1.
    ax.auto_scale_xyz([minbound, maxbound], 
                      [minbound, maxbound], 
                      [minbound, maxbound])
    #ax.grid(False)
    for a in ax.w_xaxis.get_ticklabels():
        a.set_visible(False)
    for a in ax.w_yaxis.get_ticklabels():
        a.set_visible(False)
    for a in ax.w_zaxis.get_ticklabels():
        a.set_visible(False)

#-------- Main --------------

def main():

    if 0:
        rf2=np.fromfile('slr180.rho',dtype=np.int16)
        rho_slr=np.fromfile('msrf2_slr.rho',dtype=np.int32)
        theta_slr=np.fromfile('msrf2_slr.theta',dtype=np.int32)
        rho_sinc=np.fromfile('msrf2_sinc.rho',dtype=np.int32)
        theta_sinc=np.fromfile('msrf2_sinc.theta',dtype=np.int32)

        fig_sinc=plt.figure(figsize=(15, 5))
        add_quit_button()

        ax1 = plt.subplot(1,2,1)
        plt.plot(rho_sinc)
        plt.title('sinc rho')
        ax1.grid(ls='--')
        ax1.autoscale(True,axis='x',tight=True);

        ax2 = plt.subplot(1,2,2)
        plt.plot(theta_sinc)
        plt.title('sinc theta')
        ax2.grid(ls='--')
        ax2.autoscale(True,axis='x',tight=True);

        fig_slr = plt.figure(figsize=(15, 5))
        add_quit_button()

        ax1 = plt.subplot(1,2,1)
        plt.plot(rho_slr)
        plt.title('slr rho')
        ax1.grid(ls='--')
        ax1.autoscale(True,axis='x',tight=True);

        ax2 = plt.subplot(1,2,2)
        plt.plot(theta_slr)
        plt.title('slr theta')
        ax2.grid(ls='--')
        ax2.autoscale(True,axis='x',tight=True);
    
    # Bloch sim results

    bloch_rho   = np.fromfile('bloch.rho',dtype=np.double)
    bloch_theta = np.fromfile('bloch.theta',dtype=np.double)
    print 'max(bloch_rho) = ',np.max(bloch_rho)
    print 'max(bloch_theta) = ',np.max(bloch_theta)

    f=open('Mxyz.dat')
    d=f.readlines()
    npts = len(d)
    print 'npts = ',npts
    x=map(lambda s: np.array(s.strip().split(' ')).astype(float),d)
    Mxyz=np.array(x)
    print 'min(Mxyz) = ',np.min(Mxyz), ', max(Mxyz) = ',np.max(Mxyz)
    
    #----- Plot Mxyz components separately ----------#

    fig_mxyz=plt.figure(figsize=(15, 5))
    add_quit_button()

    ax1 = plt.subplot(1,3,1)
    plt.plot(bloch_rho)
    plt.title('$A(t)$')
    ax1.grid(ls='--')
    ax1.autoscale(True,axis='x',tight=True);
    ax2 = plt.subplot(1,3,2)
    plt.plot(bloch_theta)
    plt.title('$\phi(t)$')
    ax2.grid(ls='--')
    ax2.autoscale(True,axis='x',tight=True);
    ax3 = plt.subplot(1,3,3)
    plt.plot(Mxyz)
    plt.title('$M_{xyz}$')
    ax3.grid(ls='--')
    ax3.autoscale(True,axis='x',tight=True);

    #----- Plot Mxyz on Bloch Sphere ----------#

    fig_bloch=plt.figure(figsize=(5, 5))
    add_quit_button()
    ax_bloch = fig_bloch.gca(projection='3d')
    Mxyz = np.swapaxes(Mxyz,0,1)
    ax_bloch.plot(Mxyz[0], Mxyz[1], Mxyz[2], label='parametric curve')
    plt.title('$M_{xyz}$')
    plot_sphere(fig_bloch)

    plt.show()

if __name__ == "__main__":
    main()
