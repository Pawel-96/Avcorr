#! /usr/bin/python3

#converting ASCII lightcone into 3D HDF5 file: [ra,dec,z] -> [X,Y,Z], imposing angular and radial cuts
#also creating randoms (positions of spheres/ellipsoid centers for CiC)
#redshift converted into comoving distance
#Attention: randoms assume constant catalog radial density/completeness

import numpy as np
import math
from astropy.cosmology import LambdaCDM
import h5py
import os.path

deg2rad=np.pi/180. #for converting degrees->radians
pihalf=.5*np.pi

#****************parameters***************

fname='../Angular/Data/Lightcone_P18_CENTRAL_BOX1_60x60_z_0p2_0p3.txt' #lightcone file to cut
ramin,ramax=[0.,60.] #right-ascension lower and upper cut
decmin,decmax=[-30.,30.] #declination lower and upper cut
zmin,zmax=[0.2,0.3] #redshift lower and upper cut

col_ra,col_dec,col_z=[0,1,4] #columns in lightcone file with right-scaension,declination and redshift
#using col_z corresponding with observed redshift will create catalog with RSD
maxR=100. #maximal size of kernel in CiC method - randoms will be not closer than maxR to catalog boundaries

#cosmology:
Om_M=0.3111
Om_L=1.-Om_M
h=0.6766 #H0/(100 km/s/Mpc)

#output:
fname_out='Data/Lightcone_P18_CENTRAL_BOX1_3Dcut.hdf5'
fname_random='Randoms/Random_Lightcone_P18_CENTRAL_BOX1_3Dcut.hdf5'
nmultiply=2 #multiplication factor: how many times randoms will be denser than data (approximately)
dpos='coordinates' #dataset name for coordinates

#*****************************************

cosmology=LambdaCDM(H0=100.*h, Om0=Om_M, Ode0=Om_L)


def Make_cuts(ra,dec,zz): #cutting and converting to cartesian coords
    print('Making cuts')
    nobj=len(ra) #number of objects
    indices=np.where((ra>ramin) & (ra<ramax) & (dec>decmin) & (dec<decmax) & (zz>zmin) & (zz<zmax))

    #cutting and converting to math spherical coords
    phi=ra[indices]*deg2rad
    theta=pihalf - dec[indices]*deg2rad
    zz_cut=zz[indices]

    DC=cosmology.comoving_distance(zz_cut).value #comoving distances [Mpc]

    coords=np.empty((3,len(indices[0]))) #for output data
    coords[0,:]=DC*np.cos(phi)*np.sin(theta)
    coords[1,:]=DC*np.sin(phi)*np.sin(theta)
    coords[2,:]=DC*np.cos(theta)
    
    return coords



def Dist_sph(ra1,dec1,ra2,dec2,rr): #distance between 2 points on sphere ra,dec [rad] at radial dist rr 
    return rr*np.arccos(np.sin(ra1)*np.sin(ra2)+np.cos(ra1)*np.cos(ra2)*np.cos(dec2-dec1))



def Make_randoms(nobj): #creating randoms with nobj objects
    
    #radial ranges:
    DCmin=cosmology.comoving_distance(zmin).value+maxR
    DCmax=cosmology.comoving_distance(zmax).value-maxR

    DC_rand=np.random.uniform(DCmin,DCmax,nobj) #[Mpc]
    ra_rand=np.random.uniform(ramin*deg2rad,ramax*deg2rad,nobj) #[rad]

    #sin(dec) has to be uniform:
    sindec_min=min(np.sin(decmin*deg2rad),np.sin(decmax*deg2rad))
    sindec_max=max(np.sin(decmin*deg2rad),np.sin(decmax*deg2rad))
    sin_dec_rand=np.random.uniform(sindec_min,sindec_max,nobj)
    dec_rand=np.arcsin(sin_dec_rand) #[rad]
    
    #now cutting out dist<maxR edges:
    #declination of point on ra=ramin line which is closests to center of drawn circle
    dec_raminclose=np.arctan(np.tan(dec_rand)*np.cos(ramin*deg2rad - ra_rand)) #[rad]

	#declination of point on ra=ramax line which is closests to center of drawn circle
    dec_ramaxclose=np.arctan(np.tan(dec_rand)*np.cos(ramax*deg2rad - ra_rand)) #[rad]
    
    #dist between random point and [ramin,dec_raminclose]
    dist_1=Dist_sph(ra_rand,dec_rand,ramin*deg2rad*np.ones(nobj),dec_raminclose,DC_rand)
    #dist between random point and [ramax,dec_ramaxclose]
    dist_2=Dist_sph(ra_rand,dec_rand,ramax*deg2rad*np.ones(nobj),dec_ramaxclose,DC_rand)
    
    dist_decmin=Dist_sph(ra_rand, dec_rand, ra_rand, decmin*deg2rad*np.ones(nobj), DC_rand)
    dist_decmax=Dist_sph(ra_rand, dec_rand, ra_rand, decmax*deg2rad*np.ones(nobj), DC_rand)

    indices=np.where((dist_1>maxR) & (dist_2>maxR) & (dist_decmin>maxR) & (dist_decmax>maxR)) #far enough
    DC_rand_cut=DC_rand[indices]
    phi_rand_cut=ra_rand[indices]
    theta_rand_cut=pihalf - dec_rand[indices]
    
    coords=np.empty((3,len(indices[0]))) #for output data
    coords[0,:]=DC_rand_cut*np.cos(phi_rand_cut)*np.sin(theta_rand_cut)
    coords[1,:]=DC_rand_cut*np.sin(phi_rand_cut)*np.sin(theta_rand_cut)
    coords[2,:]=DC_rand_cut*np.cos(theta_rand_cut)
    
    return coords
    




#****************main****************
def main():
    #reading data
    if os.path.isfile(fname)!=True:
        print('File: %s does not exist:/' % (fname))
        return 0
    
    print('Reading data: %s' %(fname))
    RA,DEC,ZZ=np.genfromtxt(fname,usecols=(col_ra,col_dec,col_z),unpack=True)
    coords=Make_cuts(RA,DEC,ZZ)

    #making output
    fout=h5py.File(fname_out,'w')
    fout.create_dataset(dpos, data=coords, compression="gzip", compression_opts=9)
    fout.close()
    print('Successful conversion: %s -> %s' % (fname,fname_out))

    #creating randoms
    nrand=int(1.1*nmultiply*len(coords[0])) #number of randoms
    coords_rand=Make_randoms(nrand)
    foutr=h5py.File(fname_random,'w')
    foutr.create_dataset(dpos, data=coords_rand, compression="gzip", compression_opts=9)
    foutr.close()

    print('Successfuly created %d randoms in %s' % (len(coords_rand[0]),fname_random))
    return 0

main()