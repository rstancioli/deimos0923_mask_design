import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

# For the 2nd mask on any cluster, we don't need to check zphot and
# set priorities: for the slits we simply read in what dsim output
# when making the first mask, and echo the lines that reflect
# unselected slits. In general we DO still need to check the phot
# catalog to make sure the guide/alignment stars exist, as they may
# differ from the 1st mask....but not for this cluster

# user: don't forget to edit slitPA!

if len(sys.argv)!=2:
    sys.stderr.write('usage: mkmask2.py dsimoutputfile > newdsiminputfile\n')
    sys.exit(1)
infilename = sys.argv[1]    

maskRA = 120.453
maskDEC = 36.463
slitPA=-58
guideRA,guideDEC,guidemag = (120.45428,36.49947,13.48)
cosdec = np.cos(np.radians(guideDEC))
aligninfo = [[120.31574,36.50203,16.96],
                 [120.28645,36.47283,16.72],
                 [120.31067,36.47048,15.39],
                 [120.59353,36.42568,16.45],
                 [120.60123,36.45279,18.87],
                 [120.57664,36.45969,16.03]]


def raformatter(ra):
    ra /= 15 # degrees to hours
    hrs = int(ra)
    minutesleft = 60*(ra-hrs)
    minutes = int(minutesleft)
    seconds = 60*(minutesleft-minutes)
    return '%02d:%02d:%06.3f' % (hrs,minutes,seconds)

def decformatter(dec):
    isneg = dec<0 # save this fact for later
    dec = abs(dec)
    deg = int(dec)
    minutesleft = 60*(dec-deg)
    minutes = int(minutesleft)
    seconds = 60*(minutesleft-minutes)
    if isneg:
        return '-%02d:%02d:%05.2f' % (deg,minutes,seconds)
    else:
        return '+%02d:%02d:%05.2f' % (deg,minutes,seconds)


# get photometry, just to check guide/alignment stars
photdata = pd.read_csv('phot.csv')
rmag = photdata['rMeanKronMag'].values
photID = photdata['objID'].values
photRA = photdata['raMean'].values
photDEC = photdata['decMean'].values

# print the output hdr
print('# OBJNAME         RA          DEC        EQX   MAG band PCODE LIST SEL? PA L1 L2')
# start writing some ds9 region files for later inspection
unmatched = open('unmatched.reg','w')
mag999 = open('mag999.reg','w')
slits = open('slits2.reg','w')


# find the guide star using photometric cat only
distance = ((photRA-guideRA)*cosdec)**2 + (photDEC-guideDEC)**2
nearestobj = np.argmin(distance)
success = distance[nearestobj]<1e-3 and rmag[nearestobj]<17
if not success:
    sys.stderr.write('%f %f\n'%(distance[nearestobj],rmag[nearestobj]))
    sys.stderr.write('Did not find guide star, allowing user override\n')
    #sys.exit(1)
myRA = guideRA #photRA[nearestobj] # use phot catalog RA/DEC for...
myDEC = guideDEC #photDEC[nearestobj] # ...consistency w/slits!
rastr = raformatter(myRA) 
decstr = decformatter(myDEC) 
# but use user-supplied mag b/c may be saturated in catalog
name = str(photID[nearestobj])[-10:] # see below for name construction comment
# the '-1 0 1' means: -1 for guide star; 0 1 for use this dammit 
outstr = '%-15s %-13s %-12s 2000.0 %5.2f r -1 0 1' %  (name,rastr,decstr,guidemag)
print(outstr)
# special mark on ds9 region file
slits.write('fk5;box %f %f 0.002 0.002 # color=red\n'%(myRA,myDEC))
    
# also confirm the alignment stars
for thisRA,thisDEC,thismag in aligninfo:
    distance = ((photRA-thisRA)*cosdec)**2 + (photDEC-thisDEC)**2
    nearestobj = np.argmin(distance)
    success = distance[nearestobj]<1e-3 and rmag[nearestobj]<20.3
    if not success:
        sys.stderr.write('%f %f\n'%(distance[nearestobj],rmag[nearestobj]))
        sys.stderr.write('%f %f\n'%(thisRA,thisDEC))
        sys.stderr.write('Did not find alignment star, skipping\n')
        sys.exit(1)
    myRA = photRA[nearestobj] # use phot catalog RA/DEC for...
    myDEC = photDEC[nearestobj] # ...consistency w/slits!
    rastr = raformatter(myRA) 
    decstr = decformatter(myDEC) 
    # but use user-supplied mag b/c may be saturated in catalog
    name = str(photID[nearestobj])[-10:] # see below for name construction
    # the '-2 0 1' means: -2 for guide star; 0 1 for use this dammit 
    outstr = '%-15s %-13s %-12s 2000.0 %5.2f r -2 0 1' %  (name,rastr,decstr,thismag)
    print(outstr)
    # special mark on ds9 region file
    slits.write('fk5;box %f %f 0.002 0.002 # color=red\n'%(myRA,myDEC))

# now read the dsim output and echo the non-selected slits
infile = open(infilename,'r')
lines = infile.readlines()
echo = False
for line in lines:
    # skip down to "Non-Selected Objects"
    if echo:
        # filter out empty lines just for good form
        fields = line.split()
        if len(fields) > 5:
            print(line[:-1]) # remove the extra \n
            # also save in reg file for inspection
            # (w/THIS mask's alignment stars)
            ra = fields[1]
            dec = fields[2]
            priority = fields[6]
            slits.write('fk5;circle %s %s 0.001 # color=green text={%s}\n'%(ra,dec,priority))

    else:
        if line[0:14] == '# Non-Selected':
            echo = True
