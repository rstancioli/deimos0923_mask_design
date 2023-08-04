import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import ascii, votable
import json
from types import SimpleNamespace

degtol = 0.001 # for matching targets w/archival z's

if len(sys.argv)==3:
    cluster_dir = sys.argv[1] + '/'
    mask_number = sys.argv[2]
    io_dir = cluster_dir + mask_number + '/'
else:
    print(sys.argv)
    sys.stderr.write('Usage: prioritize.py <cluster_dir> <mask_number> \n')
    sys.stderr.write(""""
        # to include special targets, add this to input.json:
        "specialtargets": [
        [
            0.948284,
            9.962229,
            15.765,
            5000
        ],
        [
            1.050969,
            10.170061,
            16.073,
            5000
        ]
        ]
        # where [226.96435,57.84257,18.02,5000] = [ra dec mag priority]
    """)
    sys.exit(1)

fp = "{}input.json".format(io_dir)
with open(fp, 'r') as f:
    json_data = json.load(f)
v = SimpleNamespace(**json_data)

guide = "{}guide.xml".format(io_dir)
try:
    votb = (votable.parse_single_table(source=guide)).to_table()
    v.guideRA,v.guideDEC,v.guidemag = list(votb.iterrows('_RAJ2000','_DEJ2000','rmag'))[0]
    sys.stderr.write("{} {} \n".format(v.guideRA, v.guidemag))
except:
    raise SystemExit('Invalid guide file.')

align = "{}align.xml".format(io_dir)
try:
    votb2 = (votable.parse_single_table(source=align)).to_table()
    v.aligninfo = list(votb2.iterrows('_RAJ2000','_DEJ2000','rmag'))
except:
    raise SystemExit('Invalid align file.')

v.cosdec = np.cos(np.radians(v.guideDEC))
print(v)

### OLD INPUT METHOD
# # RMJ 1219: USER MUST CHANGE THESE
# v.maskRA = 184.723 
# v.maskDEC = 50.900 
# v.zclust=0.535
# v.slitPA=32
# v.BCGmag = 21.23
# v.guideRA,v.guideDEC,v.guidemag = (184.82564,50.87693,16.62)
# v.cosdec = np.cos(np.radians(v.guideDEC))
# v.aligninfo = [[184.94876,50.91130,16.08], # maybe off SE corner, but great if it's on!
#                  [184.91292,50.92500,18.04],
#                  [184.88859,50.92163,17.09], 
#                  [184.82199,50.96310,15.41],#possibly in E vignetted area
#                  [184.54643,50.85060,17.12],
#                  [184.54997,50.87731,16.06]]#possibly in W vignetted area
# # RMJ 0003: USER MUST CHANGE THESE
# v.maskRA = 0.96952219625
# v.maskDEC = 10.06685767
# v.zclust=0.37
# v.slitPA=-2
# v.BCGmag = 18.288
# v.guideRA,v.guideDEC,v.guidemag = (001.08345800, +10.07273700, 16.396)
# v.cosdec = np.cos(np.radians(v.guideDEC))
# v.aligninfo = [[000.94828400, +09.96222900, 15.765],
#             [001.05096900, +10.17006100, 16.073],
#             [001.09680700, +10.13814800, 16.484],
#             [001.06405100, +10.16007700, 18.052],
#             [000.90873600, +09.97891300, 18.173],
#             [001.06154100, +10.17188100, 19.690]]
# # Special Targets 
# v.specialtargets = []
# example:   [226.96435,57.84257,18.02,5000],# ra dec mag priority

                
# # load known redshifts
# knowndata = np.loadtxt('archivalz.txt')
# knownra = knowndata[:,0]
# knowndec = knowndata[:,1]

# load known redshifts (new_method)
knowndata = ascii.read('{}archivalz.csv'.format(cluster_dir))
knownra = list(knowndata['RA'].data)
knowndec = list(knowndata['DEC'].data)

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
    
# get redshifts and likelihood of being in cluster
df = pd.read_csv('{}zphot.csv'.format(cluster_dir))
goodmask = df['z_phot0']>-999
z = df['z_phot0'][goodmask].values
zerr = df['z_photErr'][goodmask].values
# v. few of these, prevent div by zero
zerr[zerr==0] = 0.5
zID = df['objID'][goodmask].values
ra = df['raMean'][goodmask].values
dec = df['decMean'][goodmask].values
Lclust = np.exp(-0.5*(z-v.zclust)**2/zerr**2)/zerr
#print(Lclust)
#plt.hist(Lclust,bins=50)
#plt.show()


# get photometry, to prioritize bright objecs
# Using mean catalog
photdata = pd.read_csv('{}phot.csv'.format(cluster_dir))
rmag = photdata['rMeanKronMag'].values
# # Using stacked catalog
# photdata = pd.read_csv('phot_stack.csv')
# rmag = photdata['rKronMag'].values

photID = photdata['objID'].values
photRA = photdata['raMean'].values
photDEC = photdata['decMean'].values

# start writing some ds9 region files for later inspection
unmatched = open('{}unmatched.reg'.format(io_dir),'w')
mag999 = open('{}mag999.reg'.format(io_dir),'w')
slits = open('{}targets.reg'.format(io_dir),'w')

# print the output hdr
targets = open('{}targets.txt'.format(io_dir),'w')
targets.write('# OBJNAME         RA          DEC        EQX   MAG band PCODE LIST SEL? PA L1 L2\n')
#excludefile = open('knownz.reg','w')


# find the guide star using photometric cat only
distance = ((photRA-v.guideRA)*v.cosdec)**2 + (photDEC-v.guideDEC)**2
nearestobj = np.argmin(distance)
success = distance[nearestobj]<1e-3 and rmag[nearestobj]<17
if not success:
    sys.stderr.write('%f %f\n'%(distance[nearestobj],rmag[nearestobj]))
    sys.stderr.write('%f %f\n'%(photRA[nearestobj],photDEC[nearestobj]))
    sys.stderr.write('Did not find guide star in PS phot cat, so Im trusting your info!!\n')
    #sys.exit(1)
myRA = v.guideRA #photRA[nearestobj] # use phot catalog RA/DEC for...
myDEC = v.guideDEC #photDEC[nearestobj] # ...consistency w/slits!
rastr = raformatter(myRA) 
decstr = decformatter(myDEC) 
# use user-supplied mag b/c may be saturated in catalog
name = str(photID[nearestobj])[-10:] # see below for name construction comment
# the '-1 0 1' means: -1 for guide star; 0 1 for use this dammit 
outstr = '%-15s %-13s %-12s 2000.0 %5.2f r -1 0 1\n' %  (name,rastr,decstr,v.guidemag)
targets.write(outstr)
# special mark on ds9 region file
slits.write('fk5;box %f %f 0.002 0.002 # color=red\n'%(myRA,myDEC))
    
# also confirm the alignment stars
for thisRA,thisDEC,thismag in v.aligninfo:
    distance = ((photRA-thisRA)*v.cosdec)**2 + (photDEC-thisDEC)**2
    nearestobj = np.argmin(distance)
    success = distance[nearestobj]<1e-3 and rmag[nearestobj]<21.2
    if not success:
        sys.stderr.write('%f %f\n'%(distance[nearestobj],rmag[nearestobj]))
        sys.stderr.write('%f %f\n'%(thisRA,thisDEC))
        sys.stderr.write('Did not find alignment star, skipping\n')
        continue
    myRA = photRA[nearestobj] # use phot catalog RA/DEC for...
    myDEC = photDEC[nearestobj] # ...consistency w/slits!
    rastr = raformatter(myRA) 
    decstr = decformatter(myDEC) 
    # but use user-supplied mag b/c may be saturated in catalog
    name = str(photID[nearestobj])[-10:] # see below for name construction
    # the '-2 0 1' means: -2 for alignment star; 0 1 for use this dammit
    outstr = '%-15s %-13s %-12s 2000.0 %5.2f r -2 0 1\n' %  (name,rastr,decstr,thismag)
    targets.write(outstr)
    # special mark on ds9 region file
    slits.write('fk5;box %f %f 0.002 0.002 # color=red\n'%(myRA,myDEC))


# special targets.
nspecial = 0
for thisRA,thisDEC,thismag,thispri in v.specialtargets:
    # deprecated: look for this targ in PS catalog
    #distance = ((photRA-thisRA)*v.cosdec)**2 + (photDEC-thisDEC)**2
    #nearestobj = np.argmin(distance)
    #success = distance[nearestobj]<1e-3 and rmag[nearestobj]<20
    #if not success:
    #    sys.stderr.write('%f %f\n'%(distance[nearestobj],rmag[nearestobj]))
    #    sys.stderr.write('%f %f\n'%(thisRA,thisDEC))
    #    sys.stderr.write('Did not find special target, DOING IT ANYWAY\n')
    #myRA = photRA[nearestobj] # use phot catalog RA/DEC for...
    #myDEC = photDEC[nearestobj] # ...consistency w/slits!
    ## but override if given special code..
    #if thispri == 5001:
    myRA = thisRA
    myDEC = thisDEC
    rastr = raformatter(myRA) 
    decstr = decformatter(myDEC) 
    # use user-supplied mag b/c may be saturated in catalog
    name = 'special%d' % (nspecial)
    nspecial += 1
    # old:name=str(photID[nearestobj])[-10:] # see below for name construction
    outstr = '%-15s %-13s %-12s 2000.0 %5.2f r %d 1 0 %d\n' %  (name,rastr,decstr,thismag,thispri,v.slitPA)
    targets.write(outstr)
    # add to ds9 region file
    slits.write('fk5;circle %f %f 0.001 # color=green text={%d}\n'%(myRA,myDEC,thispri))


# now match the 2 galaxy cats, starting w/highest priority objects, and
# output the matched info
nslitcandidates = 0
for i in range(len(Lclust)):
    priority = Lclust[i]
    mag = 24 # default if missing
    if priority<0.01:
        continue
    # Upgrade TBD: match based on RA/DEC?
    domatch = photID==zID[i]
    if domatch.sum() != 1:
        # No match: most likely outside the footprint of the phot download.
        # Keep a record for inspection
        # but do NOT go on to write a slit candidate
        unmatched.write('fk5;circle %f %f 0.001 # color=red\n'%(ra[i],dec[i]))
        continue;
    j = domatch.argmax()
    mag = rmag[j]
    if mag <-900:
        # too faint a target. Keep a record for inspection
        # but do NOT go on to write a slit candidate
        mag999.write('fk5;circle %f %f 0.001 # color=yellow\n'%(ra[i],dec[i]))
        continue

    # if we got here we have a decent target, but let's reject it
    # if it has an archival redshift AND is near ctr of field
    myRA = photRA[j]   # use phot cat coords for consistency w/stars!
    myDEC = photDEC[j]
    # check for near ctr of field
    distfromctr = np.sqrt(((myRA-v.maskRA)*v.cosdec)**2+(myDEC-v.maskDEC)**2)
    if distfromctr < 0.083: # 0.083 deg = 5' = 1 chip
        # check for match w/known z
        distfromknown = np.sqrt(((myRA-knownra)*v.cosdec)**2+(myDEC-knowndec)**2)
        if distfromknown.min() < degtol:
            sys.stderr.write('skipping known z at %f %f\n' % (myRA,myDEC))
            continue
    
    # if we got here we have a slit candidate.
    nslitcandidates += 1
    # prioritize further based on the magnitude
    priority *= 24-np.clip(mag,v.BCGmag,None)
    # scale the priority to better match the documented example
    priority *= 100
    # exclusions...this way of doing it is deprecated
    #if int(priority) in exclude:
    #    sys.stderr.write('excluding %d\n' % (int(priority)))
    #    # add to ds9 region file
    #    excludefile.write('fk5;circle %f %f 0.001 # color=red text={%d}\n'%(myRA,myDEC,priority))
    #    continue
    # preselections: also deprecated
    #preselect = 0
    #if int(priority) in preselections:
    #    sys.stderr.write('preselecting %d\n' % (int(priority)))
    #    preselect = 1
        
    # come up with a 16-char-max name for this obj. For Pan-STARRS,
    # use the last 10 digits of the objid (the 1st 8 seem to be the
    # same for a given area of sky)
    name = str(zID[i])[-10:]
    # adjustments to slit length? *****TBD******
    outstr = '%-15s %-13s %-12s 2000.0 %5.2f r %d 1 0 %d\n' %  (name,raformatter(myRA),decformatter(myDEC),mag,priority,v.slitPA)
    color = 'green'   
    targets.write(outstr)
    # Add to region file for later inspection
    slits.write('fk5;circle %f %f 0.001 # color=%s text={%d}\n'%(myRA,myDEC,color,priority))

if nslitcandidates <100:
    sys.stderr.write('WARNING: few slit candidates, something probably went wrong.\n')

# add sky slits
#outstr = '01010101 15:09:00 +57:59:22 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010102 15:09:03.28 +57:59:22.4 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010103 15:08:52.46 +57:59:31.4 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010104 15:08:44.698 +57:57:39.13 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)

