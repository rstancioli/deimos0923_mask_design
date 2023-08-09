import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from astropy.io import ascii, votable
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import norm
import json
from types import SimpleNamespace

degtol = 0.001 # for matching targets w/archival z's

def project_onto_merger_axis(coord, center, pa=38):
    """
    coord is a SkyCoord object of the point to be projected.
    center is a SkyCoord object of the origin.
    pa is the angle (east of north) that defines the axis
    """
    ra, dec = coord.ra.deg, coord.dec.deg
    ra_cen, dec_cen = center.ra.deg, center.dec.deg
    # cosdec = np.cos(np.radians(dec.mean()))
    cosdec = np.cos(np.radians(dec_cen))
    y = dec - dec_cen
    x = (ra - ra_cen)*cosdec
    pa = np.radians(pa)
    newx = (np.cos(pa)*x - np.sin(pa)*y)*60 # deg->arcmin
    newy = (np.sin(pa)*x + np.cos(pa)*y)*60 # deg->arcmin
    return newx, newy

def project_list(wcs_list, center, pa):
    x_list = []
    y_list = []
    for point in wcs_list:
        x, y = project_onto_merger_axis(point, center, pa)
        x_list.append(x)
        y_list.append(y)
    return np.array(x_list), np.array(y_list)

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

def read_target_file(io_dir, filename):
    with open('{}{}'.format(io_dir, filename),'r') as targets:
        lines = targets.readlines()

    align = []
    targets = []
    for line in lines:
        line_data = line.strip().split()
        if len(line_data)==9:
            align.append(line_data)
        elif len(line_data)==10:
            targets.append(line_data)
    return align, targets

# prepare column names to read and write csv files
columns = ['OBJNAME', 'RA', 'DEC', 'EQX', 'MAG', 'band', 'PCODE', 'LIST', 'SEL?', 'PA', 'L1', 'L2']
dtype = []
for column in columns:
    if column=='SEL?':
        dtype.append('i4')
    else:
        dtype.append('<U12')

# Read input
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

# load known redshifts (new_method)
knowndata = ascii.read('{}archivalz.csv'.format(cluster_dir))
knownra = list(knowndata['RA'].data)
knowndec = list(knowndata['DEC'].data)

# load selected targets in other masks
selected = []
# selectedra = []
# selecteddec = []
if v.previous_masks:
    for mask_path in v.previous_masks:
        _, targets = read_target_file('', mask_path)
        df_targets = Table(np.array(targets), names=columns[:-2], dtype=dtype[:-2])
        selected += list(df_targets[df_targets['SEL?']==1]['OBJNAME'])
    
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

# calculate distance from axis
# kernel_width is the standard deviation of the gaussian around the merger axis in arcmin.
# The f_dist factor will multiply the priority.
if v.axis['kernel_width']:
    sys.stderr.write('Prioritizing cluster axis.\n')
    sc_list =  SkyCoord(ra=ra, dec=dec, unit=u.deg)
    cl_center = SkyCoord(v.axis['ra_center'], v.axis['dec_center'], unit=u.degree)
    dist_axis = project_list(sc_list, cl_center, v.axis['axisPA'])
    f_dist = norm.pdf(dist_axis[0], scale=v.axis['kernel_width'])
else:
    sys.stderr.write('Not prioritizing cluster axis.\n')

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
# careful not to repeat special targets! This is not being checked automatically
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
n_skipped = 0
n_excluded = 0
for i in range(len(Lclust)):
    # come up with a 16-char-max name for this obj. For Pan-STARRS,
    # use the last 10 digits of the objid (the 1st 8 seem to be the
    # same for a given area of sky)
    name = str(zID[i])[-10:]

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

    # If this is mask2, reject targets from mask1
    if name in selected:
        # print('Target has been excluded!!!!')
        n_skipped += 1
        continue

    if name in v.exclude_target:
        # print('Target has been excluded!!!!')
        n_excluded += 1
        continue

    # if we got here we have a decent target, but let's reject it
    # if it has an archival redshift AND is near ctr of field
    myRA = photRA[j]   # use phot cat coords for consistency w/stars!
    myDEC = photDEC[j]
    # check for near ctr of field
    distfromctr = np.sqrt(((myRA-v.maskRA)*v.cosdec)**2+(myDEC-v.maskDEC)**2)
    if distfromctr < 0.083: # 0.083 deg = 5' = 1 chip
    # if distfromctr < 0.12: # 1.5 chips
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

    if v.axis['kernel_width']:
        priority *= f_dist[i]
        
    # adjustments to slit length? *****TBD******
    outstr = '%-15s %-13s %-12s 2000.0 %5.2f r %d 1 0 %d\n' %  (name,raformatter(myRA),decformatter(myDEC),mag,priority,v.slitPA)
    color = 'green'   
    targets.write(outstr)
    # Add to region file for later inspection
    slits.write('fk5;circle %f %f 0.001 # color=%s text={%d}\n'%(myRA,myDEC,color,priority))

sys.stderr.write('skipped {} targets selected in previous masks\n'.format(n_skipped))
sys.stderr.write('excluded {} targets manually\n'.format(n_excluded))


if nslitcandidates <100:
    sys.stderr.write('WARNING: few slit candidates, something probably went wrong.\n')

# Create target catalog for Aladin

targets.close()
unmatched.close()
mag999.close()
slits.close()

align, targets = read_target_file(io_dir, 'targets.txt')

### IMPORTANT
# If L1, L2 were used, we have to change:
#   (i) the 'names' argument to match the number of columns
#   (ii) create a function like read_target_file with the right conditions
df_targets = Table(np.array(targets), names=columns[:-2])
df_align = Table(np.array(align), names=columns[:-3])
# Sometimes this doesn't work, so we have to try this. I have to fix it!
# df_align = Table(np.array(align[:-1]), names=columns[:-3])

df_targets.write('{}targets.csv'.format(io_dir), format='ascii.csv', overwrite=True)
df_align.write('{}align.csv'.format(io_dir), format='ascii.csv', overwrite=True)

# add sky slits
#outstr = '01010101 15:09:00 +57:59:22 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010102 15:09:03.28 +57:59:22.4 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010103 15:08:52.46 +57:59:31.4 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)
#outstr = '01010104 15:08:44.698 +57:57:39.13 2000.0 20.00 r 5000 1 1 %d' %  (v.slitPA)
#targets.write(outstr)

