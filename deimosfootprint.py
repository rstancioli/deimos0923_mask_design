import sys
import numpy as np
from astropy.coordinates import SkyCoord, FK5
from astropy import units as u

# setup deimos geometry
myscale = 5/60/55 # deg/mm
chipwidth = 53*myscale 
chipheight = 44*myscale 
chipy = np.array([48,2,-44,-90])*myscale + chipheight/2
chipx = 10*myscale
guiderx = np.array([-25.5,-6.5])*myscale
guidery = 2*myscale
guiderwidth = np.array([15,23])*myscale
guiderheight = 39*myscale

def xy2sky(dx,dy,pa,dec):
    #print(dx,dy,pa,dec)
    cospa = np.cos(np.radians(pa))
    sinpa = np.sin(np.radians(pa))
    dra = -(dx*cospa - dy*sinpa)/np.cos(np.radians(dec))
    ddec = dx*sinpa + dy*cospa
    #print(dra,ddec)
    return dra,ddec

if len(sys.argv)==3:
    sc = SkyCoord(sys.argv[1], unit=(u.hourangle, u.deg))
    ra = sc.ra.deg
    dec = sc.dec.deg
    pa = float(sys.argv[2])
    # print("sexagesimal")
elif len(sys.argv)==4:
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    pa = float(sys.argv[3])
    # print("decimal")
else:
    print(sys.argv)
    sys.stderr.write('Usage: deimosfootprint.py ra dec PA\n OR \n deimosfootprint.py "00:03:50.5 +10:02:06.2" PA')
    sys.exit(1)
print('# ra dec pa: ',ra,dec,pa)

# main footprint
for chip in [0,1,2,3]:
    dy = chipy[chip]
    dra,ddec = xy2sky(chipx,dy,pa,dec)
    print('fk5;box %f %f %f %f %.1f # dashlist=8 3' % (ra+dra,dec+ddec,chipwidth,chipheight,pa))

# chip 3: show vignetted region
dra,ddec = xy2sky(chipx+chipwidth/2,chipy[3],pa,dec)
ra1 = ra+dra
dec1 = dec+ddec
dra,ddec = xy2sky(chipx,chipy[3]-chipheight/2,pa,dec)
ra2 = ra+dra
dec2 = dec+ddec
print('fk5;line %f %f %f %f # dash=1' % (ra1,dec1,ra2,dec2))

# chip 0: show vignetted region
dra,ddec = xy2sky(chipx+chipwidth/2,chipy[0],pa,dec)
ra1 = ra+dra
dec1 = dec+ddec
dra,ddec = xy2sky(chipx,chipy[0]+chipheight/2,pa,dec)
ra2 = ra+dra
dec2 = dec+ddec
print('fk5;line %f %f %f %f' % (ra1,dec1,ra2,dec2))

# central vignetted region
dra,ddec = xy2sky(chipx+chipwidth/2,chipy[1]+chipheight/2,pa,dec)
ra1 = ra+dra
dec1 = dec+ddec
dra,ddec = xy2sky(chipx+0.35*chipwidth,0,pa,dec)
ra2 = ra+dra
dec2 = dec+ddec
dra,ddec = xy2sky(chipx+chipwidth/2,chipy[2]-chipheight/2,pa,dec)
ra3 = ra+dra
dec3 = dec+ddec
print('fk5;line %f %f %f %f %f %f # dash=1' % (ra1,dec1,ra2,dec2,ra3,dec3))


# guider boxes
for guider in [0,1]:
    dx = guiderx[guider]
    dy = guidery+guiderheight/2
    dra,ddec = xy2sky(dx,dy,pa,dec)
    print('fk5;box %f %f %f %f %.1f' % (ra+dra,dec+ddec,guiderwidth[guider],guiderheight,pa))

print('fk5;point %f %f' % (ra, dec))
        
