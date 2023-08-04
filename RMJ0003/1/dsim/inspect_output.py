# now read the dsim output and echo the non-selected slits
infilename = "output.txt"
outfilename = "non_selected.reg"

slits = open(outfilename, 'w')
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