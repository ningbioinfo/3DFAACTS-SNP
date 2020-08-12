# make allValidpairs into bedpe file
# breaks = 10000
import sys
import os
import pandas as pd

infile = sys.argv[1]
outfile = sys.argv[2]

breaks = int(sys.argv[3]) # resolution blocks (expecting numbers)



fchr1 = ''; fstart1 = ''; fend1 = ''; fchr2 = ''; fstart2 = ''; fend2 = ''; count = 0; interactions = []
with open(infile, 'r') as hic:
    for line in hic:
        data = line.strip().split('\t')[1:3] + line.strip().split('\t')[4:6]
        schr1 = data[0]
        pos1 = int(data[1])
        sstart1 = int(pos1 - breaks/2)
        send1 = sstart1 + breaks
        schr2 = data[2]
        pos2 = int(data[3])
        sstart2 = int(pos2 - breaks/2)
        send2 = sstart2 + breaks
        #print(schr1, sstart1, send1, schr2, sstart2, send2)
        # if any end on different chr, then they are different interactions
        if schr1 != fchr1 or schr2 != fchr2:
            if fchr1 != '':
                interactions.append([fchr1, fstart1, fend1, fchr2, fstart2, fend2, count])
                fchr1 = schr1; fstart1 = sstart1; fend1 = send1; fchr2 = schr2; fstart2 = sstart2; fend2 = send2
                count = 1
            else:# the evry first line
                fchr1 = schr1; fstart1 = sstart1; fend1 = send1; fchr2 = schr2; fstart2 = sstart2; fend2 = send2
                count = 1
        # if both end are in the same block, then they are the same, add count
        elif fstart1 <= pos1 <= fend1 and fstart2 <= pos2 <= fend2:
            count += 1
        else:# otherwise, take in and replace
            interactions.append([fchr1, fstart1, fend1, fchr2, fstart2, fend2, count])
            fchr1 = schr1; fstart1 = sstart1; fend1 = send1; fchr2 = schr2; fstart2 = sstart2; fend2 = send2
            count = 1


if os.path.exists(outfile):
    os.remove(outfile)


output = pd.DataFrame(interactions)
output.to_csv(outfile, sep = '\t', header=False, index=False)
