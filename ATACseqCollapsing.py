## This is a script to collapse redundant ATAC-seq peaks from multiple samples.
## Before using this script, you should CONCATENATE the samples you want and sort them by chromosomes, then positions
## To do that, you can either using bash command `sort` or import the summits.bed into R with `import.bed` function from
## the rtracklayer package, then use `c()` and follow by `sort()` to get it sorted ,then you can export the aggregated bed file.
## Author: Ning Liu
## Date: 13 Jun 2019
import sys

## sys.argv[1] would be the input atacseq bed file
## sys.argv[2] would be the output bed file name

with open(sys.argv[1], 'r') as atacinput:
    for getfirstline in atacinput: # get the first line data.
        data = getfirstline.strip().split('\t')
        c = data[0]
        s = int(data[1])-499
        e = int(data[2])+500
        n = data[3]
        score = float(data[4])
        break

    collapsing = []
    scores = [score]
    for line in atacinput:
        data2 = line.strip().split('\t')
        c2 = data2[0]
        s2 = int(data2[1])
        e2 = int(data2[2])
        n2 = data2[3]
        score2 = float(data2[4])
        # decide to collapse with the previous range or not.
        if c2 == c and s <= s2 <= e:
            scores.append(score2) # save all the score so we can pick the max score to represent the block
            continue
        else:
            collapsing.append([c,s,e,n,max(scores)])
            c = c2
            s = s2 - 499
            e = e2 + 500
            n = n2
            score = score2

            scores = [score]
    collapsing.append([c,s,e,n,max(scores)])

with open(sys.argv[2], 'a') as outputbed:
    for i in collapsing:
        print("\t".join(str(j) for j in i), file = outputbed)
