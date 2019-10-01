from sys import argv

script, bim, fam, snps, inds = argv

#opens up bim file, reads each line and pulls out CHR and BP of each, separates by tab. Writes each to specified file
with open(bim, 'r') as b:
    stringlist = [line.split()[0] + '\t' + line.split()[3] for line in b]
with open(snps, 'w') as outfile:
    [outfile.write(x + '\n') for x in stringlist]

#opens fam file, reads each line and pulls out individual id.  Writes each line to specified file
with open(fam, 'r') as f:
    famlist = [("_").join(["1",line.split()[1]]) for line in f]
with open(inds, 'w') as outfile:
    [outfile.write(x + '\n') for x in famlist]
