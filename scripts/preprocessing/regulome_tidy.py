import sys

DATABASE_DIR= sys.argv[1]

f = open(DATABASE_DIR + '/Regulome.bed','r')
g = f.read().splitlines()
f.close()

k = open(DATABASE_DIR + '/Regulome_new.bed','w')
for line in g:
    line = line.split(' ')
    k.write('\t'.join([line[0],line[1],line[2], ' '.join(line[3:])]))
k.close()
