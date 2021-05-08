import sys
import gzip

def read_fq(file, outfile):
    saved_name = {}
    out = gzip.open(outfile, 'wt')
    line_num = 0
    for line in gzip.open(file, 'rt'):
        line = line.strip()
        #if line[0] == '@':
        if line_num % 4 == 0:
            name = line[1:-2]
            if name in saved_name.keys():
                saved_name[name] += 1
                name = name + '_' + str(saved_name[name])
                line = line[0] + name + line[-2:]
            else:
                saved_name[name] = 1
            print (line, file= out)
            #print (line, name)
        #break
        else:
            print (line, file= out)
        line_num += 1
    out.close()

read_fq(sys.argv[1], sys.argv[2])
