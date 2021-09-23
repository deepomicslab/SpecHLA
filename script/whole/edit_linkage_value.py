import sys

file = sys.argv[1]
value=sys.argv[2]
out=sys.argv[3]
o=open(out, 'w')
for line in open(file):
    array = line.strip().split()
    array[-1] = value
    new_line = ' '.join(array)
    print (new_line, file = o)
o.close()
