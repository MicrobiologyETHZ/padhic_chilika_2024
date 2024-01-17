

import sys

infile = sys.stdin
while True:
    header = infile.readline()
    if not header:
        break
    splits = header.strip().split()
    header = splits[0] + ' ' + splits[1][2:] + '/' + splits[1][0]
    
    seq = infile.readline().strip()
    plus = infile.readline().strip()
    plus = '+'
    qual = infile.readline().strip()
    print(header)
    print(seq)
    print(plus)
    print(qual)
    
