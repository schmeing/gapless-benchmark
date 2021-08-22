#!/usr/bin/env python3
from Bio import SeqIO

import gzip
import sys

def usage():
    print("create_testset.py gaps/inversions {input}.fa {output}.fa")

def main(argv):
    if 0 == len(argv):
        usage()
        sys.exit(1)
    
    if 3 != len(argv):
        print("Wrong number of arguments")
        usage()
        sys.exit(1)
        
    if argv[0] not in ["gaps","inversions"]:
        print("Unknown testset option: ", argv[0])
        usage()
        sys.exit(1)
        
    with open(argv[2], 'w') as fout:
        with gzip.open(argv[1], 'rb') if 'gz' == argv[1].rsplit('.',1)[-1] else open(argv[1], 'rU') as fin:
            for record in SeqIO.parse(fin, "fasta"):
                # Modify record
                seq = []
                interval = 75000
                end_pos = 0
                for n, cur_pos in enumerate(range(interval, len(record.seq)-interval, interval)):
                    if n < 50:
                        seq.append(str(record.seq[end_pos:cur_pos]))
                        if n < 10:
                            end_pos = cur_pos + 500
                        elif n < 20:
                            end_pos = cur_pos + 2000
                        elif n < 30:
                            end_pos = cur_pos + 5000
                        elif n < 40:
                            end_pos = cur_pos + 20000
                        else:
                            end_pos = cur_pos + 50000
                        
                        if "gaps" == argv[0]:
                            seq.append("N"*(end_pos-cur_pos))
                        elif "inversions" == argv[0]:
                            seq.append(str(record.seq[cur_pos:end_pos].reverse_complement()))

                seq.append(str(record.seq[end_pos:]))
                
                # Write out modified record
                fout.write(">{}_{}\n".format(record.description, argv[0]))
                fout.write(''.join(seq))
                fout.write('\n')

if __name__ == "__main__":
    main(sys.argv[1:])
