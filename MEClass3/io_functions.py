import sys
import os

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_tup(inp_tup, noftf, seperation, end_chr):
        noftf.write(seperation.join(str(dat) for dat in inp_tup) + end_chr)
        
def mk_output_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def print_to_log(FH, *args, **kwargs):
    FH.write(*args, **kwargs)
    FH.flush()