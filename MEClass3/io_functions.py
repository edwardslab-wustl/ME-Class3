import sys

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_tup(inp_tup, noftf, seperation, end_chr):
        noftf.write(seperation.join(str(dat) for dat in inp_tup) + end_chr)