import sys

def findnth(c,s):
    count = 0
    for i in range(len(s)-1,0,-1):
        if s[i] == 'X': count +=1
        if count == c:
            return i

pattern = "XX...........X...........................X.X.X...X.XX................X.XX.......X.XX.X...XXXX.X.XXXXXX"

suffix_len = int(sys.argv[1])
masterfile = sys.argv[2]

pattern = pattern[findnth(suffix_len,pattern):]

with open(masterfile) as mf, open(masterfile+".pats",'w') as out:
    for pat in mf:
        base = ["."] * len(pattern)
        offset = 1
        for i in range(len(pattern),0,-1):
            if pattern[i-1] == 'X':
                base[i-1] = pat[suffix_len-offset]
                offset += 1
        out.write("".join(base)+"\n")

