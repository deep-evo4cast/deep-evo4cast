
def findnth(c,s):
    count = 0
    for i in range(len(s)-1,0,-1):
        if s[i] == 'X': count +=1
        if count == c:
            return i

pattern = "XX...........X...........................X.X.X...X.XX................X.XX.......X.XX.X...XXXX.X.XXXXXX"

suffix_len = 18
pattern = pattern[findnth(suffix_len,pattern):]

with open(f'all.pats') as final, open(f'all-dense.pats','w') as out:
    for pat in final:
        offset = 0
        base = [None for _ in range(suffix_len)]
        for i in range(len(pattern)):
            if pattern[i] == 'X':
                base[offset] = pat[i]
                offset += 1
        out.write("".join(base)+"\n")

