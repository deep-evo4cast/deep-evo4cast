def findnth(c,s):
    count = 0
    for i in range(len(s)-1,0,-1):
        if s[i] == 'X': count +=1
        if count == c:
            return i

pattern = "XX...........X...........................X.X.X...X.XX................X.XX.......X.XX.X...XXXX.X.XXXXXX"

suffix_len = 18
pattern = pattern[findnth(suffix_len,pattern):]

for i ,c in enumerate(pattern):
    if c == 'X': print(f'sed s/././{i+1} last18.pats |sort | uniq > 18-{i+505+1-len(pattern)}/last18.pats')


