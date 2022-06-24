pattern = "XX...........X...........................X.X.X...X.XX................X.XX.......X.XX.X...XXXX.X.XXXXXX"
cX = 0

for i ,c in enumerate(pattern):
    cX += (c == 'X')
    if c == 'X' and cX > 9 : print(f'sed s/././{i+1} 18-{i+505+1-len(pattern)}/iterlast/full-dense.pats > 18-{i+505+1-len(pattern)}/iterlast/full-wc.pats')

