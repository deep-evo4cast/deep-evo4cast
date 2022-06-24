pattern = "XX...........X...........................X.X.X...X.XX................X.XX.......X.XX.X...XXXX.X.XXXXXX"

for i ,c in enumerate(pattern):
    if c == 'X': print(f'sed s/././{i+1} last18-matches.pats > last18-matches-wc-{i+505+1-len(pattern)}.pats')

