def even_blocks(n0, n1, m):
    """Split the range(n0, n1) into m evenly sized chunks, yeild start, end
    for each chunk."""

    r = m
    n = n1 - n0

    while r > 0:
        s = int(round(n / r))
        yield n0, n0 + s
        n -= s
        n0 += s
        r -= 1

def index_blocks(n0, n1, osc):
    """Give back a list of blocks in the range n0, n1 of around 5°, spaced at 
    the start, start + 45°, start + 90° if possible. If <= 15° of data in
    total, use all, if < 95° get as close as possible to 45° wedge spacing."""

    n = n1 - n0
    five = int(round(5 / osc))

    if n * osc <= 15:
        return [(n0, n1)]
    
    elif n * osc > 95:
        n45 = int(round(45 / osc))
        n90 = int(round(90 / osc))
        return [(n0, n0 + five), (n45, n45 + five), (n90, n90 + five)]

    else:
        half = n0 + (n - five) // 2
        return [(n0, n0 + five), (half, half + five), (n1 - five, n1)]

    
    
