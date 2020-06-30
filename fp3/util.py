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
