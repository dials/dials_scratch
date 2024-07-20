import functools
from .util import even_blocks, index_blocks, prime_factors


def test_even_blocks():
    for b in even_blocks(0, 100, 10):
        assert b[1] - b[0] == 10

    for b in even_blocks(0, 2400, 72):
        assert b[1] - b[0] in (33, 34)


def test_index_blocks():
    assert index_blocks(0, 150, 0.1) == [(0, 150)]
    assert index_blocks(0, 3600, 0.1) == [(0, 50), (450, 500), (900, 950)]
    assert index_blocks(0, 900, 0.1) == [(0, 50), (425, 475), (850, 900)]


def test_prime_factors():
    assert prime_factors(12) == [2, 2, 3]
    assert prime_factors(2 ** 8) == [2] * 8
    assert prime_factors(2 * 3 * 4 * 5 * 6) == [2, 2, 2, 2, 3, 3, 5]

    # FIXME should probably move this to it's own function
    f = prime_factors(12)
    a = functools.reduce((lambda x, y: x * y), f[0::2])
    b = functools.reduce((lambda x, y: x * y), f[1::2])
    assert a == 6
    assert b == 2
