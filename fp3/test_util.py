from .util import even_blocks, index_blocks


def test_even_blocks():
    for b in even_blocks(0, 100, 10):
        assert b[1] - b[0] == 10

    for b in even_blocks(0, 2400, 72):
        assert b[1] - b[0] in (33, 34)

def test_index_blocks():
    assert index_blocks(0, 150, 0.1) == [(0, 150)]
    assert index_blocks(0, 3600, 0.1) == [(0, 50), (450, 500), (900, 950)]
    assert index_blocks(0, 900, 0.1) == [(0, 50), (425, 475), (850, 900)]
