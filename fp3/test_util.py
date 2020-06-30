from .util import even_blocks


def test_even_blocks():
    for b in even_blocks(0, 100, 10):
        assert b[1] - b[0] == 10

    for b in even_blocks(0, 2400, 72):
        assert b[1] - b[0] in (33, 34)
