from __future__ import print_function
from collections import OrderedDict

from bowler import *
from fissix.fixer_util import Attr, Name, ArgList, String, Comma
from fissix.pytree import type_repr, Node, Leaf  # , Node, type_repr
from bowler.types import TOKEN
from fissix.pgen2 import token  # token.COMMA
from fissix.pygram import python_symbols as syms

from libfuturize import fixer_util

import sys


def print_node(node, indent="", last=True):
    """Debugging function to print node tree"""
    if last:
        first_i = "└─"
        second_i = "  "
    else:
        first_i = "├─"
        second_i = "│ "
    prefix = indent + first_i
    if type(node) is Node:
        print(
            prefix
            + "Node[{}] prefix={} suffix={}".format(
                type_repr(node.type), repr(node.prefix), repr(node.get_suffix())
            )
        )
    else:
        print(indent + first_i + repr(node))
    indent = indent + second_i

    children = list(node.children)
    for i, child in enumerate(node.children):
        print_node(child, indent, last=(i + 1) == len(children))


def rewrite_print(node, capture, filename):
    """Given a print node, rewrite to a logger.debug"""

    params = capture["function_parameters"]
    # Args is either a Leaf or an arglist Node here
    arglist = params.children[1]

    # Extract the arguments from this list
    if isinstance(arglist, Node) and arglist.type == syms.arglist:
        args = [x for x in arglist.children if not x.type == token.COMMA]
        # Remove kwargs for now
        non_kwargs = [x for x in args if not x.type == syms.argument]
        # Read out things like sep here
        sep = " "

        if not len(non_kwargs) == len(args):
            first_leaf = next(node.leaves())
            print(
                "Warning: {}:{} Not handling keyword argument loggers".format(
                    filename, first_leaf.lineno
                )
            )
            return None

        # If more than one child, we need to construct a new joiner string
        if len(non_kwargs) > 1:
            # Instead of placing a new, modify the first if it's a string
            if arglist.children[0].type == token.STRING:
                value = arglist.children[0].value
                last_char = value[-1]
                new_end = sep + sep.join(["%s"] * (len(non_kwargs) - 1))
                arglist.children[0].value = value[:-1] + new_end + last_char
            else:
                arglist.insert_child(
                    0, String('"' + sep.join(["%s"] * len(non_kwargs)) + '"')
                )
                arglist.insert_child(1, Comma())
                arglist.children[2].prefix = " " + arglist.children[2].prefix

    # Use the possibly rewritten parameters in the new call
    new_node = Attr(Name("logger"), Name("debug"))
    new_node[0].prefix = node.children[0].prefix
    node.children = new_node + [params]


if __name__ == "__main__":
    path = sys.argv[1]

    (
        Query(path)
        .select_function("print")
        .modify(rewrite_print)
        .execute(silent=True, write=True)
    )

    # main()
