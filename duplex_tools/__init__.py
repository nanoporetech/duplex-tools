"""Duplex Sequencing Tools package."""
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from duplex_tools import \
    assess_split_on_adapter, split_on_adapter  # noqa: F401

modules = ["split_on_adapter", "assess_split_on_adapter"]

__version__ = '0.1.3'


def main():
    """Entry point."""
    parser = ArgumentParser(
        "Duplex Sequencing Tools",
        formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(
        title="subcommands", description="valid commands",
        help="additional help", dest="command")
    subparsers.required = True
    for module in modules:
        mod = globals()[module]
        p = subparsers.add_parser(module, parents=[mod.argparser()])
        p.set_defaults(func=mod.main)

    args = parser.parse_args()
    args.func(args)
