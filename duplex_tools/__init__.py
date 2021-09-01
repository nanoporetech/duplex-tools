"""Duplex Sequencing Tools package."""
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import logging

from duplex_tools import \
    assess_split_on_adapter, filter_pairs, pairs_from_summary, \
    split_on_adapter  # noqa: F401

modules = [
    "split_on_adapter", "assess_split_on_adapter",
    "pairs_from_summary", "filter_pairs"]

__version__ = '0.2.0'


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
    parser.add_argument(
        '--version', action='version',
        version='%(prog)s {}'.format(__version__))
    args = parser.parse_args()

    logging.basicConfig(
        format='[%(asctime)s - %(name)s] %(message)s',
        datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)
    args.func(args)


def _log_level():
    """Parser to set logging level and acquire software version/commit."""
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter, add_help=False)
    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument(
        '--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument(
        '--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only).')
    return parser


def get_named_logger(name):
    """Create a logger with a name."""
    logger = logging.getLogger('{}.{}'.format(__package__, name))
    logger.name = name
    return logger
