from .. import __version__


def register(subparsers) -> None:
    parser = subparsers.add_parser("version", help="Show version")
    parser.set_defaults(func=run)


def run(_args) -> int:
    print(__version__)
    return 0

