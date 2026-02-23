from . import bootstrap, chat, doctor, info, init, repl, version


def register_subcommands(subparsers) -> None:
    for module in (info, doctor, init, bootstrap, chat, repl, version):
        module.register(subparsers)
