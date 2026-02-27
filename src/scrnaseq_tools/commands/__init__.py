from . import bootstrap, chat, doctor, info, init, repl, scan, summary, version

# Lazy-load analysis commands (only register if importable)
_ANALYSIS_COMMANDS = []
try:
    from . import qc, pipeline_cmd, integrate, profile_cmd
    _ANALYSIS_COMMANDS = [qc, pipeline_cmd, integrate, profile_cmd]
except ImportError:
    pass


def register_subcommands(subparsers) -> None:
    for module in (info, doctor, init, bootstrap, chat, repl, scan, summary, version):
        module.register(subparsers)
    for module in _ANALYSIS_COMMANDS:
        module.register(subparsers)
