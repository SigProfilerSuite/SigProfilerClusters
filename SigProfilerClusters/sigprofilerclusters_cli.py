#!/usr/bin/env python3

import sys
from SigProfilerClusters.controllers import cli_controller


def main_function():
    commands = {
        "analysis": "Analyze clustered mutations"
    }

    if len(sys.argv) < 2 or sys.argv[1].lower() not in commands:
        print_usage(commands)
        return

    command = sys.argv[1].lower()
    args = sys.argv[2:]

    controller = cli_controller.CliController()

    if command == "analysis":
        controller.dispatch_sigProfilerClusters(args)


def print_usage(commands):
    """Prints the usage message."""
    print("Usage: SigProfilerClusters <command> [<args>]\n")
    print("Commands:")
    for cmd, desc in commands.items():
        print(f"  {cmd}: {desc}")


if __name__ == "__main__":
    main_function()
