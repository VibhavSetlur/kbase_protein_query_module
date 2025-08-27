"""Command-line interface."""

import click


@click.command()
@click.version_option()
def main() -> None:
    """KBUtilLib."""


if __name__ == "__main__":
    main(prog_name="KBUtilLib")  # pragma: no cover
