"""Tests for the CLI."""

from click.testing import CliRunner
from plab2.cli import compile
from pathlib import Path

THIS_DIR = Path(__file__).parent
THIS_PATH = THIS_DIR.joinpath("data")
PPI_FILE = THIS_PATH.joinpath('ppis.csv')
NODES_FILE = THIS_PATH.joinpath('nodes.tsv')
EDGES_FILE = THIS_PATH.joinpath('edges.tsv')


class TestCli:
    def test_compile(self):
        """Tests for subcommand compile"""
        runner = CliRunner()
        result = runner.invoke(compile, [str(PPI_FILE), str(NODES_FILE), str(EDGES_FILE), "--enrich"])
        assert result.exit_code == 0
        assert NODES_FILE.is_file()
        assert EDGES_FILE.is_file()
        NODES_FILE.unlink()
        EDGES_FILE.unlink()
        assert not NODES_FILE.is_file()
        assert not EDGES_FILE.is_file()

        help_result = runner.invoke(compile, ["--help"])
        assert "Enrich the graph with RNA and DNA molecules." in help_result.output







