#!/usr/bin/env python

"""Tests for `seqpandas` package."""

import pytest

from click.testing import CliRunner

from seqpandas import seqpandas
from seqpandas import cli


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    # The main command is now a group, so calling it without subcommands shows help
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0
    assert "SeqPandas CLI" in result.output
    assert "--help" in result.output
    # Check that subcommands are listed
    assert "read" in result.output
    assert "read-vcf" in result.output
