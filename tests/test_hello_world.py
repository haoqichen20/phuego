# /repo/tests/test_hello_world.py
from click.testing import CliRunner
from phuego.cli_app import greet

def test_greet():
    runner = CliRunner()
    result = runner.invoke(greet, ['--name', 'World'])
    assert result.exit_code == 0
    assert 'Hello, World!' in result.output
