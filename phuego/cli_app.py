# /repo/phuego/cli_app.py
import click

@click.command()
@click.option('--greeting', default='Hello', help='Greeting phrase.')
@click.option('--name', default='World', help='Name to greet.')
def greet(greeting, name):
    """Simple program that greets NAME for a total of COUNT times."""
    click.echo(f"{greeting}, {name}!")

if __name__ == '__main__':
    greet()
