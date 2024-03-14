# -*- coding: utf-8 -*-

import os
import ast
import click
import pandas as pd

    
def load_test_example():
    test_path = os.path.join(os.path.dirname(__file__), 'data/EGF_vs_Untreated_@min_15_63_240.txt')
    df = pd.read_csv(test_path, sep='\t', header=None)
    return test_path, df


# Parsing user input literally for list input.
class PythonLiteralOption(click.Option):

    def type_cast_value(self, ctx, value):
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)


# Group click options for organized help message display.
class GroupedOptions(click.Command):
    def format_options(self, ctx, formatter):
        options = sorted(self.get_params(ctx), key=lambda option: option.opts[0])
        
        # Retrieve grouped_options from the callback (the command function)
        grouped_options = getattr(self.callback, 'grouped_options', {})
        
        other_options = [opt for opt in options if opt.name not in sum(grouped_options.values(), [])]
        
        for group_name, opt_names in grouped_options.items():
            opts = [opt for opt in options if opt.name in opt_names]
            if not opts:
                continue
            
            with formatter.section(group_name):
                formatter.write_dl([o.get_help_record(ctx) for o in opts])

        if other_options:
            with formatter.section('Other Options'):
                formatter.write_dl([o.get_help_record(ctx) for o in other_options])