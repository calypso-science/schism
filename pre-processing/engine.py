from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
import sys
import os

import click
import yaml


from wrapper import cycle_sim

@click.group()
def main():
    pass


@main.command()
@click.argument("config", envvar="CONFIG", default=None)
@click.argument("cycle", envvar="CYCLE", default=None)
@click.argument("cycle_len", envvar="CYCLE_LEN", default=1)
@click.argument("cycle_len_unit", envvar="CYCLE_LEN_UNIT", default='months')
@click.option("--cycle_format", "-f", default="%Y%m%d", help="cycle format")
def main(config, cycle,cycle_len,cycle_len_unit,cycle_format):
    """Run model
    Usage: python cycle_sim.py config.yml preprocessing 20010101T00
    Args:
        config(str): yaml config file
        method(str):
        cycle(str):  cycle to run (e.g. 200101T00)
    """

    start_time = datetime.strptime(cycle,cycle_format)
    end_time = start_time+relativedelta(**{cycle_len_unit:cycle_len})

    cycle_sim(action=config, **{'timing':{'time start':start_time,'time end':end_time}})

if __name__ == "__main__":
    main()