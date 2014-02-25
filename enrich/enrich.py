from __future__ import print_function
import argparse
import json
import config_check
from enrich_error import EnrichError
from experiment import Experiment
from selection import Selection
from seqlib.basic import BasicSeqLib
from seqlib.barcodevariant import BarcodeVariantSeqLib
from seqlib.barcode import BarcodeSeqLib
from seqlib.overlap import OverlapSeqLib


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="JSON configuration file")
    parser.add_argument("-l", "--log", metavar="file", help="path to new log file")
    parser.add_argument("--no-plots", help="don't make plots", dest="plots", action="store_false", default=True)
    args = parser.parse_args()

    try:
        config = json.load(open(args.config, "U"))
    except IOError:
        raise EnrichError('Failed to open "%s"' % args.config, "enrich.py")
    except ValueError:
        raise EnrichError("Improperly formatted .json file", "enrich.py")

    if config_check.is_experiment(config):
        obj = Experiment(config)
    elif config_check.is_selection(config):
        obj = Selection(config)
    elif config_check.is_seqlib(config):
        obj = globals()[config_check.seqlib_type(config)](config)
    else:
        raise EnrichError("Unrecognized .json config", "enrich.py")

    if args.log:
        obj.enable_logging(open(args.log, "w"))

    obj.calculate()
    obj.write_data()
    if args.plots:
        obj.make_plots()
    