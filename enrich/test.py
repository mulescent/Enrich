from __future__ import print_function
import selection
import json


config = json.load(open("/Users/afrubin/Code/Enrich/enrich/example.json"))
test = selection.Selection(config)
test.enable_logging(open("/Users/afrubin/Documents/Enrich2_example/log", "w"))
test.get_counts()
test.calc_all()
#test.calc_barcode_variation()
test.write_enrichments("/Users/afrubin/Documents/Enrich2_example/")
