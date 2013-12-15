from __future__ import print_function
import selection
import json


config = json.load(open("/Users/afrubin/Code/Vanessa/config.json"))
test = selection.Selection(config)
test.get_counts()
test.calc_all()
test.calc_barcode_variation()
test.write_enrichments("/Users/afrubin/Code/Vanessa/")
test.barcode_map.write_variants("/Users/afrubin/Code/Vanessa/bcv.txt")