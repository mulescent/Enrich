from __future__ import print_function
import experiment
import json


config = json.load(open("/Users/afrubin/Code/Enrich/enrich/example.json"))
test = experiment.Experiment(config)
test.enable_logging(open("/Users/afrubin/Documents/Enrich2_example/log", "w"))
test.calc_selection_scores()
test.conditions['normal'][0].save_data("/Users/afrubin/Documents/Enrich2_example/new")
#test.calc_barcode_variation()
