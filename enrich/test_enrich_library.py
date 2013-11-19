from __future__ import print_function
import enrich_library

def test_set_wt(seq):
    lib = enrich_library.EnrichLibrary()
    print("Testing set_wt\nSEQ: %s" % seq)
    lib.set_wt(seq, coding=True)
    print("DNA:", lib.wt_dna)
    print("PRO:", lib.wt_protein)
    lib.set_wt(seq, coding=False)
    print("DNA:", lib.wt_dna)
    print("PRO:", lib.wt_protein)

def test_call_variants(wt, variant):
    lib = enrich_library.EnrichLibrary()
    print("Testing call_variants\nWT:  %s\nVAR: %s" % (wt, variant))
    lib.set_wt(wt, coding=True)
    result = lib.call_variants(variant)
    for x in result:
        print("\t%s" % x)

if __name__ == "__main__":
    test_set_wt("AAATTTCCCGGG")
    test_set_wt("AAATTTCCCGG")
    test_set_wt("AAATTCCCGGG")
    test_set_wt("AAATTTCCCGGN")
    test_set_wt("AAATTNCCCGGG")

    test_call_variants("AAATTTCCCGGG", "AAATTTCCCGGG")
    test_call_variants("AAATTTCCCGGG", "CAATTTCCCGGG")
