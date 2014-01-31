.. include:: global.rst

Sequencing library modules
==========================

Data for each FASTQ_ file (or pair of FASTQ_ files for overlapping paired-end data) is read into its own SeqLib-family object.

.. toctree::
    :maxdepth: 4

    seqlib
    variant
    basic
    overlap
    barcode
    barcodevariant

Inheritance diagram for :py:class:`~seqlib.seqlib.SeqLib` derived classes
-------------------------------------------------------------------------
.. inheritance-diagram:: seqlib.barcodevariant.BarcodeVariantSeqLib seqlib.overlap.OverlapSeqLib seqlib.basic.BasicSeqLib
	:parts: 1
