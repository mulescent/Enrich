.. include:: global.rst

:py:mod:`~fqread` --- Manipulation of FASTQ records
===================================================

.. py:module:: fqread
	:synopsis: Manipulation of FASTQ records.

The :py:mod:`~fqread` module contains the :py:class:`~fqread.FQRead` class for storing and manipulating FASTQ_ records, and associated utility functions for reading these data from standard files.


:py:class:`~fqread.FQRead` class
--------------------------------
.. autoclass:: FQRead
    :members:
    :special-members:

Generator functions
-------------------
The :py:mod:`~fqread` module provides two generators (functions that return iterators) for reading records from FASTQ_ files. Input files are read in chunks to improve performance by minimizing disk accesses.

.. autofunction:: read_fastq

.. autofunction:: read_fastq_multi

Miscellaneous functions
-----------------------
.. autofunction:: check_fastq

.. autofunction:: fastq_filter_chastity
