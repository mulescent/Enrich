def is_experiment(config):
    """
    Return ``True`` if *config* is for an :py:class:`~experiment.Experiment`.
    """
    if 'conditions' in config.keys():
        return True
    else:
        return False


def is_selection(config):
    """
    Return ``True`` if *config* is for a :py:class:`~selection.Selection`.
    """
    if 'libraries' in config.keys():
        return True
    else:
        return False


def is_seqlib(config):
    """
    Return ``True`` if *config* is for a :py:class:`~seqlib.seqlib.SeqLib`.
    """
    if 'fastq' in config.keys():
        return True
    else:
        return False


def seqlib_type(config):
    """
    Return the class name of the appropriate :py:class:`~seqlib.seqlib.SeqLib` 
    type for *config*. Returns ``None`` if there is no appropriate 
    :py:class:`~seqlib.seqlib.SeqLib` type.
    """
    if 'barcodes' in config:
        if 'wild type' in config:
            return "BarcodeVariantSeqLib"
        else:
            return "BarcodeSeqLib"
    elif 'overlap' in config and 'wild type' in config:
        return "OverlapSeqLib"
    elif 'wild type' in config:
        return "BasicSeqLib"
    else:
        return None

