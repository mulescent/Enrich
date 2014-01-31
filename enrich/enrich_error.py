from __future__ import print_function


class EnrichError(Exception):
    """
    Exception class for errors in the Enrich2 project.
    """
    def __init__(self, value, name):
        self.value = value
        self.name = name
    def __str__(self):
        return '%s [Enrich name: "%s"]' % (repr(self.value), self.name)



class EnrichMessage(object):
    """
    Class for warning messages in the Enrich2 project.

    .. warning:: Not currently used.
    """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)
