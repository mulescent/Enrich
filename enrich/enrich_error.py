from __future__ import print_function


class EnrichError(Exception):
    def __init__(self, value, name):
        self.value = value
        self.name = name
    def __str__(self):
        return '%s [Enrich name: "%s"]' % (repr(self.value), self.name)



class EnrichMessage(object):
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)
