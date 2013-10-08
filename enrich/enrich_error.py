from __future__ import print_function

class EnrichError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)



class EnrichFileError(EnrichError):
	def __init__(self, fname):
		self.value = 'failed to process "%s"' % fname



class EnrichMessage(object):
	def __init__(self, message):
		self.message = message
	def __str__(self):
		return repr(self.message)
