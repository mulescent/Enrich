from __future__  import print_function
from aligner import Aligner

a = Aligner()

x = "TCGAACTGAAAA"
y = "AACTGA"

trace = a.align(x, y)

print(x)
print(y)
for t in trace:
	print(t)
