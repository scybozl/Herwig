import yoda
import math

aos = [ ao for ao in yoda.read("WWbBPOWHEG.yoda", asdict=False) ]
result = [ ao.mkScatter() for ao in yoda.read("WWbBPOWHEG.yoda", asdict=False)]

for i in range(0, len(aos)-3):
	for p in range(0, aos[i].numBins):
		result[i].points[p].y = aos[i].bins[p].sumW
		result[i].points[p].yErrs = [math.sqrt(aos[i].bins[p].sumW2), math.sqrt(aos[i].bins[p].sumW2)]

yoda.write(result, "WWbBPOWHEGScatter.yoda")
