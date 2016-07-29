import yoda

aosref=yoda.read("emubB_homemade.yoda")
aosHalf=yoda.read("emubB_homemade0.5.yoda")
aosDouble=yoda.read("emubB_homemade2.0.yoda")

aos = [ ao for ao in yoda.read("emubB_homemade.yoda", asdict=False) ]
aosHalf = [ao for ao in yoda.read("emubB_homemade0.5.yoda", asdict=False)]
aosDouble = [ao for ao in yoda.read("emubB_homemade2.0.yoda", asdict=False)]

result = [ ao.mkScatter() for ao in yoda.read("emubB_homemade.yoda", asdict=False)]

for i in range(0, len(aos)-3):
	for p in range(0, aos[i].numBins):
		Half=aosHalf[i].bins[p].sumW
		Double=aosDouble[i].bins[p].sumW

		errP = max(Half, Double) - aos[i].bins[p].sumW
		errM = aos[i].bins[p].sumW - min(Half, Double)
		result[i].points[p].y = aos[i].bins[p].sumW
		result[i].points[p].yErrs = [max(errM, 0), max(errP, 0)]

yoda.write(result, "emubB_homemade_scales.yoda")
