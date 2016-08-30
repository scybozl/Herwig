rm fifo.hepmc
mkfifo fifo.hepmc
Herwig run tT_homemade.run -N 1000 -x setupfifo.in &
rivet -a ATLAS_2014_I1304688_custom fifo.hepmc

## Herwig at LO without shower, or hadronization
Herwig read tT_homemade.in
Herwig run tT_homemade.run -N 100000 -x setupfile.in
rivet -a ATLAS_2014_I1304688_custom tT_homemade-setupfile.in.hepmc
cp Rivet.yoda tT-Analyses/ATLAS_I1304688_LO.yoda

rm tT_homemade-setupfile.in.hepmc

## Herwig at LO with shower, no hadronization
Herwig run tT_homemade.run -N 100000 -x setupfile_shower.in
rivet -a ATLAS_2014_I1304688_custom tT_homemade-setupfile_shower.in.hepmc
cp Rivet.yoda tT-Analyses/ATLAS_I1304688_LO_shower.yoda

rm tT_homemade-setupfile_shower.in.hepmc

## Herwig at LO with shower and hadronization

#Herwig run tT_homemade.run -N 100000 -x setupfile_shower+had.in
#rivet -a ATLAS_2014_I1304688_custom tT_homemade-setupfile_shower+had.in.hepmc
#cp Rivet.yoda tT-Analyses/ATLAS_I1304688_LO_shower_had.yoda



## Herwig at NLO without shower, or hadronization

Herwig read tT_homemadeNLO.in
Herwig run tT_homemadeNLO.run -N 100000 -x setupfile.in
rivet -a ATLAS_2014_I1304688_custom tT_homemadeNLO-setupfile.in.hepmc
cp Rivet.yoda tT-Analyses/ATLAS_I1304688_NLO.yoda

rm tT_homemadeNLO-setupfile.in.hepmc

## Herwig at NLO with shower, no hadronization

Herwig run tT_homemadeNLO.run -N 100000 -x setupfile_shower.in
rivet -a ATLAS_2014_I1304688_custom tT_homemadeNLO-setupfile_shower.in.hepmc
cp Rivet.yoda tT-Analyses/ATLAS_I1304688_NLO_shower.yoda

rm tT_homemadeNLO-setupfile_shower.in.hepmc

## Herwig at NLO with shower and hadronization

#Herwig run tT_homemadeNLO.run -N 100000 -x setupfile_shower+had.in
#rivet -a ATLAS_2014_I1304688_custom tT_homemadeNLO-setupfile_shower+had.in.hepmc
#cp Rivet.yoda tT-Analyses/ATLAS_I1304688_NLO_shower_had.yoda


