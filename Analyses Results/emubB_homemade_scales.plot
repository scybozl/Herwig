# BEGIN PLOT /ATLAS_sherpa/Eta1
LegendXPos=0.65
Legend=1
XLabel=$\eta_{\mathrm{b}_1}$
Title=$W^+W^-b \bar{b}$ cross-section vs. first jet rapidity
YLabel=$\mathrm{d}\sigma/{d\eta}$ [pb]
# END PLOT

.*emubB_homemade_scales.yoda/ATLAS_sherpa/Eta1::ErrorBars=0
.*WWbBSherpa2DScatter/ATLAS_sherpa/Eta1::ErrorBars=1

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/Eta1
LineStyle=dashed
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
ErrorBars=0
Title=Herwig LO
# END HISTOGRAM

# BEGIN HISTO1D WWbBSherpa2DScatter.yoda/ATLAS_sherpa/Eta1
Title=Sherpa LO
ErrorBars=1
# END HISTO1D

# BEGIN PLOT /ATLAS_sherpa/HT
LegendXPos=0.65
Legend=1
XLabel=$H_T$
Title=$W^+W^-b \bar{b}$ cross-section vs. final state pT sum
YLabel=$\mathrm{d}\sigma/{dH_T}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/HT
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/HT
Title=Sherpa LO
ErrorBars=1
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/Mlb
LegendXPos=0.65
Legend=1
XLabel=$m_{lb}$
Title=$W^+W^-b \bar{b}$ cross-section vs. invariant lepton + b system mass
YLabel=$\mathrm{d}\sigma/{dm_{lb}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/Mlb
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/Mlb
Title=Sherpa LO
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/Phiemu
LegendXPos=0.65
Legend=1
XLabel=$\phi_{e^+mu^-}$
Title=$W^+W^-b \bar{b}$ cross-section vs. leptons azimuthal angle
YLabel=$\mathrm{d}\sigma/{d\phi_{e^+mu^-}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/Phiemu
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/Phiemu
Title=Sherpa LO
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/Pt1
LegendXPos=0.65
Legend=1
XLabel=$p_{\mathrm{T, b_1}}$
Title=$W^+W^-b \bar{b}$ cross-section vs. first jet $p_\mathrm{T}$
YLabel=$\mathrm{d}\sigma/{dp_{\mathrm{T, b_1}}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/Pt1
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/Pt1
Title=Sherpa LO
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/PtMiss
LegendXPos=0.65
Legend=1
XLabel=$p_{\mathrm{T, miss}}$
Title=$W^+W^-b \bar{b}$ cross-section vs. missing transverse momentum
YLabel=$\mathrm{d}\sigma/{dp_{\mathrm{T, miss}}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/PtMiss
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/PtMiss
Title=Sherpa LO
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/deltaRb
LegendXPos=0.65
Legend=1
XLabel=$\Delta R_{b\bar{b}}$
Title=$W^+W^-b \bar{b}$ cross-section vs. b-jets angular distance
YLabel=$\mathrm{d}\sigma/{d\Delta R_{b,\bar{b}}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/deltaRb
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/deltaRb
Title=Sherpa LO
# END HISTOGRAM

# BEGIN PLOT /ATLAS_sherpa/deltaRl
LegendXPos=0.65
Legend=1
XLabel=$\Delta R_{l^+l^-}$
Title=$W^+W^-b \bar{b}$ cross-section vs. leptons angular distance
YLabel=$\mathrm{d}\sigma/{d\Delta R_{l^+l^-}}$ [pb]
# END PLOT

# BEGIN HISTOGRAM emubB_homemade_scales.yoda/ATLAS_sherpa/deltaRl
ErrorBands=1
ErrorBandColor=Gray
ErrorBandOpacity=0.1
Title=<Herwig LO>
# END HISTOGRAM

# BEGIN HISTOGRAM WWbBSherpa2DScatter.yoda/ATLAS_sherpa/deltaRl
Title=Sherpa LO
# END HISTOGRAM
