A Data Analysis Algorithm for Automated Relative Quantification of Stable Isotope-Labeled Metabolites
http://sdrv.ms/17Z61DX


##Introduction:
Quantification by liquid chromatography-mass spectrometry (LC-MS) requires special consideration due to the potential for run-to-run variability of retention times and matrix effects such as ionization suppression during the commonly used electrospray ionization. Several different strategies have been employed for quantitative LCMS of metabolomic samples. Arguably, the most precise method involves the addition of an isotopically-labeled internal standard for every compound of interest (e.g. 2H, 13C, or 15N-labeled). This approach is expensive, but warranted in certain targeted studies. Another strategy for relative quantification, as opposed to absolute quantification, is chemical labeling, which has proven to be useful for quantification in proteomics (e.g., isotope-coded affinity tags). Though relative quantification by labeling has seen limited use for metabolomics due to the lack of a single functional group present in all metabolites, there has been none-the-less a number of recent reports of effective labeling schemes (Guo, 2007; Huang, 2008; Lamos, 2007; Shortreed, 2006; Yang, 2007). Labeled metabolites coelute from the chromatographic separation and appear in the mass spectrum as pairs of peaks with a characteristic mass difference. The peak intensity ratio for each pair yields the relative concentration. Such labeling strategies have a number of advantages including: improved quantitative precision, increased ability for molecular identification, and enhanced detection sensitivity. A major limitation of this strategy is lack of software tools for global identification of labeled compounds and calculation of the peak intensity ratios. Here we report a new software tool for automated identification of isotopically labeled metabolites and quantification of their relative concentration.

##Methods:
The software described here is compatible with XCMS, a widely used and freely available software program for processing LC-MS data for metabolite profiling. Both the new software and XCMS are written in the ‘R’ language. Data processing begins by applying a number of standard XCMS functions as follows: define the working directory, find peaks, group peaks, correct retention time, and refine peak groupings. The new code then scans the output for pairs of peaks with identical retention time and selected isotopic mass difference. Data from peak pairs (retention time, observed m/z, neutral unlabelled mass, peak intensity, peak area, peak ratio) is written to a second file. Extracted ion chromatograms and m/z plots are automatically produced for each pair.

##Preliminary Results:
The identification of isotope-labeled metabolites and calculation of relative abundances using the newly developed software is illustrated through analysis of samples of fatty acids that were obtained by saponifying lipids extracted from egg yolks. Fatty acids were reacted with either the heavy(D9)- or the light(D0)-isotopic form of cholamine, which contains a fixed-charged quaternary ammonium group. The heavy and light labeled samples were subsequently mixed in a 1:1 ratio, separated by reverse-phase HPLC and analyzed in positive-ion mode ESI-TOF-MS. A mixture of all samples is used as the control, which guarantees that all compounds are represented, even if they appear in only one of the samples. Our strategy for identifying labeled metabolites is to analyze a 1:1 mixture of heavy and light labeled control samples. Labeled metabolite pairs coelute, have equal intensity and display a 9Da mass difference. Differences in retention time, intensity or mass shift can be used to eliminate false positives. The data from all runs (control vs. control and sample vs. control) are processed using XCMS software. Nearly 1000 features were identified. The new software then identifies candidate peak pairs (retention times +/-5s and mass difference 9.0565+/-0.0050Da in this case, although other values can be specified). Eighty-five peak pairs satisfied these criteria. Extracted ion chromatograms and mass spectra are automatically generated for each pair to facilitate rapid visual scanning of the data. The 85 pairs were further refined to 46 by selecting those with an intensity ratio of 1.00+/0.15 from the 1:1 control mixtures. The average ratio for the 46 pairs from the 1:1 control samples was 1.00. The average of the individual standard deviations for each ratio was 0.05. The relative abundances for the 46identified metabolites in the sample versus control were then calculated.

##Conclusion:
The new software algorithm automatically identifies pairs of isotopically labeled metabolites, calculates the intensity ratio and generates extracted-ion chromatograms for each pair. The total processing time for all three sample groups (24 total files, 12 Mbyte each) analyzed for this work was under 10 minutes. In the future, we plan to enhance the functionality of this program by enabling it to collate data from all sample groups into a single spreadsheet and remove satellite peak pairs that emanate from natural isotopic abundances of 13C. Furthermore, we plan to evaluate the performance of the software algorithm using data from other instrument systems (e.g. GC-MS). The algorithm was designed with compatibility to other instruments in mind. However, this functionality remains to be evaluated.


##Project supervised by:
- Dr. Michael Shortreed
- Dr. Brian Frey
- Dr. Lloyd Smith


##Reference:
- Stable Isotope Labeling
Guo, K., Ji, C. and Li, L. (2007) Stable-isotope dimethylation labeling combined with LC-ESI MS for quantification of amine-containing metabolites in biological samples. Anal. Chem., 79, 8631-8638.
Huang, X. and Regnier, F.E. (2008) Differential metabolomics using stable isotope labeling and two-dimensional gas chromatography with time-of-flight mass spectrometry. Anal. Chem., 80, 107-114.
Lamos, S.M., Shortreed, M.R., Frey, B.L., Belshaw, P.J. and Smith, L.M. (2007) Relative quantification of carboxylic acid metabolites by liquid chromatography-mass spectrometry using isotopic variants of cholamine. Anal. Chem., 79, 5143-5149.
Shortreed, M.R., Lamos, S.M., Frey, B.L., Phillips, M.F., Patel, M., Belshaw, P.J. and Smith, L.M. (2006) Ionizable isotopic labeling reagent for relative quantification of amine metabolites by mass spectrometry. Anal. Chem., 78, 6398-6403.
Yang, W.C., Adamec, J. and Regnier, F.E. (2007) Enhancement of the LC/MS analysis of fatty acids through derivatization and stable isotope coding. Anal. Chem., 79, 5150-5157.
- XCMS
Smith, C.A., Want, E.J., O'Maille, G., Abagyan, R. and Siuzdak, G. (2006) XCMS: processing mass spectrometry data for metabolite profiling using nonlinear peak alignment, matching, and identification. Anal. Chem., 78, 779-787.


##Prerequisite knowledge:
-- The R Project for Statistical Computing - visit http://www.r-project.org/
-- XCMS - visit http://metlin.scripps.edu/xcms/
-- mzXML file format - visit http://sashimi.sourceforge.net/software_glossolalia.html


##To run the tool:
R install
Packages -> set CRAN -> USA (MI)
Select Repositories -> CRAN -> Everything except for Omegahat
Install packages -> xcms and all dependent pacakages

library(xcms)
source("findPairs.R")
source("calcRatio.R")
setAll     <- xcmsSet(snthresh=10)
setInfo    <- classInfo(mode = 1, setAll)
setAllg    <- group(setAll)
pairs      <- findPairs(setAll, mz_shift=2.0067, mz_tol=0.005, rt_tol=30, light_label=28.0312)
features   <- calcRatio(setAllg, pairs=pairs$pairs, plot=TRUE, plotSpectrum=F, dir='test',
                        filter=F, filter_tol=0.15, mz_width=0.05, target_intensity=2/3, rt_tolerance=100000, mode = 1)
write.csv(pairs$pairs,       'pairs.csv')
write.csv(pairs$peaklist,    'peaklist.csv')
write.csv(features$ratio,    'ratio.csv')
write.csv(features$features, 'result.csv')
write.csv(features$final,    'final.csv')
write.csv(features$stdev,    'stdev.csv')


##Code:
The code was written back in 2008.


##Collaborators:
Paul Benton, Colin Smith, Ralf Tautenhahn from the Scripps Center for Metabolomics and Mass Spectrometry
