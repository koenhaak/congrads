--- CONGRADS -------------------------------------------------------------------------- 

Developed at DCCN (Donders Centre for Cognitive Neuroimaging), Donders Institute for 
Brain, Cognition and Behaviour. Radboud University, Nijmegen, The Netherlands

Authors: KV Haak, AF Marquand, CF Beckmann

Congrads is a data-driven method for mapping the spatial organisation of connectivity
within a pre-defined region-of-interest (i.e. connectopic mapping) and fitting the
results with a spatial statisical model that summarises the map in terms of a limited
number of trend surface model coefficients. These coefficients can then be compared
across subjects or groups using standard statistical analysis methods or used as 
features in e.g. classification. The code has been tested using fslpython (FSL 5.0.10) 
but should also work with other Python distributions that include packages numpy, 
nibabel, scipy, and networkx.

If you use CONGRADS in your research, please quote the following journal reference:

Haak KV, Marquand AF, Beckmann CF (2018) Connectopic mapping with resting-state fMRI. 
NeuroImage 170:83-94. 

This package contains:

 - congrads      [wrapper]
 - conmap.py     [connectectopic mapping code]
 - trendsurf.py  [trend surface modeling code]
 - bayesreg.py   [Bayesian regression for TSM]

See the attached LICENSE file for the terms and conditions for use, reproduction and
distribution.

--- USAGE ----------------------------------------------------------------------------- 

Congrads can be run by executing the 'congrads' wrapper script, e.g.:

$ ./congrads

Usage:
./congrads -i <input> -r <roi> -m <mask> -o <out> [options]

Compulsory arguments:
	-i			input file name [4D image]
	-r			region-of-interest [3D binary image]
	-m			mask [3D binary image]
	-o			path to output directory

Optional arguments:
	-n <int>	number of connectopic maps [default=1]
	-s			save eta-squared matrix
	-z			normalise output maps to range [0-1]
	-p			project connectopic maps into mask
	-f <int>	fit spatial model of order <int>
	-F <int>	fit spatial model of order <int>
				to pre-estimated connectopic map(s)
				compulsory arguments: -i -r -o

Note that it is possible to use the -i argument multiple times to load multiple data-
sets (e.g.: $ ./congrads -i <func1> -i <func2> ...). If so, the data will be combined
by computing the similarity matrix S (see Haak et al. 2018) for each data-set, which
are then averaged to obtain the similarity matrix S for the combined data.

Also note that the -F <int> option requires the -i -r and -o arguments but not the -m
argument. Further, the -i option should be used only once and the input file should
be either a 4D (or 3D) image of connectopic maps (cmaps or pmaps, see below).

It is also possible to run the .py scripts without the wrapper (run without arguments
for usage). 

--- OUTPUT ----------------------------------------------------------------------------

When the wrapper script is run using only the compulsory arguments, the output folder
will contain the following file:

 - <roiname>.cmaps.nii.gz

When the wrapper script is run with the following optional arguments, the output folder 
will additionally contain these files:

 -s       :: <roiname>.eta2.mat                    :: similarity matrix S
 -p       :: <roiname>.pmaps.nii.gz                :: projection of map onto mask
 -f <int> :: <roiname>.cmaps.tsm.trendcoeff.txt    :: spatial model parameter estimates
             <roiname>.cmaps.tsm.trendcoeffvar.txt :: marginal variances
             <roiname>.cmaps.tsm.explainedvar.txt  :: explained variance by model
             <roiname>.cmaps.tsm.hyp.txt           :: hyperparameters
             <roiname>.cmaps.tsm.negloglik.txt     :: negative log marginal likelihood
             <roiname>.cmaps.tsm.rmse.txt          :: standardised mean squared error
             <roiname>.cmaps.tsm.yhat.txt          :: predictive mean
             <roiname>.cmaps.tsm.ys2.txt           :: predictive variance

--- NOTES -----------------------------------------------------------------------------
	
This package involves spatial modeling procedures that differ slightly from those 
described in Haak et al. (2018) as they do not include the Gaussian Process term in 
equation (3) of that paper, and the models are fit using Bayesian Linear Regression, 
where the model hyperparameters (controlling the noise and data variance) are set 
using an emperical Bayes approach. See Marquand, Haak & Beckmann (2017) Functional
corticostriatal connection topographies predict goal-directed behaviour in humans.
Nature Human Behaviour 1(8):0146 for details. 

For optimal results it is important that the region-of-interest is defined accurately 
and that the input data have been adequately pre-processed and denoised (using e.g. 
ICA-FIX or ICA-AROMA). Generally, some level of spatial smoothing is also advised. 
In addition, because connectopic mapping involves a multivariate data analysis, the 
quality of the results will improve with both increased temporal and spatial degrees 
of freedom. Thus, longer functional MRI scans and higher spatial resulution will 
typically improve the quality of the results (so long as the SNR stays the same).

Comparisons of the spatial model parameters across subjects requires the maps to be in
the same space (e.g. group-average or MNI). This can be achieved by either performing
connectopic mapping on spatially normalised functional data or by running connectopic 
mapping on data in native space, normalising the ensuing maps, and then running the 
spatial modeling on the normalised maps (using the -F option). In case of the latter, 
the spatial normalisation should involve nearest-neighbour interpolation.

