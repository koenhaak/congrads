# usage: conmap.py [-h] -i INFILES -r ROIFILE -m MASKFILE -o OUTDIR
#                  [--nmaps NMAPS] [--save_eta2 SAVE_ETA2] [--norm NORM_FLAG]
#
# Developed at DCCN (Donders Centre for Cognitive Neuroimaging), Donders Institute
# for Brain, Cognition and Behaviour. Radboud University, Nijmegen, The Netherlands
#
# Authors: KV Haak, AF Marquand, CF Beckmann.
#
# If you use this code in your research, please quote the following journal reference:
#
# Haak KV, Marquand AF, Beckmann CF (2018) Connectopic mapping with resting-state fMRI. 
# NeuroImage 170:83-94. 

import numpy as np

def pca(X):

    from scipy.linalg import svd

    # Center X by subtracting off column means
    X -= np.mean(X,0)

    # The principal components are the eigenvectors of S = X'*X./(n-1), but computed using SVD
    [U,sigma,V] = svd(X,full_matrices=False)

    # Project X onto the principal component axes
    Y = U*sigma

    # Convert the singular values to eigenvalues 
    sigma /= np.sqrt(X.shape[0]-1)
    evals = np.square(sigma)
    
    return V, Y, evals

def corr(X,Y):
    
    Y = Y.T
    X = X.T

    R = np.zeros((X.shape[0],Y.shape[0]))
    for i in range(0,R.shape[1]):
        y = Y[i,:]
        Xm = np.reshape(np.mean(X,axis=1),(X.shape[0],1))
        ym = np.mean(y)
        r_num = np.sum((X-Xm)*(y-ym),axis=1)
        r_den = np.sqrt(np.sum((X-Xm)**2,axis=1)*np.sum((y-ym)**2))
        R[:,i] = r_num / r_den
    
    return R

def eta2(X):
    
    S = np.zeros((X.shape[0],X.shape[0]))
    for i in range(0,X.shape[0]):
        for j in range(i,X.shape[0]):
            mi = np.mean([X[i,:],X[j,:]],0) 
            mm = np.mean(mi)
            ssw = np.sum(np.square(X[i,:]-mi) + np.square(X[j,:]-mi))
            sst = np.sum(np.square(X[i,:]-mm) + np.square(X[j,:]-mm))
            S[i,j] = 1-ssw/sst
    
    S += S.T 
    S -= np.eye(S.shape[0])
    
    return S

def norm(X):
    
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import squareform

    return squareform(pdist(X))

def adjacency(X):

    from networkx import is_connected
    from networkx import from_numpy_matrix
    
    emin = 0
    emax = np.max(X)
    tol = 0.0001
    maxiter = 1000
    cntr = 0
    
    done = False
    while not done:        
        e = (emin + emax) / 2
        A = (X < e) - np.eye(X.shape[0])
        G = from_numpy_matrix(A)
        if is_connected(G):
            emax = e
            if (emax - emin) < tol:
                done = True
        else:
            emin = e       
        cntr += 1
        if cntr == maxiter:
            done = True      
    
    return A

# Main routine 	
def main(infiles,roifile,maskfile,outdir,nmaps,save_eta2=False,norm_flag=False,proj_flag=False):

	import nibabel as nib
	import sys
	import errno
	np.seterr(invalid='ignore')

	out_base_name = roifile.split('/')[-1].split('.nii')[0]

	# Load the roi
	try:
		print('Loading roi from: ' + roifile)
		roiImg = nib.load(roifile)
		roi = roiImg.get_data()
	except:
		sys.exit('Cannot open ' + roifile | '\nExiting.')
	if len(roi.shape) != 3:
		sys.exit(roifile + ' is not a 3D image\nExiting.')

	# Store the dimensions of the roi data for later use
	roidims = roi.shape
	nVoxels = np.prod(roidims)

	# Reshape roi into a vector of size nVoxels
	roi = np.reshape(roi,(nVoxels))	

	# Find the indices inside roi
	roiIndices = np.where(roi>0)

	# Load the mask
	try:
		print('Loading mask from: ' + maskfile)
		maskImg = nib.load(maskfile)
		mask = maskImg.get_data()
	except:
		sys.exit('Cannot open ' + maskfile | '\nExiting.')
	if len(mask.shape) != 3:
		sys.exit(maskfile + ' is not a 3D image\nExiting.')
	
	# Reshape the mask into a vector of size nVoxels
	mask = np.reshape(mask,(nVoxels))

	# Find the indices outside roi but inside mask
	maskIndices = np.where((roi==0) & (mask>0)) 
	
	# Initialise similarity matrix
	S = np.zeros([np.sum(roi>0),np.sum(roi>0)])

	# Loop over infiles and create average similarity matrix
	for infile in infiles:

		print('Processing ' + infile)

		# Load functional data
		try:
			dataImg = nib.load(infile)
			data = dataImg.get_data()
		except:
			sys.exit('Cannot open ' + infile | '\nExiting.')
		if len(data.shape) != 4:
			sys.exit(infile + ' is not a 4D image\nExiting.')

		# Assert absence of nans and infs
		if np.any(np.isnan(data)) or np.any(np.isinf(data)):
			sys.exit('Data contains invalid values.\nExiting.')

		# Reshape and standardise
		nFrames = data.shape[3]
		data = np.reshape(data,(nVoxels,nFrames))
		data -= np.tile(np.mean(data,1),(nFrames,1)).T
		data /= np.tile(np.std(data,1),(nFrames,1)).T

		# Gather data inside roi
		A = data[roiIndices,:][0]

		# If the roi contains invalid data it must be due to a division by 0 (stdev)
		# since the data themselves do not contain nans or infs. If so, we terminate 
		# the program and the user should define a roi covering functional data.
		if np.any(np.isnan(A)) or np.any(np.isinf(A)):
			sys.exit('ROI includes voxels without variance.\nExiting.')

		# Gather data outside roi
		B = data[maskIndices,:][0]

		# Transpose so that the data are stored in a time x space matrix
		A = A.T
		B = B.T

		# A division by 0 (stdev) can also lead to nans and infs in the mask data.
		# In this case we can simply throw a warning and ignore all voxels without
		# variance.
		keepB = ~np.isnan(B).any(axis=0) & ~np.isinf(B).any(axis=0)
		if np.any(np.isnan(B)) or np.any(np.isinf(B)):
			print('WARNING: Mask includes voxels without variance.')
	
		del data  

		# Get voxel-wise connectivity fingerprints 
		print('Computing voxel-wise connectivity fingerprints...')
		[evecs,Bhat,evals] = pca(B[:,keepB])
		R = corr(A,Bhat)

		# Construct similarity matrix of connectivity fingerprints
		print('Computing similarity matrix...')
		S += eta2(R)

	if len(infiles) > 1:
		print('Creating average similarity matrix...')
		S /= len(infiles)

	# If requested, save the similarity matrix as a matlab .mat file
	if save_eta2:
		import scipy.io
		scipy.io.savemat(outdir + "/" +  out_base_name + ".eta2", dict(S=S))

	# Compute the graph Laplacian
	print('Computing the graph Laplacian...')
	dist = norm(S)**2
	W = np.multiply(adjacency(dist),S)
	D = np.diag(np.sum(W,0))
	L = np.subtract(D,W)
	
	# Solve generalised eigenvalue problem Ly = lDy
	print('Computing the dominant ' + str(nmaps) + ' connectopic maps...')
	from scipy.linalg import eigh
	l,y = eigh(L,D,eigvals=(0,nmaps))

	# The eigenvectors have an intrinsic sign indeterminacy, which is inconvenient
	# for spatial statistical modeling. We deal with this by flipping the sign of 
	# the eigenvectors if they correlate negatively with the reference vector 
	# defined below. 
	x0,y0,z0 = np.floor(roidims[0]/2),0,0
	X,Y,Z = np.ogrid[0:roidims[0],0:roidims[1],0:roidims[2]]
	ref = np.sqrt((X-x0)**2+(Y-y0)**2+(Z-z0)**2)
	ref = np.reshape(ref,(np.prod(roidims)))
	ref = ref[np.where(roi==1)]

	# Deal with sign ambiquity and, if requested, normalize y to range between 0 and 1
	for evec in range(0,y.shape[1]):		
		y[:,evec] = np.multiply(y[:,evec],np.sign(np.corrcoef(y[:,evec],ref)[0,1]))
		if norm_flag: 
			tmp = y[:,evec] - min(y[:,evec])
			y[:,evec] = np.divide(tmp,max(tmp))

	# Store the eigenmaps as a 4D nifti image
	print('Writing connectopic maps to: ' + outdir)
	outfile = outdir + "/" + out_base_name + ".cmaps.nii.gz"
	yDat = np.zeros(shape=roidims+(nmaps,))
	yDat = np.reshape(yDat,(np.prod(roidims),nmaps))
	yDat[roiIndices,:] = y[:,1:nmaps+1]
	yDat = np.reshape(yDat,roidims+(nmaps,))
	yImg = nib.Nifti1Image(yDat,roiImg.get_affine(),roiImg.get_header())
	try:
		nib.save(yImg,outfile)
	except:
		sys.exit('Cannot save ' + outfile | '\nExiting.')

	# Optionally project eigenmaps onto mask by spatial regression
	if proj_flag:
		print('Computing projections onto mask...')
		outfile = outdir + "/" + out_base_name + ".pmaps.nii.gz"		
		YDat = np.zeros(shape=roidims+(nmaps,))
		YDat = np.reshape(YDat,(np.prod(roidims),nmaps))
		for evec in range(1,y.shape[1]):	
			X = np.vstack([np.ones(y.shape[0]),y[:,evec].T])
			beta = np.dot(np.linalg.pinv(X.T),A.T)	
			Y = np.dot(B.T,beta.T)[:,1]
			if norm_flag:
				Y -= min(Y)
				Y /= max(Y)
			YDat[maskIndices,evec-1] = Y
		print('Writing projection maps to: ' + outdir)
		YDat = np.reshape(YDat,roidims+(nmaps,))
		YImg = nib.Nifti1Image(YDat,roiImg.get_affine(),roiImg.get_header())
		try:
			nib.save(YImg,outfile)
		except:
			sys.exit('Cannot save ' + outfile | '\nExiting.')

	print("Done.")


if __name__ == "__main__":

	import argparse
	parser = argparse.ArgumentParser(description="ConGrads")
	parser.add_argument("-i",dest="infiles",help="One or more 4D images",required=True,nargs='*')
	parser.add_argument("-r",dest="roifile",help="Region-of-Interest (binary 3D nifti)",required=True)
	parser.add_argument("-m",dest="maskfile",help="Mask (binary 3D nifti)",required=True)
	parser.add_argument("-o",dest="outdir",help="Output directory",required=True)
	parser.add_argument("--nmaps",dest="nmaps",default=1,type=int,help="Number of connectopic maps")
	parser.add_argument("--save_eta2",dest="save_eta2",action="store_true",help="Store eta2 matrix")
	parser.add_argument("--norm",dest="norm_flag",action="store_true",help="Normalise maps")
	parser.add_argument("--project",dest="proj_flag",action="store_true",help="Project maps onto mask")
	args=parser.parse_args()
	main(args.infiles,args.roifile,args.maskfile,args.outdir,args.nmaps,args.save_eta2,args.norm_flag,args.proj_flag)

