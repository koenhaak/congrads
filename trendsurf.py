# usage: trendsurf.py [-h] -i FILENAME -r MASKFILE -o OUTDIR 
#                     [-b <int> BASIS]
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
#
# Adapted from https://github.com/amarquand/nispat/blob/master/nispat/trendsurf.py

import numpy as np
import nibabel as nib

def load_data(datafile,maskfile):

    dataImg = nib.load(datafile)
    data = dataImg.get_data()
    dim = data.shape
    if len(dim) <= 3:
        dim = dim + (1,)
    data = np.reshape(data,(np.prod(dim[0:3]),dim[3]))

    maskImg = nib.load(maskfile)
    mask = maskImg.get_data()
    mask = np.reshape(mask,(np.prod(dim[0:3])))
    maskIndices = np.where(mask==1)[0]
    data = data[maskIndices,:]
	
    i,j,k = np.meshgrid(np.linspace(0, dim[0]-1, dim[0]),
                        np.linspace(0, dim[1]-1, dim[1]),
                        np.linspace(0, dim[2]-1, dim[2]), indexing='ij')

    world = np.vstack((i.ravel(),j.ravel(),k.ravel(),np.ones(np.prod(i.shape),float))).T
    world = np.dot(world,dataImg.affine.T)[maskIndices,0:3]
    
    return data,world,maskIndices

def save_nifti(data,filename,examplenii,maskIndices):

    # load example image
    ex_img = nib.load(examplenii)
    ex_img.shape
    dim = ex_img.shape[0:3]
    if len(data.shape) < 2:
        nvol = 1
        data = data[:, np.newaxis]
    else:
        nvol = int(data.shape[1])

    # write data
    array_data = np.zeros((np.prod(dim),nvol))
    array_data[maskIndices,:] = data
    array_data = np.reshape(array_data,dim+(nvol,))
    array_img = nib.Nifti1Image(array_data,ex_img.affine,ex_img.header)
    nib.save(array_img,filename)

def create_basis(X,dimpoly):

    dimx = X.shape[1]
    print('Generating polynomial basis set of degree',str(dimpoly),'...')
    Phi = np.zeros((X.shape[0],X.shape[1]*dimpoly))
    colid = np.arange(0,dimx)
    for d in range(1, dimpoly+1):
        Phi[:,colid] = X**d
        colid += dimx

    return Phi

def main(filename,maskfile,outdir,basis,ard=False):

    """ :outputs: * yhat - predictive mean
                  * ys2 - predictive variance
                  * trendcoeff - coefficients from the trend surface model
                  * negloglik - Negative log marginal likelihood
                  * hyp - hyperparameters
                  * explainedvar - explained variance
                  * rmse - standardised mean squared error
                  * trendcoeffvar - marginal variances """
	
    from bayesreg import BLR
    np.seterr(invalid='ignore')

    # load data
    print("Processing data in",filename)
    Y,X,maskIndices = load_data(filename,maskfile)
    Y = np.round(10000*Y) / 10000  # truncate precision to avoid numerical probs
    if len(Y.shape) == 1:
        Y = Y[:,np.newaxis]
    N = Y.shape[1]

    # standardize responses and covariates
    mY = np.mean(Y,axis=0)
    sY = np.std(Y,axis=0)
    Yz = (Y-mY)/sY
    mX = np.mean(X,axis=0)
    sX = np.std(X,axis=0)
    Xz = (X-mX)/sX

    # create basis set and set starting hyperparamters
    Phi = create_basis(Xz,basis)
    if ard is True:
        hyp0 = np.zeros(Phi.shape[1]+1)
    else:
        hyp0 = np.zeros(2)

    # estimate the models
    yhat = np.zeros_like(Yz)
    ys2  = np.zeros_like(Yz)
    nlZ  = np.zeros(N)
    hyp  = np.zeros((N,len(hyp0)))
    rmse = np.zeros(N)
    ev   = np.zeros(N)
    m    = np.zeros((N,Phi.shape[1]))
    bs2  = np.zeros((N,Phi.shape[1]))

    for i in range(0, N):
        print("Estimating model ",i+1,"of",N)
        breg = BLR()
        hyp[i,:] = breg.estimate(hyp0,Phi,Yz[:,i],'powell')
        m[i,:] = breg.m
        nlZ[i] = breg.nlZ

        # compute marginal variances
        bs2[i] = np.sqrt(np.diag(np.linalg.inv(breg.A)))

        # compute predictions and errors
        yhat[:,i],ys2[:,i] = breg.predict(hyp[i,:],Phi,Yz[:,i],Phi)
        yhat[:,i] = yhat[:,i]*sY[i] + mY[i]
        rmse[i] = np.sqrt(np.mean((Y[:,i]-yhat[:,i])**2))
        ev[i] = 100*(1-(np.var(Y[:,i]-yhat[:,i])/np.var(Y[:,i])))

        print("Variance explained =",ev[i],"% RMSE =",rmse[i])

    print("Mean (std) variance explained =",ev.mean(),"(",ev.std(),")")
    print("Mean (std) RMSE =",rmse.mean(),"(",rmse.std(),")")

    # Write output
    print("Writing output ...")
    out_base_name = outdir + "/" + filename.split('/')[-1].split('.nii')[0]
    np.savetxt(out_base_name + ".tsm.trendcoeff.txt", m, delimiter='\t', fmt='%5.8f')
    np.savetxt(out_base_name + ".tsm.negloglik.txt", nlZ, delimiter='\t', fmt='%5.8f')
    np.savetxt(out_base_name + ".tsm.hyp.txt", hyp, delimiter='\t', fmt='%5.8f')
    np.savetxt(out_base_name + ".tsm.explainedvar.txt", ev, delimiter='\t', fmt='%5.8f')
    np.savetxt(out_base_name + ".tsm.rmse.txt", rmse, delimiter='\t', fmt='%5.8f')
    np.savetxt(out_base_name + ".tsm.trendcoeffvar.txt", bs2, delimiter='\t', fmt='%5.8f')    
    save_nifti(yhat, out_base_name + '.tsm.yhat.nii.gz', filename, maskIndices)
    save_nifti(ys2, out_base_name + '.tsm.ys2.nii.gz', filename, maskIndices)
   

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Trend Surface Model")
    parser.add_argument("-i", help="input file", dest="filename",required=True)
    parser.add_argument("-r", help="mask file", dest="maskfile",required=True)
    parser.add_argument("-o", help="path to output directory", dest="outdir",required=True)
    parser.add_argument("-b <int>", help="model order", type=int, dest="basis", default=3)
    args = parser.parse_args()
    main(args.filename, args.maskfile, args.outdir, args.basis)
    

