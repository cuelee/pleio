import numpy as np

def write_mat(X,fname):
    np.savetxt(fname=fname,X=X,fmt='%1.5f',delimiter=' ',newline='\n',header='',footer='',comments='#')
    return();

def load_mat(fname):
    X = np.loadtxt(fname=fname);
    return(X);

