import numpy as np

def LS_chi(betas, stders, cor):
	bes = list(map(float,betas))
	C = np.matrix(cor, dtype=float)
	stds_np = np.matrix( np.diag( list(map(float,stders)) ) )
	
	V = stds_np.dot(C).dot(stds_np)
	Vinv = np.linalg.inv(V)
	ones = np.matrix(list(map(float,[1]*len(bes))))
	
	newv = 1 / (ones.dot(Vinv).dot(ones.transpose()))
	newx = np.matrix(ones).dot(Vinv).dot(np.matrix(bes).transpose()) / (ones.dot(Vinv).dot(ones.transpose()))	
	# newstd = np.sqrt(newv)
	chisqs = newx**2/newv
	return(float(chisqs)) 

def LS_apply(LS_chi, U, R, n, X, row_wise):
    if(row_wise == True):
        niter = X.shape[0];
        lr_vec = [0.0]*niter;
        stders = [1]*n;
        cor = np.diag([1]*n);
        for i in range(niter):
            beta = X[i,:];
            lr_vec[i] = LS_chi(betas = beta, stders = stders, cor = cor);
    elif(row_wise == False):
        niter = X.shape[1];
        lr_vec = [0.0]*niter;
        for i in range(niter):
            beta = X[:,i];
            lr_vec[i] = LS_chi(betas = beta, stders = stders, cor = cor);
    else:
        quit("FATAL ERROR")
    return(lr_vec)
