import numpy as np
import susieinf
import pandas as pd
import pyreadr

def load_data(file_path):
    """
    Load data from an RDS file and return the data as a numpy array
    along with the number of rows and columns.
    """
    result = pyreadr.read_r(file_path)
    # Assuming there's only one object in the RDS file
    data = list(result.values())[0]
    return data.values, data.shape[0], data.shape[1]

# Input your data file path here
data_file_path = "/Users/alexmccreight/Columbia/data/X4"

verbose = True

np.random.seed(21)
Ltrue = 5
ssq = 0.01
sigmasq = 1

# Load your custom X matrix and get dimensions
X, n, p = load_data(data_file_path)

print(f"Loaded data with dimensions: {n} rows (samples) and {p} columns (features)")

# Standardize X
X = X - X.mean(axis=0)
X = X / X.std(axis=0)

LD = (X.T.dot(X))/n
b = np.zeros(p)
inds = np.random.choice(p,size=Ltrue,replace=False)
b[inds] = np.random.normal(size=Ltrue) * np.sqrt(ssq)
order = np.argsort(inds)
print('True effects: %s' % str(inds[order]))
print('Effect sizes: %s\n' % str(b[inds[order]]))

# Generate data without infinitesimal effects
print('###### Simulating without infinitesimal effects ######')
effects = X.dot(b)
y = effects + np.random.normal(size=n) * np.sqrt(sigmasq)
print('Total fraction of variance explained by SNPs: %f\n' % (np.var(effects)/np.var(y)))
z = (X.T.dot(y))/np.sqrt(n)
meansq = y.T.dot(y)/n

print('SuSiE with moments-estimated sigma^2 and tau^2')
L = 10
output = susieinf.susie(z,meansq,n,L,LD=LD,method='moments',verbose=verbose, maxiter = 100)
PIP = 1-(1-output['PIP']).prod(axis=1)
cred = susieinf.cred(output['PIP'],LD=LD)
print('SNPs with PIP > 0.1: ' + str(np.nonzero(PIP > 0.1)[0]))
print('Corresponding PIPs: %s' % str(np.round(PIP[np.nonzero(PIP > 0.1)],2)))
print('Credible sets: ' + str(cred))
print('Estimated (sigma^2,tau^2): (%f,%f)\n' % (output['sigmasq'],output['tausq']))