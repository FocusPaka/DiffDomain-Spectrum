import sys
import gzip
import numpy as np
import pandas as pd
from collections import defaultdict
# import TracyWidom as tw
# from skrmt.ensemble.spectral_law import WignerSemicircleDistribution
# from scipy.stats import chisquare
import hicstraw 
import warnings
warnings.filterwarnings('ignore')
from joblib import Parallel, delayed


# =============================================================================
# DiffDomain-Spectrum core module.
# Merged from the former: functions.py + spectrum_test_updated.py + WignerSemicircle.py
# =============================================================================

# --- Wigner semicircle distribution (was WignerSemicircle.py) ---

class WignerSemicircle:
    def __init__(self, R=2.0, loc=0.0, scale=1.0):
        self.R = R
        self.loc = loc
        self.scale = scale

    def pdf(self, x):
        x = (np.asarray(x) - self.loc) / self.scale
        mask = np.abs(x) <= self.R
        return np.where(mask, 2 / (np.pi * self.R**2) * np.sqrt(self.R**2 - x**2) / self.scale, 0)

    def cdf(self, x):
        x = (np.asarray(x) - self.loc) / self.scale
        x = np.clip(x, -self.R, self.R)
        return 0.5 + (x * np.sqrt(self.R**2 - x**2) + self.R**2 * np.arcsin(x / self.R)) / (np.pi * self.R**2)
    
    def ppf(self, p):
        p = np.clip(p, 1e-10, 1 - 1e-10)
        grid = np.linspace(-self.R, self.R, 2000)
        cdf_vals = self.cdf(grid)
        return np.interp(p, cdf_vals, grid) * self.scale + self.loc

    def fit_loc_scale(self, data):
        xmin, xmax = np.min(data), np.max(data)
        self.loc = 0.5 * (xmin + xmax)
        self.scale = 0.5 * (xmax - xmin) / self.R
        return self.loc, self.scale

    def rvs(self, size=1, random_state=None):
        if random_state is not None:
            np.random.seed(random_state)
        u = np.random.beta(1.5, 1.5, size=size)
        x = 2 * self.R * u - self.R
        return self.loc + self.scale * x


# --- Full-spectrum test, parallel Monte-Carlo (was spectrum_test_updated.py) ---

def spectrum_test_parallel(Mat, N=10000, n_jobs=-1, random_state=None, fit_loc_scale=False):
    nd = Mat.shape[0]
    Mat = Mat / np.sqrt(nd)
    w = np.linalg.eigvalsh(Mat)
    w_sorted = np.sort(w)
    n1 = len(w_sorted)

    wsd = WignerSemicircle(R=2.0)
    if fit_loc_scale:
        wsd.fit_loc_scale(w_sorted)
     
    p1 = np.clip(wsd.cdf(w_sorted), 1e-10, 1 - 1e-10)
    p2 = 1 - p1
    ranks = np.arange(1, n1 + 1)
    h = np.log(p1) / (n1 - ranks + 0.5) + np.log(p2) / (ranks - 0.5)
    ZA_obs = -np.sum(h)

    lower_bound = max(-2,w_sorted[0])  
    upper_bound = min(w_sorted[-1],2)
    if upper_bound - lower_bound < 1e-8:
        return ZA_obs, 1.0, n1
   
    def single_mc(seed):
        rng = np.random.default_rng(seed)
        sim = wsd.rvs(n1)

        if fit_loc_scale:  
            wsd_sim = WignerSemicircle(R=wsd.R)
            wsd_sim.fit_loc_scale(sim)
        else:  # plug-in
            wsd_sim = wsd

        sim = np.sort(sim)
        p1_sim = np.clip(wsd_sim.cdf(sim), 1e-10, 1 - 1e-10)
        p2_sim = 1 - p1_sim
        h_sim = np.log(p1_sim) / (n1 - ranks + 0.5) + np.log(p2_sim) / (ranks - 0.5)
        ZA_sim = -np.sum(h_sim)
        return ZA_sim

    seeds = np.arange(N) + (0 if random_state is None else random_state)
    ZA_sims = Parallel(n_jobs=n_jobs, prefer="processes")(
        delayed(single_mc)(s) for s in seeds
    )

    ZA_sims = np.array(ZA_sims)
    p_value = np.mean(ZA_sims >= ZA_obs)

    return ZA_obs, p_value, n1


# --- Data loading, matrix extraction, normalization (was functions.py) ---

# Part 1: Load TAD list
def finftypes(filepath):
    """
    filepath can be stdin(standard input), normal txt, and gz file
    """
    if(filepath == 'stdin'):
        fin = sys.stdin
    elif(filepath[-3:] == '.gz'):
        fin = gzip.open(filepath)
    else:
        fin = open(filepath,'r')
    return fin


def loadtads(inpath, sep, chrnum, reso, min_nbin): 
    tadb = []
    fin = finftypes(inpath)
    next(fin)
    for line in fin:
        col = line.rstrip().split(sep)
        tadb.append(col[:3])
    fin.close()

    tadb = pd.DataFrame(tadb)
    tadb.iloc[:,0] = tadb.iloc[:,0].str.replace('chr', '')
    tadb.iloc[:,1:3] = tadb.iloc[:,1:3].astype(int)

    if chrnum != 'ALL':
        ind = tadb[0] != chrnum    # 只去掉染色体为chrnum的行
        tadb = tadb[ind] #tadb = [_ for _ in tadb if _[0]==chrnum]

    tadb = tadb[tadb.iloc[:,1] < tadb.iloc[:,2]].copy()
    # select domain with at least min_nbin
    nbins = np.ceil(tadb[2] /reso) - np.floor(tadb[1] / reso) 
    print(tadb.shape)
    # tadb = tadb[nbins>=min_nbin] #huadm
    # print(tadb.shape)

    return tadb


# Part 2: Extract matrix from Hi-C dataset
def makewindow2(start, end, reso):
    if start >= end :
        print('The end of TAD should be bigger than its start')
        return None

    k1 = start // reso  
    k2 = end // reso 

    wins = [ _*reso for _ in range(k1, k2+1)]
    return wins


def contact_matrix_from_hic(chrn, start, end, reso, fhic, hicnorm):
    domwin = makewindow2(start, end, reso)
    nb = len(domwin)
    domwin_dict = {domwin[_]:_ for _ in range(nb)}
    mat = np.empty(shape=(nb, nb))
    mat[:] = np.nan
    region = "{0}:{1}:{2}".format(chrn, domwin[0], domwin[-1])
    region2 = "{0}:{1}-{2}".format(chrn, domwin[0], domwin[-1])
    print(region)

    if fhic.endswith('.hic'):
        try:
            el = hicstraw.straw('observed',hicnorm, fhic, region, region, 'BP', reso)
            
            for i in range(len(el)):
                bin0 = el[i].binX
                bin1 = el[i].binY
                k=domwin_dict[bin0]
                l=domwin_dict[bin1]
                if k == l:
                    mat[k, l] = el[i].counts
                else:
                    mat[k, l] = el[i].counts
                    mat[l, k] = el[i].counts
            print(mat)

            return mat
        except IOError:
            print("Sorry,{%s} does't exist."% fhic )
            mat = None
            return mat

    
    elif fhic.endswith('cool') : # .cool or .mcool
        import cooler
        import subprocess as sp

        if fhic.endswith('.cool'):
            hicfile = fhic
        elif fhic.endswith('.mcool'):
            hicfile = f'{fhic}::resolutions/{reso}'
        
        c = cooler.Cooler(f'{hicfile}')

        if domwin[0] < c.chromsizes[chrn] and domwin[-1] < c.chromsizes[chrn]:
                    mat = c.matrix(balance=False).fetch(region2)
                    print(mat)
        else:
            mat = None
        return mat
        
    # load a sparse matrix with three columns
    else:
        print('Trying to sparse the files as three columns by "\t" ')
        data = pd.read_table(fhic,sep='\t')
        # filter interactiosn by only keeping those within the given TAD region
        ind0 = data.iloc[:, 0] >= int(start) 
        ind1 = data.iloc[:, 1] <= int(end)
        ind = ind0 & ind1
        data = data.loc[ind, ]

        el =[[],[],[]]
        el[0] = data[data.columns[0]].tolist()
        el[1] = data[data.columns[1]].tolist()
        el[2] = data[data.columns[2]].tolist()
        
        for i in range(len(el[2])):
            bin0 = el[0][i]
            bin1 = el[1][i]
            k=domwin_dict[bin0]
            l=domwin_dict[bin1]
            if k == l:
                mat[k, l] = el[2][i]
            else:
                mat[k, l] = el[2][i]
                mat[l, k] = el[2][i]
       
        return mat



# Part 3: Compare two domains
def compute_nbins(start, end, reso):
    domwin = makewindow2(start, end, reso)
    # create index for domwin
    nb = len(domwin)
    return nb

def extractKdiagonalCsrMatrix(spsCsrMat):
    nonZeroIndex = spsCsrMat.nonzero()
    contactByDistance = defaultdict(list)
    for _ in range(len(nonZeroIndex[0])):
        rowIndex = nonZeroIndex[0][_]
        colIndex = nonZeroIndex[1][_]
        if colIndex >= rowIndex:
            dist = colIndex - rowIndex
            contactByDistance[dist].append(spsCsrMat[rowIndex, colIndex])

    return contactByDistance

def normDiffbyMeanSD(D):
    # log transformation
    D = np.log(D)
    contactByDistanceDiff = extractKdiagonalCsrMatrix(D)
    # imputation of nan and inf
    # part0: get the median and maximum for each off-diagonal
    a, b = defaultdict(list), defaultdict(list)  # a: median, b: max, m: mean, sd: std
    m, sd = defaultdict(float), defaultdict(float)
    all_nan_k_list = [] 

    for k, val in contactByDistanceDiff.items():
        indnan = np.isnan(val)
        indinf = np.isinf(val)
        if np.all(indnan):
            all_nan_k_list.append(k)
        elif np.any(indnan) and np.any(indinf):
            val = np.array(val)
            ind = np.logical_or(indnan, indinf)
            val1 = val[np.logical_not(ind)]
            a[k] = np.median(val1)
            b[k] = np.max(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        elif np.any(indnan) and not np.any(indinf):
            val = np.array(val)
            val1 = val[np.logical_not(indnan)]
            a[k] = np.median(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        elif not np.any(indnan) and np.any(indinf):
            val = np.array(val)
            val1 = val[np.logical_not(indinf)]
            b[k] = np.max(val1)
            m[k] = np.mean(val1)
            sd[k] = np.std(val1)
        else:
            m[k] = np.mean(val)
            sd[k] = np.std(val)

    #part1: impute nan and inf
    indnan = np.isnan(D)
    indnan.astype(int)
    indr, indc = np.nonzero(indnan)
    for _ in range(len(indr)):
        k = abs(indr[_] - indc[_])
        if k in m.keys():
            D[indr[_], indc[_]] = a[k]

    posinf = D == np.inf
    neginf = D == -np.inf

    indr, indc = np.nonzero(posinf)
    for _ in range(len(indr)):
        k = abs(indr[_] - indc[_])
        D[indr[_], indc[_]] = b[k]

    indr, indc = np.nonzero(neginf)
    for _ in range(len(indr)):
        k = abs(indr[_] - indc[_])
        D[indr[_], indc[_]] = -b[k]

    # subtract mean and dividing by sd
    for k, v in sd.items():
        if np.isnan(v) or v == 0:
            sd[k] = 1

    # print(m.keys())
    nd = D.shape[0]
    for i in range(nd):
        for j in range(nd):
            k = abs(i-j)
            if k in m.keys():
                D[i, j] = (D[i,j] - m[k]) / sd[k]
            if k in all_nan_k_list:
                # impute by a random variable sampled from a N(0,1)
                D[i, j]=D[j, i] = np.random.randn()
    D = np.clip(D, -3, 3)

    return D, all_nan_k_list



def visualization(chrn, start, end, reso, hicnorm, fhic0, fhic1):

        mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
        mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)
        return mat0,mat1



def comp2domins_by_spectrum(chrn, start, end, reso, hicnorm, fhic0, fhic1, min_nbin, f, N):
    mat0 = contact_matrix_from_hic(chrn, start, end, reso, fhic0, hicnorm)
    mat1 = contact_matrix_from_hic(chrn, start, end, reso, fhic1, hicnorm)

    if not mat0 is None and not mat1 is None:
        # remove rows that have more than half np.nan
        # nbins = compute_nbins(start, end, reso)
        nbins = mat0.shape[0]
        mat0 = np.where(mat0==0,np.nan,mat0)
        mat1 = np.where(mat1==0,np.nan,mat1)
        ind0 = np.sum(np.isnan(mat0), axis=0) < nbins * (1-float(f))  
        ind1 = np.sum(np.isnan(mat1), axis=0) < nbins * (1-float(f))
        ind = ind0 & ind1

        indarray = np.array(ind)
        mat0rmna = mat0[indarray, :]
        mat0rmna = mat0rmna[:, indarray]
        mat1rmna = mat1[indarray, :]
        mat1rmna = mat1rmna[:, indarray]
        print('The number of bins after removing nan rows/columns:', ind.sum())

        if ind.sum() >= min_nbin:
            # compute the differnece matrix
            Diffmat = mat0rmna / mat1rmna
            # remove the nan and inf values

            Diffmatnorm,all_nan_k_list = normDiffbyMeanSD(D=Diffmat)
            result = spectrum_test_parallel(Diffmatnorm,N=N)
            domname = '%s:%s-%s' % (chrn, start, end)
            result = [chrn, start, end, domname, result[1], result[2], result[0]]
            print('The number of nan k:', len(all_nan_k_list))

        else:
            domname = '%s:%s-%s' % (chrn, start, end)
            print('The matrix is too spase or too small at this resolution to be calculated !')
            result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
    else:
        domname = '%s:%s-%s' % (chrn, start, end)
        result = [chrn, start, end, domname, np.nan, np.nan, np.nan]
        print('The matrixes failed to be load!')

    return result
