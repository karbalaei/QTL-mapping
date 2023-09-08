import pandas as pd
import xarray as xr
from numpy.linalg import cholesky
from pandas_plink import read_plink1_bin
import os
import numpy as np
from cellregmap import run_association, run_interaction, estimate_betas
from itertools import repeat
from joblib import Parallel, delayed
from PIL import ImageDraw, Image
import numpy as np
from pathlib import Path
from time import sleep, time
from multiprocessing import cpu_count
import os.path


###### functions ######
from limix._bits import dask, numpy, pandas, xarray



def quantile_gaussianize(X, axis=1, inplace=False):
    r'''Normalize a sequence of values via rank and Normal c.d.f.

    It defaults to column-wise normalization.

    Parameters
    ----------
    X : array_like
        Array of values.
    axis : int, optional
        Axis value. Defaults to `1`.
    inplace : bool, optional
        Defaults to `False`.

    Returns
    -------
    array_like
        Gaussian-normalized values.

    Examples
    --------
    .. doctest::

        >>> from limix.qc import quantile_gaussianize
        >>> from numpy import array_str
        >>>
        >>> qg = quantile_gaussianize([-1, 0, 2])
        >>> print(qg) # doctest: +FLOAT_CMP
        [-0.67448975  0.          0.67448975]
    '''
    from numpy import issubdtype, integer, asarray

    if hasattr(X, 'dtype') and issubdtype(X.dtype, integer):
        raise ValueError('Integer type is not supported.')

    if isinstance(X, (tuple, list)):
        if inplace:
            raise ValueError("Can't use `inplace=True` for {}.".format(type(X)))
        X = asarray(X, float)

    if numpy.is_array(X):
        X = _qg_numpy(X, axis, inplace)
    elif pandas.is_series(X):
        X = _qg_pandas_series(X, axis, inplace)
    elif pandas.is_dataframe(X):
        X = _qg_pandas_dataframe(X, axis, inplace)
    elif dask.is_array(X):
        X = _qg_dask_array(X, axis, inplace)
    elif dask.is_series(X):
        raise NotImplementedError()
    elif dask.is_dataframe(X):
        X = _qg_dask_dataframe(X, axis, inplace)
    elif xarray.is_dataarray(X):
        X = _qg_xarray_dataarray(X, axis, inplace)
    else:
        raise NotImplementedError()

    return X


def _qg_numpy(X, axis, inplace):
    from scipy.stats import norm
    from numpy import isfinite
    from numpy.ma import masked_invalid
    from numpy_sugar import nanrankdata
    from numpy import apply_along_axis
    import warnings

    orig_shape = X.shape
    if X.ndim == 1:
        X = X.reshape(orig_shape + (1,))

    if not inplace:
        X = X.copy()

    D = X.swapaxes(1, axis)
    D = masked_invalid(D)
    D *= -1

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        D = nanrankdata(D)
        D = D / (isfinite(D).sum(axis=0) + 1)
        D = apply_along_axis(norm.isf, 0, D)

    D = D.swapaxes(1, axis)
    X[:] = D

    return X.reshape(orig_shape)


def _qg_pandas_series(X, axis, inplace):
    if not inplace:
        X = X.copy()

    a = X.to_numpy()
    _qg_numpy(a, axis, True)
    X[:] = a

    return X


def _qg_pandas_dataframe(x, axis, inplace):
    if not inplace:
        x = x.copy()

    a = x.to_numpy()
    _qg_numpy(a, axis, True)
    x[:] = a

    return x


def _qg_dask_array(x, axis, inplace):
    import dask.array as da
    from scipy.stats import norm
    from numpy_sugar import nanrankdata

    if inplace:
        raise NotImplementedError()

    x = x.swapaxes(1, axis)

    x = dask.array_shape_reveal(x)
    shape = da.compute(*x.shape)
    x = da.ma.masked_array(x)
    x *= -1
    x = da.apply_along_axis(_dask_apply, 0, x, nanrankdata, shape[0])
    x = x / (da.isfinite(x).sum(axis=0) + 1)
    x = da.apply_along_axis(_dask_apply, 0, x, norm.isf, shape[0])

    return x.swapaxes(1, axis)


def _qg_dask_dataframe(x, axis, inplace):
    if inplace:
        raise NotImplementedError()

    d = x.to_dask_array(lengths=True)
    orig_chunks = d.chunks
    d = _qg_dask_array(d, axis, False).rechunk(orig_chunks)
    return d.to_dask_dataframe(columns=x.columns, index=x.index)


def _qg_xarray_dataarray(X, axis, inplace):
    if not inplace:
        X = X.copy(deep=True)

    data = X.data

    if dask.is_array(data):
        data = _qg_dask_array(data, axis, inplace)
    else:
        data = _qg_numpy(data, axis, inplace)

    X.data = data
    return X


def _dask_apply(x, func1d, length):
    from numpy import resize
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        x = func1d(x)

    return resize(x, length)



os.getcwd()
os.chdir('/data/fgoes1/RK/CellRegMap')
os.getcwd()


## genotype_individual_id are donor IDs, as found in the genotype matrix (G) and GRM covariance (K)
## phenotype_sample_id are cell IDs, as found in the scRNA-seq phenotype vector (y) and cell context covariance (C)

sample_mapping = pd.read_csv('excit_sample_mapping_file.csv', dtype={'genotype_individual_id': str, 'phenotype_sample_id': str})


sample_mapping.head()


## extract unique individuals
donors = sample_mapping['genotype_individual_id'].unique()
donors.sort()
print('Number of unique donors: {}'.format(len(donors)))


#kinship file
kinship_file =  'Kinship_matrix_05_geno.txt'
K = pd.read_csv( kinship_file, sep='	', index_col=0)
assert all(K.columns == K.index) #symmetric matrix, donors x donors

K


K = xr.DataArray(K.values, dims=['sample_0', 'sample_1'], coords={'sample_0': K.columns, 'sample_1': K.index})

K = K.sortby('sample_0').sortby('sample_1')
donors = sorted(set(list(K.sample_0.values)).intersection(donors))
print('Number of donors after kinship intersection: {}'.format(len(donors)))


## and decompose such as K = hK @ hK.T (using Cholesky decomposition)
hK = cholesky(K.values)
hK = xr.DataArray(hK, dims=['sample', 'col'], coords={'sample': K.sample_0.values})
assert all(hK.sample.values == K.sample_0.values)


## use sel from xarray to expand hK (using the sample mapping file)
hK_expanded = hK.sel(sample=sample_mapping['genotype_individual_id'].values)
assert all(hK_expanded.sample.values == sample_mapping['genotype_individual_id'].values)


######################################
############### Genotypes ############
######################################
plink_file = 'MDD_control_sc_maf_filtered_05_geno.bed'


G = read_plink1_bin(plink_file)


# Cordinates (PCAs)
E_file = 'excit_CellRegMap_PCAs.txt'
E = pd.read_csv(E_file, index_col=0 ,  sep='	')

E = xr.DataArray(E.values, dims=['cell', 'pc'], coords={'cell': E.index.values, 'pc': E.columns.values})

E.head()


#Annotation file
annotation_file = 'excit_annotation.txt'
anno_df = pd.read_csv(annotation_file, sep='	')
anno_df.head(2)

#make definition
def cis_snp_selection(feature_id, annotation_df, G, window_size):
    feature = anno_df[anno_df['feature_id']==feature_id].squeeze()
    chrom = str(feature['chromosome'])
    start = feature['start']
    end = feature['end']
    # make robust to features self-specified back-to-front
    lowest = min([start,end])
    highest = max([start,end])
    # for cis, we sequentially add snps that fall within each region
    G = G.where((G.chrom == str(chrom)) & (G.pos > (lowest-window_size)) & (G.pos < (highest+window_size)), drop=True)
    return G



phenotype_file = 'excit_GE.txt'
phenotype = pd.read_csv(phenotype_file, index_col=0, sep='	')
phenotype = phenotype.astype(float)

phenotype.head()


gene_list = np.array(phenotype.index)
gene_list2 = gene_list[17000:17999]

phenotype = xr.DataArray(phenotype.values, dims=['trait', 'cell'], coords={'trait': phenotype.index.values, 'cell': phenotype.columns.values})

phenotype = phenotype.sel(cell=sample_mapping['phenotype_sample_id'].values)
assert all(phenotype.cell.values == sample_mapping['phenotype_sample_id'].values)

phenotype.head()

# window size (cis)
w = 1000000




W= pd.read_csv('excit_CellRegMap_covariate.txt', index_col=0, sep='	')


results = pd.DataFrame ()

def CellRegMap_interaction(gene,anno_df, G, w, W, E , hK , results):
    print(gene)
    G_sel = cis_snp_selection(gene, anno_df, G, w)
    G_expanded = G_sel.sel(sample=sample_mapping['genotype_individual_id'].values)
    y = phenotype.sel(trait=gene)
    y= quantile_gaussianize(y)
    pv = run_association(y=y, G=G_expanded, W=W, E=E, hK=hK_expanded)[0]
    s= gene
    n= pv.size
    Id= np.array(list(repeat(s, n)))
    A0= G_expanded.a0.values
    A1= G_expanded.a1.values
    SNP= G_expanded.snp.values
    Chromosome= G_expanded.chrom.values
    Position = G_expanded.pos.values
    df = (pd.DataFrame ({'ID':Id , 'SNP': SNP ,'Chromosome': Chromosome, 'Position': Position , 'A0': A0 , 'A1':A1  ,'Pvalues':pv}))
    df.to_csv('/data/fgoes1/RK/CellRegMap/excit/'+ gene + '_excit_Association.csv',sep=',',index=False)


Parallel(n_jobs=int(cpu_count()), prefer='processes')(
     delayed(CellRegMap_interaction)(gene = gene ,anno_df = anno_df, G= G, w= w, W = W, E=E , hK=hK_expanded , results = results)
     for gene in gene_list2)
