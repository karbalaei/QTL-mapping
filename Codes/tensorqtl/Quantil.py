from limix._bits import dask, numpy, pandas, xarray
from numpy import issubdtype, integer, asarray

def quantile_gaussianize(X, axis=1, inplace=False):
    r"""Normalize a sequence of values via rank and Normal c.d.f.
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
    """
    if hasattr(X, "dtype") and issubdtype(X.dtype, integer):
        raise ValueError("Integer type is not supported.")
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
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        D = nanrankdata(D)
        D = D / (isfinite(D).sum(axis=0) + 1)
        D = apply_along_axis(norm.isf, 0, D)
    D = D.swapaxes(1, axis)
    X[:] = D
    return X.reshape(orig_shape)