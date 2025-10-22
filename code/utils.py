# -*- coding: utf-8 -*-
"""
Created on Wed Jul 1 12:09:43 2020

@author: shuwentan-po

Codes that I find useful in futher codings.

To keep things simple this should only import modules from the python standard
library or numpy and scipy. 

For now the space_related class contains gsw module, will make a seperate .py

Inspired by Jesse Cusack

"""

import numpy as np
import scipy.signal as sig
import scipy.io as io
import scipy.stats as stats


class Bunch(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        

# general handy code:
def extract_num_from_str(txt):
    l = ''
    for s in txt:
        if s.isdigit() or s == '.' or s == '-':
            try:
                l += s
            except ValueError:
                pass
    return float(l)


# data_bin:
def get_blocks(var, z, z_span, z_center):
    """
    return var & z within z_center +- z_span, and indexes
    both var and z are 1d
    """
    var_blocks = [] 
    z_blocks = [] 
    var_blocks_idx = [] 
    for i in range(len(z_center)):
        mask_block = (z>=z_center[i]-z_span) & (z<=z_center[i]+z_span)
        var_blocks.append(var[mask_block])
        z_blocks.append(z[mask_block])
        var_blocks_idx.append(np.where(mask_block==1)[0])

    return var_blocks, z_blocks, var_blocks_idx



def mean_blocks(var_blocks):
    """
    compute mean value within blocks
    """
    var_mean = np.zeros(len(var_blocks),) + np.nan
    for i in range(len(var_blocks)):
        var_mean[i] = var_blocks[i].mean()    

    return var_mean



# nan_related:

def get_nans_blocks_length(self, method = True):
    """
    https://stackoverflow.com/questions/15062205/finding-start-and-stops-of-consecutive-values-block-in-python-numpy-pandas
    method = True: Returns 1D length of np.nan s block in sequence depth wise (last axis).
    method = False: Returns 1D length of ~np.nan s block in sequence depth wise (last axis).
    """
    if method:
        nan_mask = np.isnan(self)
    else:
        nan_mask = ~np.isnan(self)
    start_nans_mask = np.concatenate((np.resize(nan_mask[...,0],self.shape[:-1]+(1,)),
                                    np.logical_and(np.logical_not(nan_mask[...,:-1]), nan_mask[...,1:])
                                    ), axis=self.ndim-1)
    stop_nans_mask = np.concatenate((np.logical_and(nan_mask[...,:-1], np.logical_not(nan_mask[...,1:])),
                                np.resize(nan_mask[...,-1], self.shape[:-1]+(1,))
                                ), axis=self.ndim-1)

    start_idxs = np.where(start_nans_mask)
    stop_idxs = np.where(stop_nans_mask)

    return stop_idxs[-1] - start_idxs[-1] + 1, start_idxs[0], stop_idxs[0]



def nan_index_in_var(var):
    """ 
    Find leading and trailing indexes for nans in a 1-d array.

    Parameters
    ----------
    var : array_like

    Returns
    -------
    leading_nans_idx : numpy array
    trailing_nans_idx : numpy array
    middle_nans_idx : numpy array

    """

    if len(var.shape)>1:
        raise ValueError('var must be a 1-d array')
    else:
        mask = np.isnan(var)
        leading_nans = mask.argmin()
        trailing_nans = mask[::-1].argmin()

        leading_nans_idx = np.arange(leading_nans)
        trailing_nans_idx = np.arange(var.size - trailing_nans, var.size)

        mask[leading_nans_idx] = 0
        mask[trailing_nans_idx] = 0
        middle_nans_idx = np.where(mask==1)[0]

    return leading_nans_idx, trailing_nans_idx, middle_nans_idx
    



# time_related:
    
def sec_in_yr():
    """ return how many seconds in one year (365 days in a year convension) """
    return 3600*24*365



def datenum_to_datetime64(t, method=1):
    """
    Convert a MATLAB datenums into python datetime64['D']
    https://stackoverflow.com/questions/13965740/converting-matlabs-datenum-format-to-python

    Parameters
    ----------
    datenum : array_like
        MATLAB datenumber which is the number of days since 0000-01-00.

    method = 1 : return dt [D], day resolution
    method = 0 : return dt [us], higher resolution  

    Returns
    -------
    dt [us] : ndarray
        tim_sec = (time - time[0]).astype('timedelta64[s]').astype('float')

    """

    origin = np.datetime64('0000-01-01', 'D') - np.timedelta64(1, 'D')
    delta = np.timedelta64(1,'us') * 86400e6 

    if method:
       dt = t * np.timedelta64(1, 'D') + origin   
    else:
       dt = t * delta + origin

    return dt



# space_related:

def m_per_nm():
    """return how many meters in one nautical mile"""
    return 1851.85



def track_dis(trk_lon, trk_lat):
    """return the distance in meters between two geographical points"""
    import gsw
    trk_dis = np.zeros(len(trk_lon), )
    for i in np.arange(1, len(trk_lon)):
        trk_dis[i] = gsw.distance([trk_lon[i], trk_lon[i-1]], [trk_lat[i], trk_lat[i-1]], p=0, axis=0)

    return trk_dis 



def fine_trk(trk_lon, trk_lat, n):
    """
    Return a finer track, distance between old track points divided by n.
    
    Parameters
    ----------
    trk_lon : array_like
    trk_lat : array_like
    
    Returns
    -------
    trk_lon_new : nnumpy array
    trk_lat_new : nnumpy array

    """

    trk_lon_new = []
    trk_lat_new = []
    for i in range(len(trk_lon)-1):
        trk_lon_ = np.linspace(trk_lon[i], trk_lon[i+1], n)[0:-1]
        trk_lat_ = np.linspace(trk_lat[i], trk_lat[i+1], n)[0:-1]
        trk_lon_new = np.append(trk_lon_new, trk_lon_)
        trk_lat_new = np.append(trk_lat_new, trk_lat_)

    return trk_lon_new, trk_lat_new 


# filters-related:

def butter(cutoff, fs, btype="low", order=4):
    """Return Butterworth filter coefficients. See scipy.signal.butter for a
    more thorough documentation.

    Parameters
    ----------
    cutoff : array
        Cutoff frequency, e.g. roughly speaking, the frequency at which the
        filter acts. Units should be same as for fs paramter.
    fs : float
        Sampling frequency of signal. Units should be same as for cutoff
        parameter.
    btype : {‘lowpass’, ‘highpass’, ‘bandpass’, ‘bandstop’}, optional
        Default is 'low'.
    order : optional, int
        Default is 4. The order of the Butterworth filter.

    Returns
    -------
    b : numpy array
        Filter b coefficients.
    a : numpy array
        Filter a coefficients.

    """
    cutoff = np.asarray(cutoff)
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = sig.butter(order, normal_cutoff, btype=btype, analog=False)

    return b, a   



def butter_filter(x, cutoff, fs, btype="low", order=4, **kwargs):
    """Apply Butterworth filter to data using scipy.signal.filtfilt.

    Parameters
    ----------
    x : array
        The data to be filtered. Should be evenly sampled.
    cutoff : array
        Cutoff frequency, e.g. roughly speaking, the frequency at which the
        filter acts. Units should be same as for fs paramter.
    fs : float
        Sampling frequency of signal. Units should be same as for cutoff
        parameter.
    btype : optional, string
        Default is 'low'. Filter type can be 'low', 'high' or 'band'.
    order : optional, int
        Default is 4. The order of the Butterworth filter.

    Returns
    -------
    y : numpy array
        The filtered data.

    """
    b, a = butter(cutoff, fs, btype=btype, order=order)
    y = sig.filtfilt(b, a, x, **kwargs)

    return y


