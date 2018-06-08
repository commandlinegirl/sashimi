def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    Adapted from:
    http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    see also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if window_len<3:
        return x

    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    print("Smoothing window")
    print(w)
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


def smooth_salmon_output(df, smooth_raw_output, window_len, smoothing_strategy):
    '''
    Input df columns: Chromosome Start End NumReads TPM
    Returns: df with appended 'Smoothed_TPM' column
    '''
    #TODO: do not smooth over chromosome boundaries, for now it's OK, since
    # TPM at the boundaries are all 0
    if smooth_raw_output:
        if smoothing_strategy == 'medianfilter':
            df['Smoothed_TPM'] = signal.medfilt(df['TPM'].tolist(), window_len)
        else:
            df['Smoothed_TPM'] = smooth(df['TPM'].tolist(), window_len, smoothing_strategy)[int(window_len/2):(int(window_len/2)*-1)]
    else:
        df['Smoothed_TPM'] = df['TPM']
    return df
