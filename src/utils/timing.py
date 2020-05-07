import time


def current_time(scale='s', timescale=1):
    if scale == 's':
        timescale = 1
    elif scale == 'ms':
        timescale = 1000
    return int(round(time.time()*timescale))
