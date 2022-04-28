def load_PEERNGA_record(filepath):
    '''
    Load record in .at2 format (PEER NGA Databases)
    Input:
        filepath : file path for the file to be load

    Returns:

        acc : vector wit the acceleration time series
        dt : time step
        npts : number of points in record
        eqname : string with year_name_station_component info
    '''

    import numpy as np

    with open(filepath) as fp:
        line = next(fp)
        line = next(fp).split(',')
        year = (line[1].split('/'))[2]
        eqname = (year + '_' + line[0].strip() + '_' +
                  line[2].strip() + '_comp_' + line[3].strip())
        line = next(fp)
        line = next(fp).split(',')
        npts = int(line[0].split('=')[1])
        dt = float(line[1].split('=')[1].split()[0])
        acc = np.array([p for l in fp for p in l.split()]).astype(float)

    return acc, dt, npts, eqname
