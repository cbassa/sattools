import csv
from beyond.io.tle import Tle
from astropy import units as u


def tle_from_file(filename):
    '''
    Returns: TLE:Tle, Object name:str
    '''

    with open(filename, 'r') as f:
        lines = f.readlines()

    return Tle(''.join(lines)), lines[0].strip()


def read_site_file(filename):
    sites = {}

    with open(filename, newline='') as f:
        reader = csv.reader(f,
                            delimiter=' ',
                            quoting=csv.QUOTE_NONE,
                            skipinitialspace=True)
        for row in reader:
            if row[0][0] == '#':
                continue

            sites[int(row[0])] = {'symbol': row[1],
                                  'latitude': float(row[2]) * u.deg,
                                  'longitude': float(row[3]) * u.deg,
                                  'altitude': float(row[4]) * u.m,
                                  'name': ' '.join(row[5:])}
    return sites
