#Modified from agdc/agdc-v2-examples/notebooks/06_simple_change_detection_using_annual_mean_NDVI.ipynb
#Retrieve the NBAR, FC and PQ data 

import os, sys
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import fiona
from fiona.crs import from_epsg
import pandas
from pyproj import Proj, transform
import xarray

import datacube
from datacube.helpers import ga_pq_fuser
from datacube.storage import masking

class NoDataError(Exception):
    pass

def get_last_oid(csv):
    """Get last OID from CSV so we can continue where we left off"""
    df = pandas.read_csv(csv)
    return int(df.tail(1).OID)
    
def sort_data(data):
    data = xarray.concat(data, dim='time')
    time_sorted = data.time.argsort()
    data = data.isel(time=time_sorted)
    return data

def get_data(query, mask_components, pnbars, pfcs, pqas, pnbar_measurements, pfc_measurements):
    nbars = []
    fcs = []
    for pnbar, pfc, pqa in zip(pnbars,pfcs,pqas):

        #Load the NBAR, FC and corresponding PQ
        nbar = dc.load(product=pnbar, measurements=pnbar_measurements, **query)
        fc = dc.load(product=pfc, measurements=pfc_measurements, **query)
        pq = dc.load(product=pqa, fuse_func=ga_pq_fuser, **query)

        #Apply the PQ masks to the data
        try:
            cloud_free = masking.make_mask(pq, **mask_components)
            good_data = cloud_free.pixelquality.loc[query['time'][0]:query['time'][1 ]]
            nbar = nbar.where(good_data)
            fc = fc.where(good_data)
            del cloud_free, good_data
        except ValueError:
            continue

        nbars.append(nbar)
        fcs.append(fc)
        del nbar, fc, pq

    if not nbars:
        raise NoDataError
        
    #Concatenate data from different sensors together and sort so that observations are sorted by time rather
    # than sensor
    nbar = sort_data(nbars)
    fc = sort_data(fcs)

    return nbar, fc


if __name__ == '__main__':

    dc = datacube.Datacube(app='koala-obs')

    ws = os.path.join(os.path.abspath('..'), 'koala_obs')

    inshp = 'koala_obs_1995_2009.shp'
    errshp = 'koala_err_area.shp'  #In case there's a dodgy point or two
    outcsv = 'koala_data_area.csv' #Output csv

    resolution = 25
    window = (4,4)

    sensors = [
        'ls7',
        'ls5'
    ] 

    pnbar_measurements=['red', 'nir']
    pfc_measurements=['BS']

    query = {
        'group_by':'solar_day'

    }

    #PQA Mask
    mask_components = {
        'cloud_acca':'no_cloud',
        'cloud_shadow_acca' :'no_cloud_shadow',
        'cloud_fmask' :'no_cloud',
        'cloud_shadow_fmask' : 'no_cloud_shadow',
        'contiguous':True,
        'red_saturated' : False,
        'nir_saturated' : False,
    }

    pnbars = ['{}_nbar_albers'.format(s) for s in sensors]
    pfcs = ['{}_fc_albers'.format(s)   for s in sensors]
    pqas = ['{}_pq_albers'.format(s)   for s in sensors]

    # Get last OBJECTID so we can pick up where we left off after a crash or VDI session getting killed
    last_oid = -1
    try:
        last_oid = get_last_oid(outcsv)
    except (FileNotFoundError,TypeError):  # file doesn't exist, or exists but has no data
        with open(outcsv, 'w') as out:
            #out.write('time,NDVI_MEAN,NDVI_MEDIAN,SAVI_MEAN,SAVI_MEDIAN,OID,YEAR,longitude,latitude' + '\n')
            #out.write('time,RED,NIR,BS,NDVI,SAVI,OID,YEAR,longitude,latitude' + '\n')
            out.write('time,NDVI,SAVI,OID,YEAR,longitude,latitude' + '\n')

    with fiona.open(inshp, 'r') as source, open(outcsv, 'a') as out:

        recs = filter(lambda r:  
                      int(r['properties']['OBJECTID']) > last_oid, 
                      source)

        driver = source.driver
        crs = source.crs
        schema = source.schema.copy()
        schema['properties']['ORIGID']=schema['properties']['OBJECTID']

        from_srs = Proj(crs)
        to_srs = Proj(init='EPSG:3577')
        
        mode = 'a' if os.path.exists(errshp) else 'w'
        with fiona.open(errshp, mode,
                        driver=driver,
                        crs=crs,
                        schema=schema) as err:

            for i,rec in enumerate(recs):

                x,y = transform(from_srs, to_srs, *rec['geometry']['coordinates'])

                oid = rec['properties']['OBJECTID']
                year = rec['properties']['DYEAR']
                                               
                print('%s\t%s\t%s\t (%s,%s)'%(year, i, oid, x,y))

                query['time'] = ('%s-01-01'%year, '%s-12-31'%year)
                query['crs'] = 'EPSG:3577'
                query['x'] = (x-resolution*(window[0]-1)/2, x+resolution*(window[0]-1)/2)
                query['y'] = (y-resolution*(window[1]-1)/2, y+resolution*(window[1]-1)/2)

                try:
                    nbar, fc = get_data(query, mask_components, pnbars, pfcs, pqas, pnbar_measurements, pfc_measurements)

                    #Getting some unmasked "-999" (NBAR) and "-1" FC invalid values coming through
                    nbar = nbar.where((nbar.red >= 0) & (nbar.nir >= 0) & (fc.BS >= 0)).dropna('time', how = 'any')
                    fc = fc.where((nbar.red >= 0) & (nbar.nir >= 0) & (fc.BS >= 0)).dropna('time', how = 'any')

                    red, nir, lfac = nbar.red/10000, nbar.nir/10000, fc.BS/100

                    ndvi = ((nir-red)/(nir+red))
                    savi = ((nir-red)/(nir+red+lfac))*(1+lfac)

                    ndvi = ndvi.where(ndvi >= 0)
                    savi = savi.where(savi >= 0)
                    
                    #ndvimean = ndvi.mean(dim=('y', 'x')).to_dataset(name='NDVI_MEAN')
                    #ndvimedian = ndvi.median(dim=('y', 'x')).to_dataset(name='NDVI_MEDIAN')
                    #savimean = savi.mean(dim=('y', 'x')).to_dataset(name='SAVI_MEAN')
                    #savimedian = savi.median(dim=('y', 'x')).to_dataset(name='SAVI_MEDIAN')
                    #ds = xarray.merge([ndvimean, ndvimedian,savimean,savimedian])

                    #red = red.median(dim=('y', 'x')).to_dataset(name='RED')
                    #nir = nir.median(dim=('y', 'x')).to_dataset(name='NIR')
                    #lfac = lfac.median(dim=('y', 'x')).to_dataset(name='BS')
                    ndvi = ndvi.median(dim=('y', 'x')).to_dataset(name='NDVI')
                    savi = savi.median(dim=('y', 'x')).to_dataset(name='SAVI')

                    #ds = xarray.merge([red, nir,lfac, ndvi, savi])
                    ds = xarray.merge([ndvi, savi])

                    ds['OID'] = oid
                    ds['YEAR'] = year
                    ds['longitude'], ds['latitude'] = rec['geometry']['coordinates']
                    
                    df=ds.to_dataframe()
                    df.to_csv(out, index=True, header=False)
                    out.flush()
                    os.fsync(out.fileno())
                    
                except (NoDataError, ValueError):
                    print('No data for rec {} ({},{})'.format(oid,x,y))
                    err.write(rec)
                    err.flush()
                    continue
                except (Exception):
                    import traceback
                    print('Unhandled exception for rec {} ({})'.format(oid,year))
                    traceback.print_exc()
                    continue