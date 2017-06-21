#Modified from agdc/agdc-v2-examples/notebooks/06_simple_change_detection_using_annual_mean_NDVI.ipynb
#Retrieve the NBAR, FC and PQ data 

import os, sys
import warnings
warnings.simplefilter("ignore", RuntimeWarning)

import fiona
import pandas
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

def get_data(query, mask_components, pnbars, pfcs, pqas):
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

    year = sys.argv[1]
    dc = datacube.Datacube(app='koala-obs')

    ws = os.path.join(os.path.abspath('..'), 'koala_obs')

    inshp = 'koala_obs_1995_2009.shp'
    errshp = 'koala_err_%s.shp'%year  #In case there's a dodgy point or two
    outcsv = 'koala_data_%s.csv'%year #Output csv

    sensors = [
        'ls7',
        'ls5'
    ] 

    pnbar_measurements=['red', 'nir']
    pfc_measurements=['BS']

    query = {
        'output_crs':'EPSG:4326',
        'resolution':(-0.00025, 0.00025),
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
            out.write('time,latitude,longitude,NDVI,SAVI,OID' + '\n')

    with fiona.open(inshp, 'r') as source, open(outcsv, 'a') as out:

        recs = filter(lambda r:  
                      int(r['properties']['DYEAR']) == int(year) and int(r['properties']['OBJECTID']) > last_oid, 
                      source)

        driver = source.driver
        crs = source.crs
        schema = source.schema.copy()
        schema['properties']['ORIGID']=schema['properties']['OBJECTID']

        
        mode = 'a' if os.path.exists(errshp) else 'w'
        with fiona.open(errshp, mode,
                        driver=driver,
                        crs=crs,
                        schema=schema) as err:

            for i,rec in enumerate(recs):

                x,y = rec['geometry']['coordinates']
                oid = rec['properties']['OBJECTID']
                rec['properties']['ORIGID'] = oid

                print('%s\t%s\t%s'%(year,i,oid))

                query['time'] = ('%s-01-01'%year, '%s-12-31'%year)
                query['lon'] = (x-0.00025,x+0.00025)
                query['lat'] = (y-0.00025,y+0.00025)

                try:
                    nbar, fc = get_data(query, mask_components, pnbars, pfcs, pqas)
                    px,py=(int(xy) for xy in ~nbar.affine * (x,y))

                    nbar = nbar.isel(longitude=px,latitude=py).dropna('time', how = 'any')
                    fc = fc.isel(longitude=px,latitude=py).dropna('time', how = 'any')

                    red, nir, lfac = nbar.red/10000, nbar.nir/10000, fc.BS/100
                    ndvi = ((nir-red)/(nir+red))
                    savi = ((nir-red)/(nir+red+lfac))*(1+lfac)

                    ndvids=ndvi.to_dataset(name='NDVI')
                    savids=savi.to_dataset(name='SAVI')
                    ds = xarray.merge([ndvids,savids])
                    ds['OID'] = oid
                    
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
