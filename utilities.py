from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy import coordinates
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import matplotlib.gridspec as gridspec
import aplpy
import pandas as pd
import os.path

def read7ms(path='/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl04.fits'):
    cat7 = fits.open(path)[1].data
    return cat7

def read4ms(path='/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_4ms/4Ms_Chandra_deg.fits'):
    cat4 = fits.open(path)[1].data
    return cat4

def ra2deg(hours, minutes, sec):
    return (hours + minutes/60. + sec/3600.) * 15

def dec2deg(sign, deg, arcmin, arcsec):
    return sign *(deg + arcmin/60. + arcsec/3600.)

def ra2hms(ra):
    rah = int(np.floor(ra/15.))
    ram = int(np.floor((ra/15.-rah)*60))
    ras = ((ra/15.-rah)*60-ram)*60
    print(str(rah) + 'h' + str(ram) + 'm' + str(ras))
    return 0

def dec2dms(sign, dec):
    decd = int(np.floor(dec))
    decm = int(np.floor((dec-decd)*60))
    decs = (((dec-decd)*60)-decm)*60
    print (str(sign*decd) + '°' + str(decm) + 'm' + str(decs))
    return 0

def goods_sect(ra, dec):
    a = [11, 12, 13, 14, 21, 22, 23, 24, 25, 31, 32, 33, 34, 35, 42, 43, 44, 45]
    sect=np.zeros(len(ra))
    im=[]
    w=[]
    for i in range(len(a)):
        path = '/Users/justin/Documents/Master_Thesis/data/surveys/goods/h_sz_sect' + str(a[i]) + '_v2.0_drz_img.fits'
        im.append(fits.open(path)[0])
        w.append(wcs.WCS(im[i].header))
        for j in range(len(ra)):
            x, y = w[i].all_world2pix(ra[j], dec[j], 1)
            if (x > 0.) & (x < 8192.) & (y > 0.) & (y < 8192.):
                sect[j] = a[i]
    return sect

def hugs_sect(ra, dec):
    sect=np.zeros(len(ra))
    path = '/Users/justin/Documents/Master_Thesis/data/surveys/hugs/HUGS_GOODSS_K_naturalseeing/HUGS_GOODSS_K_naturalseeing.fits'
    im = fits.open(path)[0]
    w = wcs.WCS(im.header)
    for j in range(len(ra)):
        x, y = w.all_world2pix(ra[j], dec[j], 1)
        if (x > 0.) & (x < 9373.) & (y > 0.) & (y < 10930.):
            sect[j] = 1
    return sect

def cutoutpixsize(band, arcsecsize):
    p='chandra_7ms/CDFS-7Ms-0p5to2-asca-im-bin1-01.fits'
    if band == 'GOODS':
        p='goods/h_sb_sect33_v2.0_drz_img.fits'
    if band == 'CANDELS':
        p='candels/hlsp_candels_hst_wfc3_gs-tot_f105w_v1.0_drz.fits'
    if band == 'HUGS':
        p='hugs/HUGS_GOODSS_K_naturalseeing/HUGS_GOODSS_K_naturalseeing.fits'
    if band == 'Spitzer':
        p='spitzer/SIMPLE_ECDFS_IRAC_ch1.fits'
    path = '/Users/justin/Documents/Master_Thesis/data/surveys/' + p
    q = fits.open(path)[0]
    v = wcs.WCS(q.header)
    ra = 53.08
    dec = -27.80
    center = coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    co = Cutout2D(q.data, center, size=[arcsecsize,arcsecsize]*u.arcsec, wcs=v)
    hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
    return hdu.header['NAXIS1']

def cutout(Weight=False, RMS=False, Instr='All', arcsec=10., catalogue='main'):
    if catalogue == 'main':
        path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl04.fits'
    if catalogue == 'secondary':
        path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl05.fits'
    if catalogue == 'extra':
        path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/short_table.fits'
    cat7ms = fits.open(path)[1].data
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    big_table = fits.open(path)[1].data
    sect = big_table.field('GOODS_sect')
    sect2 = big_table.field('HUGS_sect')
    #misx = big_table.field('GOODS_x_pixel')
    #misy = big_table.field('GOODS_y_pixel')
    #pix_chandra = cutoutpixsize('Chandra', arcsec)
    pix_goods = cutoutpixsize('GOODS', arcsec)
    #pix_candels = cutoutpixsize('CANDELS', arcsec)
    #pix_hugs = cutoutpixsize('HUGS', arcsec)
    #pix_spitzer = cutoutpixsize('Spitzer', arcsec)
    #chandra and candels
    if (Instr=='All') | (Instr=='Chandra+CANDELS'):
        q = []
        if Weight==False:
            p = ['chandra_4ms/CDFS-4Ms-0p5to2-asca-im-bin1.fits', 'chandra_4ms/CDFS-4Ms-0p5to8-asca-im-bin1.fits', 'chandra_4ms/CDFS-4Ms-2to8-asca-im-bin1.fits', 'chandra_7ms/CDFS-7Ms-0p5to2-asca-im-bin1-01.fits', 'chandra_7ms/CDFS-7Ms-0p5to7-asca-im-bin1-01.fits', 'chandra_7ms/CDFS-7Ms-2to7-asca-im-bin1-01.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f105w_v1.0_drz.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits']
        else:
            p = ['chandra_4ms/CDFS-4Ms-0p5to2-bin1.emap.fits', 'chandra_4ms/CDFS-4Ms-0p5to8-bin1.emap.fits', 'chandra_4ms/CDFS-4Ms-2to8-bin1.emap.fits', 'chandra_7ms/CDFS-7Ms-0p5to2-bin1-01.emap.fits', 'chandra_7ms/CDFS-7Ms-0p5to7-bin1-01.emap.fits', 'chandra_7ms/CDFS-7Ms-2to7-bin1-01.emap.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f105w_v1.0_wht.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_wht.fits', 'candels/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_wht.fits']
        o = ['chandra_4ms/cutout_chandra4_05_to_2_', 'chandra_4ms/cutout_chandra4_05_to_8_', 'chandra_4ms/cutout_chandra4_2_to_8_', 'chandra_7ms/cutout_chandra7_05_to_2_', 'chandra_7ms/cutout_chandra7_05_to_7_', 'chandra_7ms/cutout_chandra7_2_to_7_', 'candels/cutout_candels_105_', 'candels/cutout_candels_125_', 'candels/cutout_candels_160_']
        v = []
        for l in range(len(p)):
            path = '/Users/justin/Documents/Master_Thesis/data/surveys/' + p[l]
            q.append(fits.open(path)[0])
            v.append(wcs.WCS(q[l].header))
            for i in range(cat7ms.size):
                ra = cat7ms.field('RA')[i]
                dec = cat7ms.field('DEC')[i]
                center = coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
                co = Cutout2D(q[l].data, center, size=[arcsec,arcsec]*u.arcsec, wcs=v[l])
                hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
                if l > 5:
                    hdu.header['CD1_1'] = -1.666667E-5 
                    hdu.header['CD1_2'] = 0.
                    hdu.header['CD2_1'] = 0.
                    hdu.header['CD2_2'] = 1.666667E-5 
                if catalogue == 'main':
                    if Weight==False:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + str(i+1) + '.fits'
                    else:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + str(i+1) + 'wht.fits'
                if catalogue == 'secondary':
                    if Weight==False:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + 'x' + str(i+1) + '.fits'
                    else:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + 'x' + str(i+1) + 'wht.fits'
                if catalogue == 'extra':
                    if Weight==False:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + 'n' + str(i+1) + '.fits'
                    else:
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + o[l] + 'n' + str(i+1) + 'wht.fits'
                
                os.system('rm -f ' + outpath)
                hdu.writeto(outpath)

    #goods
    if (Instr=='All') | (Instr=='GOODS'):
        a = [11, 12, 13, 14, 21, 22, 23, 24, 25, 31, 32, 33, 34, 35, 42, 43, 44, 45]
        b = ['b', 'i', 'v', 'z']
        im = []
        w = []
        for k in range(len(b)):
            del im[:]
            del w[:]
            for i in range(len(a)):
                if Weight==False:
                    path = '/Users/justin/Documents/Master_Thesis/data/surveys/goods/h_s' + b[k] + '_sect' + str(a[i]) + '_v2.0_drz_img.fits'
                else:
                    path = '/Users/justin/Documents/Master_Thesis/data/surveys/goods/h_s' + b[k] + '_sect' + str(a[i]) + '_v2.0_wht_img.fits'
                im.append(fits.open(path)[0])
                w.append(wcs.WCS(im[i].header))

            for i in range(cat7ms.size):
                ra = cat7ms.field('RA')[i]
                dec = cat7ms.field('DEC')[i]
                for j in range(len(im)):
                    x, y = w[j].all_world2pix(ra, dec, 1)
                    if (x > 0.) & (x < 8192.) & (y > 0.) & (y < 8192.) & (i!=542):
                        center = coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
                        co = Cutout2D(im[j].data, center, size=[arcsec,arcsec]*u.arcsec, wcs=w[j], mode='trim')
                        hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
                        if (hdu.header['NAXIS1']!=pix_goods) | (hdu.header['NAXIS2']!=pix_goods):
                            newx = (2.*hdu.header['NAXIS1'])-(pix_goods)
                            newy = (2.*hdu.header['NAXIS2'])-(pix_goods)
                            co = Cutout2D(im[j].data, center, size=[newy,newx], wcs=w[j], mode='trim')
                            hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
                        hdu.header['CD1_1'] = -8.3333334575000E-6
                        hdu.header['CD1_2'] = 0.
                        hdu.header['CD2_1'] = 0.
                        hdu.header['CD2_2'] = 8.3333334575000E-6
                        #hdu.header.remove('PC1_1')
                        #hdu.header.remove('PC2_2')
                        if catalogue == 'main':
                            if Weight==False:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + str(i+1) + '.fits'
                            else:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + str(i+1) + 'wht.fits'
                        if catalogue == 'secondary':
                            if Weight==False:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + 'x' + str(i+1) + '.fits'
                            else:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + 'x' + str(i+1) + 'wht.fits'
                        if catalogue == 'extra':
                            if Weight==False:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + 'n' + str(i+1) + '.fits'
                            else:
                                outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/goods/cutout_goods_' + b[k] + 'n' + str(i+1) + 'wht.fits'
                        os.system('rm -f ' + outpath)
                        hdu.writeto(outpath)

    #hugs
    if (Instr=='All') | (Instr=='HUGS'):
        if (Weight==False) & (RMS==False):
            path = '/Users/justin/Documents/Master_Thesis/data/surveys/hugs/HUGS_GOODSS_K_naturalseeing/HUGS_GOODSS_K_naturalseeing.fits'
        if (Weight==True) & (RMS==False):
            path = '/Users/justin/Documents/Master_Thesis/data/surveys/hugs/HUGS_GOODSS_K_naturalseeing/HUGS_GOODSS_K_naturalseeing.wht.fits'
        if (Weight==False) & (RMS==True):
            path = '/Users/justin/Documents/Master_Thesis/data/surveys/hugs/HUGS_GOODSS_K_naturalseeing/HUGS_GOODSS_K_naturalseeing.rms.fits'    
        ima = fits.open(path)[0]
        wc = wcs.WCS(ima.header)
        for i in range(cat7ms.size):
            ra = cat7ms.field('RA')[i]
            dec = cat7ms.field('DEC')[i]
            x, y = wc.all_world2pix(ra, dec, 1)
            if (x > 0.) & (x < 9373.) & (y > 0.) & (y < 10930.):
                center = coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
                co = Cutout2D(ima.data, center, size=[arcsec,arcsec]*u.arcsec, wcs=wc)
                hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
                hdu.header['CD1_1'] = -2.959764788102E-05
                hdu.header['CD1_2'] = 0.
                hdu.header['CD2_1'] = 0.
                hdu.header['CD2_2'] = 2.959764788102E-05
                if catalogue == 'main':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + str(i+1) + 'rms.fits'
                if catalogue == 'secondary':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'x' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'x' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'x' + str(i+1) + 'rms.fits'
                if catalogue == 'extra':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'n' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'n' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/hugs/cutout_hugs_' + 'n' + str(i+1) + 'rms.fits'
                os.system('rm -f ' + outpath)
                hdu.writeto(outpath)
    #spitzer
    if (Instr=='All') | (Instr=='Spitzer'):

        if (Weight==False) & (RMS==False):
            path = ['/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch1.fits', '/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch2.fits']
        if (Weight==True) & (RMS==False):
            path = ['/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch1_exp.fits', '/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch2_exp.fits']
        if (Weight==False) & (RMS==True):
            path = ['/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch1_rms.fits', '/Users/justin/Documents/Master_Thesis/data/surveys/spitzer/SIMPLE_ECDFS_IRAC_ch2_rms.fits']
        for j in range(len(path)):
            ima = fits.open(path[j])[0]
            for i in range(cat7ms.size):
                ra = cat7ms.field('RA')[i]
                dec = cat7ms.field('DEC')[i]
                wc = wcs.WCS(ima.header)
                center = coordinates.SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
                co = Cutout2D(ima.data, center, size=[arcsec,arcsec]*u.arcsec, wcs=wc)
                hdu = fits.PrimaryHDU(data=co.data, header=co.wcs.to_header())
                hdu.header['CD1_1'] = -1.6666666666700E-4 
                hdu.header['CD1_2'] = 0.
                hdu.header['CD2_1'] = 0.
                hdu.header['CD2_2'] = 1.66666666667000E-4 
                if catalogue == 'main':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + str(i+1) + 'rms.fits'
                if catalogue == 'secondary':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'x' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'x' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'x' + str(i+1) + 'rms.fits'
                if catalogue == 'extra':
                    if (Weight==False) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'n' + str(i+1) + '.fits'
                    if (Weight==True) & (RMS==False):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'n' + str(i+1) + 'wht.fits'
                    if (Weight==False) & (RMS==True):
                        outpath = '/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/spitzer/cutout_spitzer_ch' + str(j+1) + '_' + 'n' + str(i+1) + 'rms.fits'
                os.system('rm -f ' + outpath)
                hdu.writeto(outpath)
    return 0
    
def make4msregion():
    path4ms = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_4ms/4MS_Chandra.fits'
    cat4 = fits.open(path4ms)[1].data
    ID = cat4.field('col1')
    RAh = cat4.field('col2')
    RAmin = cat4.field('col3')
    RAs = cat4.field('col4')
    DECdeg = cat4.field('col5')
    DECarcmin = cat4.field('col6')
    DECarcsec = cat4.field('col7')
    sigmax = cat4.field('col10')
    RA_deg = ra2deg(RAh, RAmin, RAs)
    DEC_deg = dec2deg(DECdeg/np.absolute(DECdeg), np.absolute(DECdeg), DECarcmin, DECarcsec)
    sigmax_deg = sigmax/3600.
    file4 = open('/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_4ms/region_file.fits', 'w')
    file4.write("global color=red\n")
    for i in range(ID.size):
        #print('FK5;circle(' + str(RA_deg[i]) + ',' + str(DEC_deg[i]) + ',' + str(sigmax_deg[i]) + ');')
        file4.write("FK5;circle(%f,%f,%f)\n" % (RA_deg[i], DEC_deg[i], sigmax_deg[i]))
        file4.write("FK5;text %f %f {                     ID %s} # color=yellow, font={helvetica 8 normal roman}\n" % (RA_deg[i], DEC_deg[i], ID[i]))
    file4.close()
    return 0

def make7msregion(text=True, out='region_file'):
    path7ms = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl04.fits'
    cat7 = fits.open(path7ms)[1].data
    ID = cat7.field('ID')
    RA = cat7.field('RA')
    DEC = cat7.field('DEC')
    SIGMAX = cat7.field('SIGMAX')
    SIGMAX_deg = SIGMAX/3600.
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/' + out + '.fits'
    file7 = open(path, 'w')
    file7.write("global color=blue\n")
    for i in range(RA.size):
        #print('FK5;circle(' + str(RA[i]) + ',' + str(DEC[i]) + ',' + str(SIGMAX_deg[i]) + ');')
        #file7.write("FK5;circle(%f,%f,%f)\n" % (RA[i], DEC[i], SIGMAX_deg[i]))
        file7.write("FK5;point(%f,%f) # point=circle 2 \n" % (RA[i], DEC[i]))
        if text == True:
            file7.write("FK5;text %f %f {                     ID %s} # color=yellow, font={helvetica 8 normal roman}\n" % (RA[i], DEC[i], ID[i]))
    file7.close()
    return 0

def cutoutpdf(sources='all', spitzer=False, aperture=False,  outpath='/Users/justin/Documents/Master_Thesis/presentation/plots/cutout_stamps/classic/cutout_row_', filetype='.pdf'):
    path2 = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl05.fits'
    cat2 = fits.open(path2)[1].data
    if sources == 'extra':
        path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/short_table.fits'
    else:
        path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    cat = fits.open(path)[1].data
    if sources == 'secondary':
        ID = cat2.field('ID')
        RA = cat2.field('RA')
        DEC = cat2.field('DEC')
        PosErr =cat2.field('SIGMAX')
        hugs = hugs_sect(RA, DEC)
        goods = goods_sect(RA, DEC)
    elif sources == 'extra':
        ID = cat.field('ID')
        RA = cat.field('RA')
        DEC = cat.field('DEC')
        PosErr = np.zeros(len(cat))
        hugs = hugs_sect(RA, DEC)
        goods = goods_sect(RA, DEC)
    else:
        ID = cat.field('ID')
        RA = cat.field('RA')
        DEC = cat.field('DEC')
        PosErr =cat.field('SIGMAX')
        hugs = cat.field('HUGS_sect')
        goods = cat.field('GOODS_sect')
    check = []
    typ = []
    title = []
    if spitzer == True:
        check = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5]
        typ = ['chandra_7ms/cutout_chandra7_2_to_7_', 'chandra_7ms/cutout_chandra7_05_to_2_', 'goods/cutout_goods_b', 'goods/cutout_goods_v', 'goods/cutout_goods_i', 'goods/cutout_goods_z', 'candels/cutout_candels_105_', 'candels/cutout_candels_125_', 'candels/cutout_candels_160_', 'hugs/cutout_hugs_', 'spitzer/cutout_spitzer_ch1_', 'spitzer/cutout_spitzer_ch2_']
        title = ['Hard X-ray', 'Soft X-ray', 'B', 'V', 'i', 'z', 'Y', 'J', 'H', 'K', '3.6 micron', '4.5 micron']
    else:
        check = [1, 1, 2, 2, 2, 2, 3, 3, 3, 4]
        typ = ['chandra_7ms/cutout_chandra7_2_to_7_', 'chandra_7ms/cutout_chandra7_05_to_2_', 'goods/cutout_goods_b', 'goods/cutout_goods_v', 'goods/cutout_goods_i', 'goods/cutout_goods_z', 'candels/cutout_candels_105_', 'candels/cutout_candels_125_', 'candels/cutout_candels_160_', 'hugs/cutout_hugs_']
        title = ['Hard X-ray', 'Soft X-ray', 'B', 'V', 'i', 'z', 'Y', 'J', 'H', 'K']

    line = []
    path = []

    for i in range(len(ID)):
        if (sources == 'extra') | (sources == 'secondary') | (sources == 'all') | (sources == i) | ((sources == 'main') & (cat.field('Visual_Flag')[i] > 3)):
            num_plots = 10
            if spitzer == True:
                num_plots += 2
            fig = plt.figure(facecolor='w', figsize=(num_plots*1.5,2.), dpi=600)
            #fig = plt.figure(facecolor='w', figsize=(4.5,8), dpi=600)
            gs=gridspec.GridSpec(1,num_plots)
            #gs=gridspec.GridSpec(4,3)
            del path[:]
            for k in range(0,num_plots):
                print(k)
                if sources == 'secondary':
                    path.append('/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + typ[k] + 'x' + str(i+1) + '.fits')
                elif sources == 'extra':
                    path.append('/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + typ[k] + 'n' + str(i+1) + '.fits')
                else: 
                    path.append('/Users/justin/Documents/Master_Thesis/data/surveys/cutouts/' + typ[k] + str(i+1) + '.fits')
                #line.append(fits.open(path[j])[0])
                ax=plt.subplot(gs[0,k])
                #ax=plt.subplot(gs[int(np.floor(k/3.)),k%3])
                ax.set_aspect('equal')
                ax.axes.get_xaxis().set_ticks([]); ax.axes.get_yaxis().set_ticks([])
                plt.subplots_adjust(left=0.01, bottom=0.13, right=0.99, top=0.87, wspace=0, hspace=0)
                ax.set_title(title[k], size=12, fontname='serif')
                if ((check[k] == 2) & ((goods[i] == 0) | (i == 542))) | ((check[k] == 4) & (hugs[i] == 0)):  
                    print(str(i+1))
                else:
                    f1=aplpy.FITSFigure(path[k], figure=fig, subplot=list(ax.get_position().bounds))
                    if check[k]==1:
                        f1.show_grayscale(invert=True, stretch='linear', vmin=0., smooth=1., interpolation='nearest')
                    else:
                        f1.show_grayscale(invert=True, stretch='linear', vmin=0., interpolation='nearest')
                    #f1.add_label(0.5, 0.9, title[k], relative=True, size='large', style='normal', family='serif')
                    if k == 0:
                        f1.add_scalebar(1./3600.)
                        f1.scalebar.show(1./3600.)
                        f1.scalebar.set_corner('bottom left')
                        f1.scalebar.set_linewidth(3)
                        f1.scalebar.set_color('black')
                        f1.scalebar.set_font(size='medium', style='normal')
                        f1.scalebar.set_label('1 arcsec')
                        if sources == 'secondary':
                            f1.add_label(0.2, 0.9, 'X' + str(ID[i]), relative=True, size='large', style='normal', weight='bold', family='serif')
                        elif sources == 'extra':
                            f1.add_label(0.2, 0.9, 'N' + str(ID[i]), relative=True, size='large', style='normal', weight='bold', family='serif')
                        else:
                            #pass
                            f1.add_label(0.2, 0.9, str(ID[i]), relative=True, size='large', style='normal', weight='bold', family='serif')
                    f1.axis_labels.hide_x(); f1.axis_labels.hide_y()
                    f1.tick_labels.hide_x(); f1.tick_labels.hide_y()
                    f1.ticks.show_x(); f1.ticks.show_y()
                    if aperture == True:
                        lw=0.7
                        if (k == 2):
                            if cat.field('GOODS_b_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('GOODS_b_RA')[i], cat.field('GOODS_b_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('GOODS_b_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('GOODS_b_RA')[i], cat.field('GOODS_b_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('GOODS_b_FluxAper')[i] == 0) & (cat.field('GOODS_b_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 3):
                            if cat.field('GOODS_v_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('GOODS_v_RA')[i], cat.field('GOODS_v_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('GOODS_v_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('GOODS_v_RA')[i], cat.field('GOODS_v_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('GOODS_v_FluxAper')[i] == 0) & (cat.field('GOODS_v_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 4):
                            if cat.field('GOODS_i_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('GOODS_i_RA')[i], cat.field('GOODS_i_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('GOODS_i_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('GOODS_i_RA')[i], cat.field('GOODS_i_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('GOODS_i_FluxAper')[i] == 0) & (cat.field('GOODS_i_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 5):
                            if cat.field('GOODS_z_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('GOODS_z_RA')[i], cat.field('GOODS_z_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('GOODS_z_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('GOODS_z_RA')[i], cat.field('GOODS_z_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('GOODS_z_FluxAper')[i] == 0) & (cat.field('GOODS_z_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 6):
                            if cat.field('CANDELS_105_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('CANDELS_105_RA')[i], cat.field('CANDELS_105_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('CANDELS_105_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('CANDELS_105_RA')[i], cat.field('CANDELS_105_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('CANDELS_105_FluxAper')[i] == 0) & (cat.field('CANDELS_105_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 7):
                            if cat.field('CANDELS_125_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('CANDELS_125_RA')[i], cat.field('CANDELS_125_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('CANDELS_125_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('CANDELS_125_RA')[i], cat.field('CANDELS_125_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('CANDELS_125_FluxAper')[i] == 0) & (cat.field('CANDELS_125_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 8):
                            if cat.field('CANDELS_160_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('CANDELS_160_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('CANDELS_160_FluxAper')[i] == 0) & (cat.field('CANDELS_160_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 9):
                            if cat.field('HUGS_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('HUGS_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('HUGS_FluxAper')[i] == 0) & (cat.field('HUGS_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Aperture')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 10):
                            if cat.field('Spitzer_ch1_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('Spitzer_ch1_RA')[i], cat.field('Spitzer_ch1_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('Spitzer_ch1_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('Spitzer_ch1_RA')[i], cat.field('Spitzer_ch1_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('Spitzer_ch1_FluxAper')[i] == 0) & (cat.field('Spitzer_ch1_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    if cat.field('Spitzer_Aper_Flag')[i] > 10.:
                                        f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], (cat.field('Spitzer_Aper_Flag')[i]-10.)/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                    else:
                                        f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    if cat.field('Spitzer_Aper_Flag')[i] > 10.:
                                        f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], (cat.field('Spitzer_Aper_Flag')[i]-10.)/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                    else:
                                        f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                        if (k == 11):
                            if cat.field('Spitzer_ch2_FluxAper')[i] < 0:
                                f1.show_circles(cat.field('Spitzer_ch2_RA')[i], cat.field('Spitzer_ch2_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='blue', linewidth=lw, linestyle='solid')
                            if cat.field('Spitzer_ch2_FluxAper')[i] > 0:
                                f1.show_circles(cat.field('Spitzer_ch2_RA')[i], cat.field('Spitzer_ch2_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='solid')
                            if (cat.field('Spitzer_ch2_FluxAper')[i] == 0) & (cat.field('Spitzer_ch2_FluxErrAper')[i] != 0):
                                if (i==100) | (i==392) | (i==452) | (i==545):
                                    if cat.field('Spitzer_Aper_Flag')[i] > 10.:
                                        f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], (cat.field('Spitzer_Aper_Flag')[i]-10.)/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                    else:
                                        f1.show_circles(cat.field('HUGS_RA')[i], cat.field('HUGS_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                else:
                                    if cat.field('Spitzer_Aper_Flag')[i] > 10.:
                                        f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], (cat.field('Spitzer_Aper_Flag')[i]-10.)/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                                    else:
                                        f1.show_circles(cat.field('CANDELS_160_RA')[i], cat.field('CANDELS_160_Dec')[i], cat.field('Spitzer_Aper_Flag')[i]/7200., edgecolor='red', linewidth=lw, linestyle='dashed')
                    if (sources != 'secondary') & (sources != 'extra'):
                        if ((k == 2) & (cat.field('GOODS_b_Usable')[i]!=0)) | ((k == 3) & (cat.field('GOODS_v_Usable')[i]!=0)) | ((k == 4) & (cat.field('GOODS_i_Usable')[i]!=0)) | ((k == 5) & (cat.field('GOODS_z_Usable')[i]!=0)) | ((k == 6) & (cat.field('CANDELS_105_Usable')[i]!=0)) | ((k == 7) & (cat.field('CANDELS_125_Usable')[i]!=0)) | ((k == 8) & (cat.field('CANDELS_160_Usable')[i]!=0)) | ((k == 9) & (cat.field('HUGS_Usable')[i]!=0)):
                            f1.show_circles(RA, DEC, PosErr/3600., edgecolor='black', linewidth=0.5)
                        if (k < 2) | (k > 9):
                            f1.show_circles(RA, DEC, PosErr/3600., edgecolor='black', linewidth=0.5)
                    else:
                        f1.show_circles(RA, DEC, PosErr/3600., edgecolor='black', linewidth=0.5)
                
                    #f1.set_tick_labels_font(size='x-small')
                    #f1.set_axis_labels_font(size='small')
                    #f1.hide_yaxis_label()
                    #f1.hide_ytick_labels()
            if sources == 'secondary':
                out = outpath + 'x' + str(i+1) + filetype
            elif sources == 'extra':
                out = outpath + 'n' + str(i+1) + filetype
            else:
                out = outpath + str(i+1) + filetype
            fig.savefig(out, bbox_inches='tight', dpi=300)
            plt.clf()
    
    return 0

def mergecatalogues(criterion='closest', pixel_width=20., arcsec_radius=1.):
    path2 = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    cat = fits.open(path2)[1].data
    ra = cat.field('RA')
    dec = cat.field('DEC')
    misx = cat.field('GOODS_x_pixel')
    misy = cat.field('GOODS_y_pixel')
    a = ['b', 'v', 'i', 'z']
    t2 = []
    for j in range(len(a)):
        for i in range(1008):
            path = '/Users/justin/Documents/Master_Thesis/data/catalogues/cutouts/goods/cat_goods_' + a[j] + str(i+1) + '.cat'
            k=0
            if (misx[i] > pixel_width) & (misy[i] > pixel_width):
                try:
                    t1 = pd.read_csv(path, delim_whitespace=True, header=None, names=['#', 'x', 'y', 'RA', 'Dec', 'Mag', 'MagErr', 'Flux', 'FluxErr'])
                    #t1 = fits.open(path)[1].data
                    radius = ((t1['RA']-ra[i]).pow(2) + (t1['Dec']-dec[i]).pow(2)).pow(0.5)
                    n = -1
                    for l in range(radius.size):
                        if radius[l] < arcsec_radius/3600.:
                            if (n == -1):
                                n = l
                            else:
                                if criterion == 'closest':
                                    if (radius[l] < radius[n]): #closest
                                        n = l
                                if criterion == 'brightest':
                                    if (t1['Mag'][l] < t1['Mag'][n]): #brightest
                                        n = l
                    if n != -1:
                        t1 = t1.iloc[[n],:]
                        t1['ID'] = i+1
                    else:
                        error_trigger
                    k=1
                except:
                    print('file ' + str(i+1) + ' has failed opening or has no source')
                    pass
            if k==1:
                try:
                    t2[j] = t2[j].append(t1, ignore_index=True, verify_integrity=False)
                except:
                    t2.append(t1)
                    pass
        if criterion == 'closest':
            path = '/Users/justin/Documents/Master_Thesis/data/catalogues/cutouts/goods/cat_goods_' + a[j] + 'full_closest.cat'
        if criterion == 'brightest':
            path = '/Users/justin/Documents/Master_Thesis/data/catalogues/cutouts/goods/cat_goods_' + a[j] + 'full_brightest.cat'
        t2[j].to_csv(path_or_buf=path, sep=' ')
    return 0

def magmagplot(criterion='closest', band='b', separation='1', unit='magnitude'):
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/crossmatchs/goods/crossmatch_goods_' + band + separation + '_' + criterion + '.fits'
    cat1 = fits.open(path)[1].data
    if unit == 'magnitude':
        mag1a = cat1.field('col7_1') #7 mag, 8 magerr, 9 flux, 10 fluxerr
        err1a = cat1.field('col8_1')
        mag1b = cat1.field('col42') #42 mag(auto), 43 magerr, 44 flux, 45 fluxerr
        err1b = cat1.field('col43')
    if unit == 'flux':
        mag1a = cat1.field('col9_1') #7 mag, 8 magerr, 9 flux, 10 fluxerr
        err1a = cat1.field('col10_1')
        mag1b = cat1.field('col44') #42 mag(auto), 43 magerr, 44 flux, 45 fluxerr
        err1b = cat1.field('col45')
    f = plt.figure()
    f.subplots_adjust(wspace=0.35)
    ax1 = f.add_subplot(121)
    plt.errorbar(mag1a, mag1b, err1a, err1b, ls='', capsize=0.3, capthick=0.3, elinewidth=0.3)
    ax1.plot(mag1a, mag1b, marker='.', ms=1, ls='', color='red')
    x0,x1 = ax1.get_xlim()
    y0,y1 = ax1.get_ylim()
    z0 = min(x0, y0)
    z1 = max(x1, y1)
    ax1.plot([z0, z1], [z0, z1], color='black', marker='', ls='-', lw=0.3)
    #ax1.set_aspect('equal')
    plt.axis([z0, z1, z0, z1])
    ax1.set_aspect((z1-z0)/(z1-z0))
    if unit == 'magnitude':
        plt.ylabel('Mag$_\mathrm{off}$')
        plt.xlabel('Mag')
    if unit == 'flux':
        plt.ylabel('Flux$_\mathrm{off}$')
        plt.xlabel('Flux')
    ax2 = f.add_subplot(122)
    plt.errorbar(mag1a, mag1b-mag1a, err1a, err1b, ls='', capsize=0.3, capthick=0.3, elinewidth=0.3)
    ax2.plot(mag1a, mag1b-mag1a, marker='.', ms=1, ls='', color='red')
    x0,x1 = ax2.get_xlim()
    y0,y1 = ax2.get_ylim()
    ax2.plot([x0, x1], [0, 0], color='black', marker='', ls='-', lw=0.3)
    ax2.set_aspect((x1-x0)/(y1-y0))
    if unit == 'magnitude':
        plt.ylabel('Mag$_\mathrm{off}$ - Mag')
        plt.xlabel('Mag')
    if unit == 'flux':
        plt.ylabel('Flux$_\mathrm{off}$ - Flux')
        plt.xlabel('Flux')
    if unit == 'magnitude':
        path = '/Users/justin/Documents/Master_Thesis/presentation/plots/goods_mag/goods_' + band + separation + '_' + criterion + '_mag_mag.eps'
    if unit == 'flux':
        path = '/Users/justin/Documents/Master_Thesis/presentation/plots/goods_mag/goods_' + band + separation + '_' + criterion + '_flux_flux.eps'
    plt.savefig(path, format='eps', dpi=1000)
    return 0

def pos_offset_plot():
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    cat = fits.open(path)[1].data
    chandra_RA = cat.field('RA')
    chandra_Dec = cat.field('Dec')
    goods_b_RA = cat.field('GOODS_b_RA')
    goods_b_Dec = cat.field('GOODS_b_Dec')
    goods_v_RA = cat.field('GOODS_v_RA')
    goods_v_Dec = cat.field('GOODS_v_Dec')
    goods_i_RA = cat.field('GOODS_i_RA')
    goods_i_Dec = cat.field('GOODS_i_Dec')
    goods_z_RA = cat.field('GOODS_z_RA')
    goods_z_Dec = cat.field('GOODS_z_Dec')
    b_RA_offset = 0.
    v_RA_offset = 0.
    i_RA_offset = 0.
    z_RA_offset = 0.
    ref_RA_offset = 0.
    b_Dec_offset = 0.
    v_Dec_offset = 0.
    i_Dec_offset = 0.
    z_Dec_offset = 0.
    ref_Dec_offset = 0.
    ref_RA_var = 0.
    ref_Dec_var = 0.
    n_b=0.
    n_v=0.
    n_i=0.
    n_z=0.
    n_ref=0.
    for i in range(len(cat)):
        if cat.field('GOODS_b_Flag')[i] ==1:
            b_RA_offset += (chandra_RA[i]-goods_b_RA[i])*3600.
            b_Dec_offset += (chandra_Dec[i]-goods_b_Dec[i])*3600.
            n_b += 1.
        if cat.field('GOODS_v_Flag')[i] ==1:
            v_RA_offset += (chandra_RA[i]-goods_v_RA[i])*3600.
            v_Dec_offset += (chandra_Dec[i]-goods_v_Dec[i])*3600.
            n_v += 1.
        if cat.field('GOODS_i_Flag')[i] ==1:
            i_RA_offset += (chandra_RA[i]-goods_i_RA[i])*3600.
            i_Dec_offset += (chandra_Dec[i]-goods_i_Dec[i])*3600.
            n_i += 1.
        if cat.field('GOODS_z_Flag')[i] ==1:
            z_RA_offset += (chandra_RA[i]-goods_z_RA[i])*3600.
            z_Dec_offset += (chandra_Dec[i]-goods_z_Dec[i])*3600.
            n_z += 1.
        if cat.field('Visual_Flag')[i] > 3:
            if (i==100) | (i==392) | (i==452) | (i==545):
                ref_RA_offset += (chandra_RA[i]-cat.field('HUGS_RA')[i])*3600.
                ref_Dec_offset += (chandra_Dec[i]-cat.field('HUGS_Dec')[i])*3600.
                n_ref += 1.
            else:
                ref_RA_offset += (chandra_RA[i]-cat.field('CANDELS_160_RA')[i])*3600.
                ref_Dec_offset += (chandra_Dec[i]-cat.field('CANDELS_160_Dec')[i])*3600.
                n_ref += 1.
    b_RA_offset /= n_b
    v_RA_offset /= n_v
    i_RA_offset /= n_i
    z_RA_offset /= n_z
    ref_RA_offset /= n_ref
    b_Dec_offset /= n_b
    v_Dec_offset /= n_v
    i_Dec_offset /= n_i
    z_Dec_offset /= n_z
    ref_Dec_offset /= n_ref
    n_ref=0.
    for i in range(len(cat)):
        if cat.field('Visual_Flag')[i] > 3:
            if (i==100) | (i==392) | (i==452) | (i==545):
                ref_RA_var += (((chandra_RA[i]-cat.field('HUGS_RA')[i])*3600.)-ref_RA_offset)**2.
                ref_Dec_var += (((chandra_Dec[i]-cat.field('HUGS_Dec')[i])*3600.)-ref_Dec_offset)**2.
                n_ref += 1.
            else:
                ref_RA_var += (((chandra_RA[i]-cat.field('CANDELS_160_RA')[i])*3600.)-ref_RA_offset)**2.
                ref_Dec_var += (((chandra_Dec[i]-cat.field('CANDELS_160_Dec')[i])*3600.)-ref_Dec_offset)**2.
                n_ref += 1.
    ref_RA_var = (ref_RA_var**0.5)/n_ref
    ref_Dec_var = (ref_Dec_var**0.5)/n_ref
    print(ref_RA_offset, ref_Dec_offset, ref_RA_var, ref_Dec_var)
    
    f1 = plt.figure()
    ax1 = f1.add_subplot(111)
    ax1.plot((chandra_RA-goods_b_RA)*3600., (chandra_Dec-goods_b_Dec)*3600., marker='.', ms=1.5, ls='', color='red')
    plt.axis([-1, 1, -1, 1])
    ax1.set_aspect('equal')
    ax1.plot([-1, 1], [0, 0], color='black', marker='', ls='--', lw=0.3)
    ax1.plot([0, 0], [-1, 1], color='black', marker='', ls='--', lw=0.3)
    plt.xlabel('RA$_\mathrm{Chandra}$-RA$_\mathrm{ACS}$ in arcsec')
    plt.ylabel('Dec$_\mathrm{Chandra}$-Dec$_\mathrm{ACS}$ in arcsec')
    ax1.plot(b_RA_offset, b_Dec_offset, color='white', marker='*', ls='', ms=10)
    #plt.title('GOODS b filter')
    path = '/Users/justin/Documents/Master_Thesis/presentation/plots/pos_offset/RA_Dec_offset_goods_b.eps'
    plt.savefig(path, format='eps', dpi=1000, bbox_inches='tight')
    
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)
    ax2.plot((chandra_RA-goods_v_RA)*3600., (chandra_Dec-goods_v_Dec)*3600., marker='.', ms=1.5, ls='', color='red')
    plt.axis([-1, 1, -1, 1])
    ax2.set_aspect('equal')
    ax2.plot([-1, 1], [0, 0], color='black', marker='', ls='--', lw=0.3)
    ax2.plot([0, 0], [-1, 1], color='black', marker='', ls='--', lw=0.3)
    plt.xlabel('RA$_\mathrm{Chandra}$-RA$_\mathrm{ACS}$ in arcsec')
    plt.ylabel('Dec$_\mathrm{Chandra}$-Dec$_\mathrm{ACS}$ in arcsec')
    ax2.plot(v_RA_offset, v_Dec_offset, color='white', marker='*', ls='', ms=10)
    #plt.title('GOODS v filter')
    path = '/Users/justin/Documents/Master_Thesis/presentation/plots/pos_offset/RA_Dec_offset_goods_v.eps'
    plt.savefig(path, format='eps', dpi=1000, bbox_inches='tight')
    
    f3 = plt.figure()
    ax3 = f3.add_subplot(111)
    ax3.plot((chandra_RA-goods_i_RA)*3600., (chandra_Dec-goods_i_Dec)*3600., marker='.', ms=1.5, ls='', color='red')
    plt.axis([-1, 1, -1, 1])
    ax3.set_aspect('equal')
    ax3.plot([-1, 1], [0, 0], color='black', marker='', ls='--', lw=0.3)
    ax3.plot([0, 0], [-1, 1], color='black', marker='', ls='--', lw=0.3)
    plt.xlabel('RA$_\mathrm{Chandra}$-RA$_\mathrm{ACS}$ in arcsec')
    plt.ylabel('Dec$_\mathrm{Chandra}$-Dec$_\mathrm{ACS}$ in arcsec')
    ax3.plot(i_RA_offset, i_Dec_offset, color='white', marker='*', ls='', ms=10)
    #plt.title('GOODS i filter')
    path = '/Users/justin/Documents/Master_Thesis/presentation/plots/pos_offset/RA_Dec_offset_goods_i.eps'
    plt.savefig(path, format='eps', dpi=1000, bbox_inches='tight')
    
    f4 = plt.figure()
    ax4 = f4.add_subplot(111)
    ax4.plot((chandra_RA-goods_z_RA)*3600., (chandra_Dec-goods_z_Dec)*3600., marker='.', ms=1.5, ls='', color='red')
    plt.axis([-1, 1, -1, 1])
    ax4.set_aspect('equal')
    ax4.plot([-1, 1], [0, 0], color='black', marker='', ls='--', lw=0.3)
    ax4.plot([0, 0], [-1, 1], color='black', marker='', ls='--', lw=0.3)
    plt.xlabel('RA$_\mathrm{Chandra}$-RA$_\mathrm{ACS}$ in arcsec')
    plt.ylabel('Dec$_\mathrm{Chandra}$-Dec$_\mathrm{ACS}$ in arcsec')
    ax4.plot(z_RA_offset, z_Dec_offset, color='white', marker='*', ls='', ms=10)
    #plt.title('GOODS z filter')
    path = '/Users/justin/Documents/Master_Thesis/presentation/plots/pos_offset/RA_Dec_offset_goods_z.eps'
    plt.savefig(path, format='eps', dpi=1000, bbox_inches='tight')
    
    ref_ra=[]
    ref_dec=[]
    for i in range(len(cat)):
        if cat.field('Visual_Flag')[i] > 3:
            if (i==100) | (i==392) | (i==452) | (i==545):
                ref_ra.append((cat.field('RA')[i]-cat.field('HUGS_RA')[i])*3600.)
                ref_dec.append((cat.field('Dec')[i]-cat.field('HUGS_Dec')[i])*3600.)
            else:
                ref_ra.append((cat.field('RA')[i]-cat.field('CANDELS_160_RA')[i])*3600.)
                ref_dec.append((cat.field('Dec')[i]-cat.field('CANDELS_160_Dec')[i])*3600.)
    f5 = plt.figure()
    ax5 = f5.add_subplot(111)
    ax5.plot(ref_ra, ref_dec, marker='.', ms=3.5, ls='', color='black')
    plt.axis([-1, 1, -1, 1])
    ax5.set_aspect('equal')
    ax5.plot([-1, 1], [0, 0], color='black', marker='', ls='--', lw=0.3)
    ax5.plot([0, 0], [-1, 1], color='black', marker='', ls='--', lw=0.3)
    plt.xlabel('RA$_\mathrm{Chandra}$-RA$_\mathrm{Ref}$ in arcsec')
    plt.ylabel('Dec$_\mathrm{Chandra}$-Dec$_\mathrm{Ref}$ in arcsec')
    ax5.plot(ref_RA_offset, ref_Dec_offset, color='white', marker='*', ls='', ms=12.5)
    #plt.title('GOODS z filter')
    path = '/Users/justin/Documents/Master_Thesis/presentation/plots/pos_offset/RA_Dec_offset_candels_H.eps'
    plt.savefig(path, format='eps', dpi=1000, bbox_inches='tight')
    return 0

def get_filter_wavelengths():
    #BVizYJHK3.6/4.5, pivot, min, max, Weff, FWHM
    a = [[4332.66,3599.05,4861.28,760.18,885.20],[5963.10,4632.17,7179.45,1872.25,2281.53],[7689.71,6801.39,8630.66,1299.16,1515.12],[9177.76,8015.10,10946.45,1552.04,1531.72],[10550.29,8946.81,12129.45,2648.66,2917.03],[12486.05,10844.59,14139.24,2845.23,3005.20],[15370.33,13854.10,16999.35,2682.29,2874.18],[21460,0,0,0,3240],[35508,0,0,0,7432],[44960,0,0,0,10097]]
    return a

def get_sextractor_params():
    #write here all sextractor parameters available for candels, hugs or even spitzer data, in case it's useful later
    #goods parameters are already in the .txt files
    #pixel scale: 0"03 goods, 0"06 candels
    #psf fwhm: goods v 0"08, z 0"09 candels 0"11, 0"12, 0"18 hugs 0"4
    #zero point: candels 26.2687, 26.2303, 25.9463
    return 0

def sed_plot():
    a = get_filter_wavelengths()
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    cat = fits.open(path)[1].data
    f = plt.figure()
    for i in range(len(cat)):
        if (cat.field('GOODS_b_Flag')[i] != 0) | (cat.field('GOODS_v_Flag')[i] != 0) | (cat.field('GOODS_i_Flag')[i] != 0) | (cat.field('GOODS_z_Flag')[i] != 0):
            ax = f.add_subplot(111)
            plt.yscale('log')
            if cat.field('GOODS_b_Flag')[i]!=0: 
                ax.plot(a[0][0], cat.field('GOODS_b_Jansky')[i]*pow(10,6), color='blue')
                plt.errorbar(a[0][0], cat.field('GOODS_b_Jansky')[i]*pow(10,6), xerr=50, yerr=cat.field('GOODS_b_JanskyErr')[i]*pow(10,6), ls='', capsize=1.2, capthick=0.8, elinewidth=0.5, color='blue')
            if cat.field('GOODS_v_Flag')[i]!=0: 
                ax.plot(a[1][0], cat.field('GOODS_v_Jansky')[i]*pow(10,6), color='darkcyan')
                plt.errorbar(a[1][0], cat.field('GOODS_v_Jansky')[i]*pow(10,6), xerr=50, yerr=cat.field('GOODS_v_JanskyErr')[i]*pow(10,6), ls='', capsize=1.2, capthick=0.8, elinewidth=0.5, color='darkcyan')
            if cat.field('GOODS_i_Flag')[i]!=0: 
                ax.plot(a[2][0], cat.field('GOODS_i_Jansky')[i]*pow(10,6), color='yellowgreen')
                plt.errorbar(a[2][0], cat.field('GOODS_i_Jansky')[i]*pow(10,6), xerr=50, yerr=cat.field('GOODS_i_JanskyErr')[i]*pow(10,6), ls='', capsize=1.2, capthick=0.8, elinewidth=0.5, color='yellowgreen')   
            if cat.field('GOODS_z_Flag')[i]!=0: 
                ax.plot(a[3][0], cat.field('GOODS_z_Jansky')[i]*pow(10,6), color='orange')
                plt.errorbar(a[3][0], cat.field('GOODS_z_Jansky')[i]*pow(10,6), xerr=50, yerr=cat.field('GOODS_z_JanskyErr')[i]*pow(10,6), ls='', capsize=1.2, capthick=0.8, elinewidth=0.5, color='orange') 
            plt.xlabel('$\mathrm{\lambda\ [\AA}]$', size=13)
            plt.ylabel('$\mathrm{F_{\\nu} \ [\mu Jy}]$', size=13)
            path = '/Users/justin/Documents/Master_Thesis/presentation/plots/sed/sed' + str(i+1) + '.eps'
            plt.savefig(path, format='eps', dpi=1000)
            plt.clf()
    return 0

def luminosity(): #would be nice to get d_l(z) relation
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    big_table = fits.open(path)[1].data
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl04.txt'
    cat7ms = pd.read_csv(path, delim_whitespace=True, header=0, skiprows=[1])
    l_d = [47602.3,117256.7,59338.9,49635.1] #d_l in Mpc
    l_x = big_table.field('Soft_Flux')[526]*4*3.14159*((l_d[0]*3.08568*(10**24))**2)*((10**0.2)-(2**0.2))/(((2*(1+big_table.field('phot_z')[526]))**0.2)-((0.5*(1+big_table.field('phot_z')[526]))**0.2))
    return l_x

def line_flux():
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/big_table.fits'
    big_table = fits.open(path)[1].data
    path = '/Users/justin/Documents/Master_Thesis/data/catalogues/chandra_7ms/tbl04.txt'
    cat7ms = pd.read_csv(path, delim_whitespace=True, header=0, skiprows=[1])
    for i in [526,621,639,713]:
        path_sed = '/Users/justin/Documents/Master_Thesis/softwares/eazy/inputs/OUTPUT/photz_' + str(i+1) + '.temp_sed'
        cat_sed = pd.read_csv(path_sed, delim_whitespace=True, header=None, names=['Wavelength', 'Flux'], skiprows=[0,1])
        integral=0.
        i1=0
        i2=0
        minl=[]
        mini=[]
        for j in range(cat_sed['Flux'].size):
            if j > 1:
                if (cat_sed['Flux'][j-1] < cat_sed['Flux'][j]) & (cat_sed['Flux'][j-1] < cat_sed['Flux'][j-2]):
                    minl.append(cat_sed['Wavelength'][j-1])
                    mini.append(j-1)
        lmid=1216.*(1+big_table.field('phot_z')[i])
        for j in range(len(minl)):
            if j > 0:
                if (minl[j-1]< lmid) & (lmid < minl[j]):
                    i1 = mini[j-1]
                    i2 = mini[j]
        for j in range(i1, i2):
            newflux1=cat_sed['Flux'][j+1]*2.998*(10**18)/(cat_sed['Wavelength'][j+1]**2)
            newflux2=cat_sed['Flux'][j]*2.998*(10**18)/(cat_sed['Wavelength'][j]**2)
            integral+= (cat_sed['Wavelength'][j+1]-cat_sed['Wavelength'][j])*(newflux1+newflux2)*0.5
        print(i+1, cat_sed['Wavelength'][i2]-cat_sed['Wavelength'][i1], integral*(10**(-23)))
    return 0