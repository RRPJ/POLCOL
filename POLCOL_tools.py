import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS

import rht
import RHT_tools
from scipy.optimize import curve_fit

from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist.grid_finder import FixedLocator
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist.grid_finder import DictFormatter


def cutout(filename, position, size):
    """ Cut out a square of size size at center position position of the filename filename. Returns: the cut out filename, so it can be used in other scripts. """

    # Load the image and the WCS
    hdu = fits.open(filename)[0]
    naxis = hdu.header['NAXIS']
    wcs = WCS(hdu.header, naxis=naxis)

    # Make the cutout, including the WCS
    cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

    # Put the cutout image in the FITS HDU
    hdu.data = cutout.data

    # Update the FITS header with the cutout WCS
    hdu.header.update(cutout.wcs.to_header())

    # Write the cutout to a new FITS file
    cutout_filename = filename.split(
        '.')[0] + '_cutout.' + filename.split('.')[1]
    hdu.writeto(cutout_filename, overwrite=True)

    print("Made cutout, output: " + cutout_filename)

    # Return the output filename so we can use it in other scripts
    return cutout_filename


def setup_axes(fig, rect, theta, radius):
    """ Sets up the axes in such a way we can plot the polarplot. """

    tr = PolarAxes.PolarTransform()

    pi = np.pi
    angle_ticks = [(-4. / 8. * pi, r"$-\frac{1}{2}\pi$"),
                   (-3. / 8. * pi, r""),
                   (-2. / 8. * pi, r"$-\frac{1}{4}\pi$"),
                   (-1. / 8. * pi, r""),
                   (0. / 8. * pi, r"$0$"),
                   (1. / 8. * pi, r""),
                   (2. / 8. * pi, r"$\frac{1}{4}\pi$"),
                   (3. / 8. * pi, r""),
                   (4. / 8. * pi, r"$\frac{1}{2}\pi$"),
                   (6. / 8. * pi, r"$\frac{3}{4}\pi$"),
                   (8. / 8. * pi, r"$\pi$")]

    grid_locator1 = FixedLocator([v for v, s in angle_ticks])
    tick_formatter1 = DictFormatter(dict(angle_ticks))
    grid_locator2 = MaxNLocator(3)

    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                                      extremes=(
                                                          theta[0], theta[1], radius[0], radius[1]),
                                                      grid_locator1=grid_locator1,
                                                      grid_locator2=grid_locator2,
                                                      tick_formatter1=tick_formatter1,
                                                      tick_formatter2=None)

    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)

    # adjust axis
    # "bottom", "top", "left", "right"
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")

    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")

    # create a parasite axes
    aux_ax = ax1.get_aux_axes(tr)

    return ax1, aux_ax


def make_polarplot(image_filename):
    """ Makes a polarplot of the image_filename. The RHT data is extracted and plotted onto a warped axes from -0.5 to 0.5 pi. This function contains a conversion of the angles hthets to be defined with respect to the RA axis. """

    ipoints, jpoints, hthets, naxis1, naxis2, wlen, smr, thresh = RHT_tools.get_RHT_data(
        image_filename)
    thetas = np.concatenate(hthets)
    thetas_norm = thetas / np.max(thetas)
    thets_arr = RHT_tools.get_thets(wlen, save=False) / np.pi  # Units of pi

    # Convert hthets to be defined w.r.t. RA axis

    w = WCS(image_filename)
    RA_hthets, DEC_hthets = w.all_pix2world(ipoints, jpoints, 0)

    header = fits.getheader(image_filename)
    rotate_by = (header['CRVAL1'] - RA_hthets) / 180  # Units of pi

    big_thets_arr_list = []

    for index, value in enumerate(rotate_by):
        big_thets_arr_list.append((thets_arr + rotate_by[index]) % 1)

    big_thets_arr_list = np.concatenate(big_thets_arr_list)

    # Convert [0, pi] to [-0.5pi, 0,5pi] as in Jelic et al. 2018
    indices_to_shift = np.where(big_thets_arr_list > 0.5)
    big_thets_arr_list[indices_to_shift] -= 1

    # Plotting

    # If you want fill-between, they need to be in order... this does that.
    xdata = big_thets_arr_list * np.pi  # Convert back to good units
    ydata = thetas_norm[np.argsort(xdata)]
    xdata = np.sort(xdata)

    fig = plt.figure(figsize=(18, 9))
    ax1, aux_ax1 = setup_axes(
        fig, 111, theta=[-0.5 * np.pi, 0.5 * np.pi], radius=[0, ydata.max() * 1.05])
    aux_ax1.plot(xdata, ydata, alpha=0.75)
#    aux_ax1.fill_between(xdata, y1=ydata, y2=0, facecolor='blue', interpolate=True, alpha=0.10)
    #plt.axhline(np.max(ydata), color='red', alpha=0.75)
    plt.savefig('polarplot_' + os.path.basename(image_filename) +
                '.png', dpi=600, bbox_inches='tight')
    plt.show()


def get_xpix_and_ypix(catalogue, image_filename):
    """ Returns ALL the pixels of the stars in catalogue on image_filename. """

    data = pd.read_csv(catalogue)
    RA = data['ra']
    DEC = data['dec']

    w = WCS(image_filename)
    xpix, ypix = w.all_world2pix(RA, DEC, 0)

    return xpix, ypix


def get_indices_of_stars(catalogue, image_filename):
    """ Returns the indices of where in the catalogue the stars are in the image. This is useful because we can then later extract polarisation angle, error and distance from these indices. """

    xpix, ypix = get_xpix_and_ypix(catalogue, image_filename)

    # Restrict to image boundaries
    header = fits.getheader(image_filename)
    xbound = header['NAXIS1']
    ybound = header['NAXIS2']

    good_indices = []

    # Find out which pixels are in the image
    for index, value in enumerate(xpix):
        if(value > 0 and value < xbound):
            if(ypix[index] > 0 and ypix[index] < ybound):
                good_indices.append(index)

    return good_indices


def restrict_xpix_and_ypix(catalogue, image_filename):
    """ Returns (xpix, ypix) of all stars in catalogue that are in image with image_filename. """

    good_indices = get_indices_of_stars(catalogue, image_filename)

    xpix, ypix = get_xpix_and_ypix(catalogue, image_filename)

    # Reduce to pixels in image
    xpix = [xpix[i] for i in good_indices]
    ypix = [ypix[i] for i in good_indices]

    return np.transpose([xpix, ypix])


def draw_squares(catalogue, image_filename, square_size):
    """ "Draws" squares around the stars found in get_xpix_and_ypix. This function outputs the pixel values that the RHT angle must be extracted from. Square_size must be odd, as we draw (square_size - 1) / 2 dots on each side. """

    if(square_size % 2 == 0):
        print("Please input an odd square_size")
        return -1

    XX = []
    YY = []

    star_pixcoords = restrict_xpix_and_ypix(catalogue, image_filename)
    rounded = np.round(star_pixcoords)

    s = (square_size - 1) / 2

    for el in rounded:
        pixel_range = np.transpose(np.linspace(
            el - s, el + s, square_size))

        xx, yy = np.meshgrid(pixel_range[0], pixel_range[1])
        # plt.plot(xx, yy, marker='.', linestyle='none')
        XX.append(xx)
        YY.append(yy)

    # plt.show()

    return XX, YY


def extract_data(catalogue, image_filename):
    """ Gets the RHT angle average from the pixels in the square around each star. Also gets the star's polarisation angle and error, aswell as its distance and error. """

    # Star polarisation angle, polarisation error and distance
    data = pd.read_csv(catalogue)
    good_indices = get_indices_of_stars(catalogue, image_filename)

    pol_angle = data['pa'][good_indices]
    err_pol_angle = data['e_pa'][good_indices]
    dist = data['r_est'][good_indices]
    # TO BE IMPLEMENTED: DISTANCE ERROR

    # NOTE: angles from squares not implemented yet.
    # XX, YY = draw_squares(catalogue, image_filename, square_size)

    # RHT angle average calculation
    ipoints, jpoints, hthets, naxis1, naxis2, wlen, smr, thresh = RHT_tools.get_RHT_data(
        image_filename)

    h = np.sum(hthets, axis=0)
    h_norm = h / np.max(h)  # Total intensity per angle found, normalised.

    thets_arr = RHT_tools.get_thets(
        wlen, save=False)  # x-axis containing angles
    thets_deg = thets_arr * 180 / np.pi

    return pol_angle, err_pol_angle, dist, thets_deg, h_norm


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def plot_dist_polangle(catalogue, image_filename):
    """ Plots the distance versus polarisation angle plot. This contains the stars, but also the RHT angle around the stars, averaged. This means that you assume that you are looking at one "zone" with a related RHT angle. """

    pol_angle, err_pol_angle, dist, thets_deg, h_norm = extract_data(
        catalogue, image_filename)

    # Shift lower half of plot since polarisation is modulo pi
    thets_shifted = thets_deg
    thets_shifted[0:59] += 180
    o = np.argsort(thets_shifted)
    thets_sorted = np.sort(thets_shifted)

    # Find normal distribution fit on data
    mu_est = thets_sorted[np.where(h_norm[o] == 1)]
    sig_est = 20  # Not sure about this one, mostly goes well
    popt, pcov = curve_fit(gaussian, thets_sorted,
                           h_norm[o], p0=[int(mu_est), sig_est])

    plt.plot(thets_sorted, h_norm[o], label='RHT angle distribution')
    plt.plot(thets_sorted, gaussian(thets_sorted,
                                    popt[0], popt[1]), label='Gaussian fit')
    plt.xlabel('Polarisation angle (degrees)')
    plt.ylabel('Intensity (normalised)')
    plt.legend()
    plt.savefig('rhtdistribution_' + os.path.basename(image_filename) +
                '.png', dpi=600, bbox_inches='tight')
    plt.show()

    # Put star as close to RHT bar as possible: (modulo pi)
    l = list(pol_angle)
    # Stars "too low" on plot
    i = [x for x in np.concatenate(np.where(popt[0] - l > 90))]
    # Stars "too high" on plot
    j = [y for y in np.concatenate(np.where(l - popt[0] > 90))]

    for index in i:
        l[index] += 180
    for index in j:
        l[index] -= 180

    plt.errorbar(dist, l, xerr=[dist - dist_err_low, dist_err_high - dist],
                 yerr=err_pol_angle, linestyle='none', marker='*', markersize=10, label='stars')
    plt.barh(popt[0], height=popt[1], width=1.1 * np.max(dist), alpha=0.5,
             facecolor='g', label='RHT angle 1-sigma range')
    plt.xlabel('Distance (pc)')
    plt.ylabel('Polarisation angle (degrees)')
    plt.legend()
    plt.savefig('distanceplot_' + os.path.basename(image_filename) +
                '.png', dpi=600, bbox_inches='tight')
    plt.show()
