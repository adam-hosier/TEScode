import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import csv
import struct
from math import sqrt
import warnings

from lmfit.models import GaussianModel, PolynomialModel, VoigtModel
from lmfit import Parameters, Model


class MultiPeakGaussian:

    """
    Fits a sum of peak-like functions (gaussians by default) to an array with
    options for plotting or disk storage of the resulting fit.
    """

    def __init__(self, arr, xs = None, num_peaks=10, resolution=20000, num_poly = 4):
        """
        arr: numpy array or list - of y values of spectra
        xs: numpy array or list - xvalues of spectra, defaults to [0, 1, ..., len(arr) - 1]
        num_peaks: integer - number of peaks to find automatically in the spectra
        resolution: integer - number of points to evaluate fit at for plotting
                              note: a large resolution can result in issues with h5 storage
        num_poly: integer - number of polynomial terms for background in spectra
        """



        self.num_peaks = num_peaks
        self.num_poly = num_poly
        self.arr = np.array(arr)
        self.xs = xs if not xs is None else np.linspace(0, self.arr.size-1, self.arr.size)
        self.resolution = resolution
        self.num_gauss = 0 #counter of number of peaks added. Needed for prefixes
        self.lit_vals, self.rel_intensities = [], [] #Used for composite peaks
        self.manual_x_locs = [] #Used for single peaks


    def poly(self, x, coefs):
        return sum(c * x ** i for i, c in enumerate(coefs))
    #BELOW IS OLD DEPRECATED FUNCTIONS FOR FITTING
    #
    # # gaussian function with a varying number of peaks based on the arguments provided
    # def gauss(self, x, *args):
    #     # background function
    #     poly, gauss = args[:self.num_poly], args[self.num_poly:]
    #     result = self.poly(x, poly)
    #
    #     # ensure that each peak has 2 arguments: amplitude and centroid
    #     if len(gauss) % 3 != 0:
    #         print("odd number of arguments")
    #         return
    #
    #     # iterate through argument pairs and add gaussian function for each pair
    #     i=0
    #     while i < len(gauss):
    #         result += gauss[i] * np.exp(-(x-gauss[i+1])**2/(2*gauss[i+2]**2))
    #         i += 3
    #
    #     return result
    #
    # @staticmethod
    # def custom_gauss(x, amplitude, center, sigma):
    #     y_pred = np.zeros(x.size)
    #     y_pred[np.abs(x-center) < 5*sigma] = amplitude * np.exp(-(x[np.abs(x-center) < 5*sigma]-center)**2 / (2*sigma**2))
    #     return y_pred

    def add_composite(self, lit_vals: list, rel_intensities = None):
        """
        Adds a composite / multipeak feature to be built later.

        lit_vals: list of floats - literature values for the wavelengths/energies
                    example: [501.1, 500.3, 502.4]
        rel_intensties: list of lists of two floats or None - describes relative intensities of lines
                        At least one entry must be None to avoid circular references.


                    example: [(1, 2), None, (0, 1/2)]
                    This will be interpreted as follows:
                    the 0th peak will have 2 times the amplitude of the 1st peak
                    The 1st peak is allowed to vary freely due to the None
                    the 2nd peak will have 1/2 times the amplitude of the 0th peak.

                    Equivalently, you could have written:
                    rel_intensities = [(1, 2), None, (1, 1)]
                    rel_intensities = [None, (0, 1/2), (1, 1)]

        """

        self.lit_vals.append(lit_vals)
        if rel_intensities == None:
            rel_intensities = [None] * len(lit_vals)
        self.rel_intensities.append(rel_intensities)

    def add_gaussian(self, x_position: float):
        """
        Finds the index of the given x-position, and adds a peak to build later

        x_postion: float - x_position to add peak later
        """

        idx = np.argmin(np.abs(self.xs - x_position))
        self.manual_x_locs.append(idx)


    def build_composites(self, model, params, function = 'gaussian'):

        """
        Used in fit() to add composite features to the model and params.

        model: lmfit.Model - model to add features to.
        params: lmfit.Parameters - Parameters object to store feature info in
        function: string - 'gaussian' or 'voigt' depending on peak type
        """

        for lit_vals, rel_intensities in zip(self.lit_vals, self.rel_intensities):
            params.add('g{self.num_gauss)_stretch', value = 1)
            for i, nrg in enumerate(lit_vals):
                params.add(f'lit{i+self.num_gauss}', value = nrg, vary = False)
                if function == 'gaussian':
                    G = GaussianModel(prefix = f'g{i+self.num_gauss}_')
                elif function == 'voigt':
                    G = VoigtModel(prefix = f'g{i+self.num_gauss}_')
                params += G.make_params(amplitude = 300, center = nrg, sigma = 3)
                params[f'g{i+self.num_gauss}_amplitude'].min = 0
                model += G

            for i, intensity in enumerate(rel_intensities):
                if i != 0:
                    params[f'g{i+self.num_gauss}_sigma'].expr = f'g{self.num_gauss}_sigma'
                    params[f'g{i+self.num_gauss}_center'].expr = \
                        f'g{self.num_gauss}_center + stretch * (lit{i+self.num_gauss} - lit{self.num_gauss})'

                if intensity != None:
                    idx, rel = intensity
                    params[f'g{i+self.num_gauss}_amplitude'].expr = f'{rel} * g{idx+self.num_gauss}_amplitude'
            self.num_gauss += len(lit_vals)

        return model, params

    def build_gaussians(self, model, params, same_sigma = False, function = 'gaussian'):

        """
        Used in fit() to add single peaks to the model and params.

        model: lmfit.Model - model to add features to.
        params: lmfit.Parameters - Parameters object to store feature info in
        same_sigma: boolean - whether or not to force all peaks to share sigma / width
        function: string - 'gaussian' or 'voigt' depending on peak type
        """

        y, x = self.arr, self.xs
        for idx in self.manual_x_locs:

            if function == 'gaussian':
                G = GaussianModel(prefix = f'g{self.num_gauss}_')
            elif function == 'voigt':
                G = VoigtModel(prefix = f'g{self.num_gauss}_')

            params += G.make_params(amplitude = y[idx] - np.median(y), center = x[idx], sigma = 2)
            params[f'g{self.num_gauss}_center'].set(min = 0)
            params[f'g{self.num_gauss}_amplitude'].set(min = 0)
            if self.num_gauss != 0 and same_sigma:
                params[f'g{self.num_gauss}_sigma'].expr = 'g0_sigma'
            model += G
            self.num_gauss += 1

        return model, params


    # find the indices of specified number of most prominent peaks in the initial array
    def get_peak_indices(self, min_up = 3, min_down = 2, smooth = True):
        """
        Helper function to automatically detect peaks in the data.

        min_up: int - min number of consecutive increasing steps to detect peak
        min_up: int - min number of consecutive decreasing steps to detect peak
        smooth: boolean - whether to find peaks based off a rolling average of width 5

        example: if min_up = 3 and min_down = 2, arr would have to increase for
        at least 3 steps in a row and then decrease for at least two values in a row
        for a peak to be found.
        """
        if smooth:
            data = [np.mean(self.arr[i:i+5]) for i in range(self.arr.size - 4)]
            data.insert(0, 0)
            data.insert(0, 0)
            data = np.array(data)
        else:
            data = self.arr

        #goes through array and finds sections that satisfy peak reqs
        result = []
        consec_up, consec_down, prev_val = 0, 0, np.inf
        for i, curr_val in enumerate(data):
            if curr_val > prev_val:
                if consec_down > 0: #was going down, now going up
                    if consec_down >= min_down and consec_up >= min_up: #peak reqs
                        result.append(i - (consec_down + 1))
                    consec_up = 0
                consec_up += 1
                consec_down = 0
            else:
                consec_down += 1
            prev_val = curr_val


        result = sorted(result, key = lambda i: self.arr[i], reverse = True)
        result = result[:self.num_peaks]
        return np.array(result)


    def fit(self, return_dict = None,  function = 'gaussian', same_sigma = False, auto_find_peaks = True):
        """
        Primary function of class, finds parameters of peaks to best fit data.

        return_dict: dict - used for threading, if given stores results in dict
        function: string - 'gaussian' or 'voigt' depending on peak type
        same_sigma: boolean - whether to force peaks to share sigma / width
        auto_find_peaks: boolean - if True, self.num_peaks are automatically identified.
        """
        x, y = self.xs, self.arr
        y[y < np.median(y) * 0.8] = np.median(y) * 0.8

        #initialize polynomial background with linear coefficients
        L = PolynomialModel(degree = 1)
        params = L.guess(y, x=x)
        for i in range(2, self.num_poly):
            params.add(f'c{i}', value = 0)

        M = PolynomialModel(degree = self.num_poly - 1)


        if auto_find_peaks:
            for idx in self.get_peak_indices():
                self.manual_x_locs.append(idx)

        M, params = self.build_gaussians(M, params, same_sigma, function = function)
        M, params = self.build_composites(M, params, function = function)


        print('Going into fit.')
        out = M.fit(y, params, x = x, max_nfev = 10**7, verbose = True)
        print('Fit Completed.')

        params.update(out.params)
        #stores vital parameters in list for later use
        popt = []
        for i in range(self.num_poly):
            popt.append(out.params[f'c{i}'].value)
        for i in range(len(self.get_peak_indices())):
            popt.append(out.params[f'g{i}_height'].value)
            popt.append(out.params[f'g{i}_center'].value)
            popt.append(out.params[f'g{i}_sigma'].value)
            popt.append(out.params[f'g{i}_fwhm'].value)
        self.popt = popt


        #stores data for each peak in its own row.
        new = []
        for i in range(len(self.get_peak_indices())):
            temp = []
            temp.append([out.params[f'g{i}_center'].value, out.params[f'g{i}_center'].stderr])
            temp.append([out.params[f'g{i}_height'].value, out.params[f'g{i}_height'].stderr])
            temp.append([out.params[f'g{i}_sigma'].value, out.params[f'g{i}_sigma'].stderr])
            temp.append([out.params[f'g{i}_fwhm'].value, out.params[f'g{i}_fwhm'].stderr])
            new.append(temp)


        #sorts gaussian data by amplitude and stores centroids in self.centroids
        new = sorted(new, key=lambda x: x[1][0], reverse=True)
        self.centroids = []
        for i in new:
            self.centroids.append(i[0])
        self.gauss_data = new

        self.x0 = np.linspace(min(x), max(x), self.resolution)
        self.y0 = out.eval(x = self.x0)
        self.plist = new

        if not return_dict is None:
            return_dict['popt'], return_dict['centroids'] = self.popt, self.centroids
            return_dict['x0'], return_dict['y0'], return_dict['xs'] = self.x0, self.y0, self.xs
            return_dict['gauss_data'] = self.gauss_data
            return_dict['rez'] = self.plist
            return_dict['output'] = out
            return_dict['params'] = out.params

        return self.x0, self.y0

    # get centroids in ascending order of x-value
    def get_centroids_ascending(self):
        return np.array(sorted(self.centroids, key=lambda x: x[0], reverse=False))

    # get centroids in descending order of x-value
    def get_centroids_descending(self):
        return np.array(sorted(self.centroids, key=lambda x: x[0], reverse=True))


    def plot_fit(self, normalize_background = True, xtransform = None,
                    xlabel = 'Channel Number (-)', title = 'Comparison', ax = None, linestyle = '-r'):

        """
        Plotting function to verify fit quality.

        normalize_background: boolean - whether to subtract polynomial background
        xtransform: (function, args) or None - function to apply to x-axis
                    example: (self.poly, (1,2,3)) would map all x-values to
                             x_new = self.poly(x, 1, 2, 3)
        xlabel: string - text on x-axis
        title: string - title text
        ax: plt.axis - if given, plotting will be done on the given axis. Use
                       if you want to plot other data on the same plot.
        linestyle: string - matplotlib.pyplot line specifier, eg: '--b', '-.k'
        """

        xFit = self.x0
        xArr = self.xs

        yFit = self.y0
        yArr = self.arr

        if normalize_background:
            yFit = yFit  - self.poly(xFit, self.popt[:self.num_poly])
            yArr = yArr
            #yArr = yArr - self.poly(xArr, self.popt[:self.num_poly])

        if xtransform:
            func, popt = xtransform
            xFit = func(xFit, *popt)
            xArr = func(xArr, *popt)


        if ax is not None:
            plt.sca(ax)
            show = False
        else:
            f, ax  = plt.subplots()
            show = True

        plt.plot(xFit, yFit, linestyle , label = 'MultiGaussian Fit', c='b')
        plt.plot(xArr, yArr, '--' + linestyle[-1], label = 'Original Data')
        plt.xlabel(xlabel)
        plt.ylabel('Intensity (ADU)')
        plt.title(title)
        plt.legend()

        if show:
            plt.show()

    def to_hdf(self, dataset):
        """
        stores fit data in header of the hdf5 dataset for later use
        """
        dataset.attrs['popt'], dataset.attrs['centroids'] = self.popt, self.centroids
        dataset.attrs['x0'], dataset.attrs['y0'], dataset.attrs['xs'] = self.x0, self.y0, self.xs
        dataset.attrs['gauss_data'] = self.gauss_data

    def read_hdf(self, dataset):
        """
        reads fit data from header of the hdf5 dataset to save runtime
        """
        try:
            self.popt, self.centroids = dataset.attrs['popt'], dataset.attrs['centroids']
            self.x0, self.y0, self.xs = dataset.attrs['x0'], dataset.attrs['y0'], dataset.attrs['xs']
            self.gauss_data = dataset.attrs['gauss_data']
            return True
        except:
            print('Could not load MultiPeakGaussian data from file.')
            return False



class SpeReader:
    """
    SPE files have a header of 4100 bits which store information about
    the image.
    Relevant positions(parameter/position/size):
    xdim/42/16bits
    ydim/656/16bits
    frames/1446/16bits
    rawdate/20/72bits
    rawtime/172/6/48bits

    The recorded intensities start from bit 4100.
    """

    def __init__(self, fname):
        self.fid = open(fname, 'rb')
        self.load_size()
        self.load_metadata()

    def load_size(self):
        """
        Load the dimensions of the data.
        """
        self.xdim = np.int64(self.read_at(42, 1, np.int16)[0])
        self.ydim = np.int64(self.read_at(656, 1, np.int16)[0])
        self.frame = np.int64(self.read_at(1446, 1, np.int16)[0])

    def load_metadata(self):
        """
        Loads metadata from the spe header and stores it in self.metadata dict
        """

        self.metadata = dict()
        self.metadata['gain'] = np.int64(self.read_at(198,1,np.int16)[0])
        self.metadata['adcrate'] = np.int64(self.read_at(190,1,np.int16)[0])
        self.metadata['adcresolution'] = np.int64(self.read_at(192,1,np.int16)[0])

        self.metadata['temp'] = np.float64(self.read_at(118,1,np.float32)[0])

        #self.metadata['offset'] = np.float64(self.read_at(78,1,np.float32)[0])
        #self.metadata['finalwl'] = np.float64(self.read_at(82,1,np.float32)[0])

        self.metadata['rawdate'] = self.read_text_at(20,72)
        self.metadata['rawtime'] = self.read_text_at(172,48)
        comments = (self.read_text_at(200+80*i,80).replace('\x00','') for i in range(5))
        self.metadata['comments'] = tuple(comments)



    def get_size(self):
        """
        Return the dimensions of the data.
        """
        return (self.frame, self.xdim, self.ydim)


    def read_text_at(self,pos,bits):
        self.fid.seek(pos)
        bytes = self.fid.read(bits // 8)
        return bytes.decode('ascii')



    def read_at(self, pos, size, ntype):
        """
        Binary reader.
        pos: starting position
        size: number of items to read
        ntype: data type
        """
        self.fid.seek(pos)
        return np.fromfile(self.fid, ntype, size)

    def load_img(self):
        """
        Load img into a numpy array.
        """

        data = self.read_at(4100, self.xdim * self.ydim * self.frame, np.uint16)
        return data.reshape((self.frame, self.xdim))

    def close(self):
        """
        Close file.
        """
        self.fid.close()

    def print_metadata(self):
        for key, value in self.metadata.items():
            print(f'{key}: {value}')


class CosmicRayFilter:
    """
    Cosmic Ray Filter class for removal of cosmic rays from spectra.
    """
    def __init__(self, filterval = 5):
        self.filterval = filterval

    def apply(self, data, combine = True):
        """
        Apply filter to input data.

        Returns the spectra with cosmic rays removed.
        """
        data = data.astype(float)
        nr_frames, camera_size = data.shape

        if nr_frames == 1:
            warnings.warn('Only 1 frame, can''t detect cosmic rays', RuntimeWarning)
            return data[0]

        for pixel in range(camera_size): # loop through each pixel
            allframes = data[:, pixel]
            cosmic_frames = []
            for frame in range(nr_frames): # loop through each frame

                testval = allframes[frame]
                testframes = np.delete(allframes, frame)
                poisson_noise = sum([sqrt(val) for val in testframes])/(nr_frames-1)
                mean = sum(testframes)/(nr_frames-1)

                # if outside range add it's index to a list
                if testval > (mean + (self.filterval*poisson_noise)):
                    cosmic_frames.append(frame)

            if len(cosmic_frames) > 0: # if cosmic ray found replace with mean value of other pixels
                non_cosmic_avg = np.mean(np.delete(allframes, cosmic_frames))
                data[cosmic_frames, pixel] = non_cosmic_avg

        if combine:
            return np.sum(data,axis = 0)
        else:
            return data


if __name__ == '__main__':
    S, C = SpeReader('./270.spe'), CosmicRayFilter(5)
    print(S.load_img().shape)
    print(S.print_metadata())
    img = C.apply(S.load_img())
    multi = MultiPeakGaussian(img, num_peaks = 35)

    plt.plot(img)
    plt.show()
    plt.close()

    multi.fit(same_sigma = True, function = 'gaussian')
    multi.plot_fit()
