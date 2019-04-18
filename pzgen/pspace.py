# import csv
# import timeit
# import os
# import pickle as pkl

import numpy as np
import pomegranate
from pomegranate import GeneralMixtureModel as GMM
from pomegranate.distributions import IndependentComponentsDistribution as ICD
from pomegranate.distributions import NormalDistribution as ND
from pomegranate.distributions import UniformDistribution as UD

import pzgen
from pzgen import util as u
# # import chippr
# # from chippr import defaults as d
# # from chippr import utils as u
# # from chippr import sim_utils as su
# # from chippr import gauss
# # from chippr import discrete
# # from chippr import multi_dist
# # from chippr import gmix
# # from chippr import catalog_plots as plots

class PSpace(object):

    def __init__(self, grid, vb=True):
        """
        Object containing catalog of photo-z interim posteriors

        Parameters
        ----------
        grid: numpy.ndarray, float
            grid over which the probability space is defined
        vb: boolean, optional
            True to print progress messages to stdout, False to suppress
        """
        self.grid = grid
        self.npts = len(self.grid)
        self.inds = range(self.npts)

    def make_space(self, bparams={}, sparams={}, oparams={}, vb=True):
        """
        Makes the continuous 2D probability distribution

        Parameters
        ----------
        specs:
        vb: boolean
            print progress to stdout?

        Returns
        -------

        Notes
        -----
        TO DO: only one outlier population at a time for now, will enable more
        TO DO: also doesn't yet include perpendicular features from passing between filter curves, should add that
        """
        hor_funcs = GMM([UD(np.array([self.grid[kk], self.grid[kk+1]])) for kk in self.inds], weights=np.ones_like(self.grid))

        eval_grid = self._make_bias(self.grid, vb=vb, bparams=**bparams)

        sigmas = self._make_scatter(eval_grid, vb=vb, sparams=**sparams)

        vert_funcs = [ND(eval_grid[kk], sigmas[kk]) for kk in self.inds]


        comb_funcs = self._make_outliers(vert_funcs, vb=vb, oparams=**oparams)

        # true n(z) in z_spec, uniform in z_phot
        # grid_amps *= true_func.evaluate(grid_means)
        p_space = [ICD([hor_funcs[kk].dist, comb_funcs[kk].dist]) for kk in self.inds]

        return(p_space)

    def _make_bias(self, x, vb, **bparams):#bvar=False, bconst=-0.003):
        """
        Introduces global redshift bias

        Parameters
        ----------
        bparams: dict
            keyword variables for bias

        Returns
        -------
        y: numpy.ndarray, float
            cental redshifts to use as Gaussian means
        """
        if 'bconst' not in bparams.keys():
            y = x
        else:
            bias = np.asarray(bparams['bconst'])
            if 'bvar' not in bparams.keys() or bparams['bvar'] == False:
                y = x + (np.ones_like(x) * bias[np.newaxis])
            else:
                y = x + ((np.ones_like(x) + x) * bias[np.newaxis])
        return(y)

    def _make_scatter(self, x, vb, **sparams):
        """
        Makes the intrinsic scatter

        Parameters
        ----------
        x: numpy.ndarray, float
            the x-coordinate values upon which to base the intrinsic scatter

        Returns
        -------
        sigmas: numpy.ndarray, float
            the intrinsic scatter values for each galaxy
        """
        sigma = np.ones_like(x)
        if 'sconst' not in sparams.keys():
            const = u.eps
        else:
            const = np.asarray(np.max(u.eps, sparams['sconst']))
        if 'svar' not in sparams.keys() or sparams['svar'] == False:
            sigmas = sigma * const[np.newaxis]
        else:
            sigmas = sigma * (np.ones_like(x) + x)
        return(sigmas)

    def _make_outliers(self, vf, vb, **oparams):

        if 'ofrac' not in oparams.keys() or oparams['ofrac'] == 0.):
            comb_funcs = vf
        else:
            try:
                assert oparams['ofrac'] > 0. and oparams['ofrac'] < 1.
            except AssertionError:
                print('Warning: outlier fraction < 0 or > 1 is being replaced by 0.')
            relfracs = np.array([oparams['ofrac'], 1. - oparams['ofrac']])
            if 'opops' not in oparams.keys()
        if len(oparams) == 0:
            comb_func = vf
        else:
            relfracs = np.array([opop.frac for opop in opops])
            assert np.sum(relfracs) < 1.
            # WILL REFACTOR THIS TO ADD SUPPORT FOR MULTIPLE OUTLIER POPULATIONS
        if self.params['catastrophic_outliers'] != '0':
            frac = self.params['outlier_fraction']
            rel_fracs = np.array([frac, 1. - frac])
            uniform_lf = discrete(np.array([self.bin_ends[0], self.bin_ends[-1]]), np.array([1.]))
            if self.params['catastrophic_outliers'] == 'uniform':
                # use_frac = np.max((0., frac-0.01))
                grid_funcs = [gmix(rel_fracs, [uniform_lf, vert_funcs[kk]], limits=(self.bin_ends[0], self.bin_ends[-1])) for kk in range(self.n_tot)]
            else:
                outlier_lf = gauss(self.params['outlier_mean'], self.params['outlier_sigma']**2)
                # in_amps = np.ones(self.n_tot)
                if self.params['catastrophic_outliers'] == 'template':
                    grid_funcs = [gmix(rel_fracs, [outlier_lf, vert_funcs[kk]], limits=(self.bin_ends[0], self.bin_ends[-1])) for kk in range(self.n_tot)]
                    # out_amps = uniform_lf.pdf(grid_means)
                elif self.params['catastrophic_outliers'] == 'training':
                    full_pdf = np.exp(u.safe_log(outlier_lf.pdf(self.z_fine)))
                    intermediate = np.dot(full_pdf, np.ones(self.n_tot) * self.dz_fine)
                    full_pdf = full_pdf / intermediate[np.newaxis]
                    # flat_pdf = np.exp(u.safe_log(uniform_lf.pdf(self.z_fine)))
                    # flat_pdf = flat_pdf / np.dot(flat_pdf, self.dz_fine)
                    # items = np.array([vert_funcs[kk].pdf(self.z_fine[kk]) for kk in range(self.n_tot)])
                    fracs = np.array([full_pdf, np.ones(self.n_tot)]).T
                    fracs = fracs * np.array([frac, 1.-frac])[np.newaxis, :]
                    grid_funcs = [gmix(fracs[kk], [uniform_lf, vert_funcs[kk]], limits=(self.bin_ends[0], self.bin_ends[-1])) for kk in range(self.n_tot)]
                    # out_funcs = [multi_dist([uniform_lfs[kk], uniform_lf]) for kk in range(self.n_tot)]
                    # out_amps = self.outlier_lf.pdf(grid_means)

                # out_amps /= np.dot(out_amps, self.bin_difs_fine)
                # in_amps *= (1. - self.params['outlier_fraction'])
                # out_amps *= self.params['outlier_fraction']
                # try:
                #     test_out_frac = np.dot(out_amps, self.bin_difs_fine)
                #     assert np.isclose(test_out_frac, self.params['outlier_fraction'])
                # except:
                #     print('outlier fraction not normalized: '+str(test_out_frac))
                # grid_funcs = [gmix(np.array([in_amps[kk], out_amps[kk]]), [grid_funcs[kk], out_funcs[kk]]) for kk in range(self.n_tot)]
                # np.append(grid_means, [self.params['outlier_mean'], self.uniform_lf.sample_one()])
            else:
                grid_funcs = vert_funcs

        return grid_funcs

    def sample_space(self, howmany=1, saveto=None, vb=True):
        pass

    def evaluate_pdfs(self, evalat=np.zeros_like(nd), vb=True):
        pass
