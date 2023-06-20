#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
from scipy.special import gammaln, digamma, polygamma

# import statsmodels.api as sm


class Noise_Model:
    def __init__(self, response, covariates, constant, distribution):
        self.response = response
        self.X = covariates  # X
        self.C = constant  # size ncells
        self.distribution = distribution
        self.ncells, self.ncovs = np.shape(self.X)

    def loglikelihood(self, params):

        if self.distribution == "Poisson":
            lam = params
            logl = self.response * np.log(lam) - lam - gammaln(self.response + 1)
            return logl

        elif self.distribution == "NB" or self.distribution == "GeneralDispersion":
            mu = params[:, 0]
            alpha = params[:, 1]
            logl = self.response * np.log(alpha * mu)
            logl -= (self.response + 1 / alpha) * np.log(1 + alpha * mu)
            logl += gammaln(self.response + 1 / alpha)
            logl -= gammaln(1 / alpha)
            logl -= gammaln(self.response + 1)
            return logl

        elif self.distribution == "ZeroPoisson":
            p = params[:, 0]
            lam = params[:, 1]
            logl = np.zeros(self.ncells)
            indices_zero = np.where(self.response < 0.1)
            indices_pos = np.where(self.response > 0.1)
            p_zero = p[indices_zero]
            logl[indices_zero] = np.log(
                p_zero + (1 - p_zero) * np.exp(-lam[indices_zero])
            )
            logl[indices_pos] = (
                np.log(1 - p[indices_pos])
                + self.response[indices_pos] * np.log(lam[indices_pos])
                - lam[indices_pos]
                - gammaln(self.response[indices_pos] + 1)
            )
            return logl

        else:
            # If none of the previous is chosen, use Normal with log link
            mu = params[:, 0]
            sigma = params[:, 1]
            logl = -np.log(2 * np.pi * sigma)
            logl -= (self.response - mu) ** 2 / sigma
            logl /= 2
            return logl

    def score_statistic(self, params):

        if self.distribution == "Poisson":
            lam = params
            return self.response - lam

        elif self.distribution == "NB":
            mu = params[:, 0]
            alpha = params[:, 1]
            u_mu = (self.response - mu) * mu / (mu + alpha * mu**2)
            u_alpha = -np.sum(
                digamma(self.response + 1 / alpha)
                - digamma(1 / alpha)
                - np.log(1 + alpha * mu)
                - -alpha * (self.response - mu) / (1 + alpha * mu)
            ) / (alpha[0] ** 2)
            return u_mu, u_alpha

        elif self.distribution == "ZeroPoisson":
            p = params[:, 0]
            lam = params[:, 1]
            eps = 0.01
            indices_zero = np.where(self.response < eps)
            indices_pos = np.where(self.response > eps)
            u_p = np.zeros(self.ncells)
            u_p[indices_zero] = (
                (1 - p[indices_zero])
                * p[indices_zero]
                * (1 - np.exp(-lam[indices_zero]))
                / (p[indices_zero] + (1 - p[indices_zero]) * np.exp(-lam[indices_zero]))
            )
            u_p[indices_pos] = -p[indices_pos]
            u_l = np.zeros(self.ncells)
            u_l[indices_zero] = (
                -np.exp(-lam[indices_zero])
                * lam[indices_zero]
                * (1 - p[indices_zero])
                / (p[indices_zero] + (1 - p[indices_zero]) * np.exp(-lam[indices_zero]))
            )
            u_l[indices_pos] = self.response[indices_pos] - lam[indices_pos]
            return u_p, u_l

        elif self.distribution == "GeneralDispersion":
            mu = params[:, 0]
            alpha = params[:, 1]
            u_mu = (self.response - mu) * mu / (mu + alpha * mu**2)
            return u_mu

        else:
            mu = params[:, 0]
            sigma = params[:, 1]
            return mu * (self.response - mu) / sigma

    def hessian_weights(self, params):

        if self.distribution == "Poisson":
            return params

        elif self.distribution == "NB":
            mu = params[:, 0]
            alpha = params[:, 1]
            w_mu = mu**2 / (mu + alpha * mu**2)
            u_alpha = -(
                digamma(self.response + 1 / alpha)
                - digamma(1 / alpha)
                - np.log(1 + alpha * mu)
                - -alpha * (self.response - mu) / (1 + alpha * mu)
            ) / (alpha**2)
            w_alpha = np.sum(
                -(2 / alpha) * u_alpha
                + 1
                / alpha**4
                * (
                    polygamma(1, self.response + 1 / alpha)
                    - polygamma(1, 1 / alpha)
                    - alpha
                    - 1 / (1 / alpha + mu)
                    + (self.response - mu) / (mu + 1 / alpha) ** 2
                )
            )
            return w_mu, w_alpha

        elif self.distribution == "ZeroPoisson":
            p = params[:, 0]
            lam = params[:, 1]
            eps = 0.01
            indices_zero = np.where(self.response < eps)
            indices_pos = np.where(self.response > eps)
            w_p = np.zeros(self.ncells)
            p_zero = p[indices_zero]
            p_pos = p[indices_pos]
            num_p = 1 - np.exp(-lam[indices_zero])
            den_p = p_zero + (1 - p_zero) * np.exp(-lam[indices_zero])
            w_p[indices_zero] = ((num_p / den_p) * p_zero * (1 - p_zero)) ** 2 - (
                1 - p_zero
            ) * (1 - 2 * p_zero) * p_zero * num_p / den_p
            w_p[indices_pos] = p[indices_pos] * (1 - p[indices_pos])
            w_l = np.zeros(self.ncells)
            num_l = (1 - p_zero) * np.exp(-lam[indices_zero])
            den_l = p_zero + (1 - p_zero) * np.exp(-lam[indices_zero])
            w_l[indices_zero] = (
                (num_l / den_l)
                * lam[indices_zero]
                * (1 + (num_l / den_l) * lam[indices_zero])
            )
            w_l[indices_pos] = lam[indices_pos]
            return w_p, w_l

        elif self.distribution == "GeneralDispersion":
            mu = params[:, 0]
            alpha = params[:, 1]
            w_mu = mu**2 / (mu + alpha * mu**2)
            return w_mu

        else:
            mu = params[:, 0]
            sigma = params[:, 1]
            return mu**2 / sigma

    def iterate(self, coef):

        if self.distribution == "Poisson":
            tmp = np.matmul(self.X, coef)
            params = np.exp(tmp + self.C)

        elif self.distribution == "NB":
            mu = np.exp(np.matmul(self.X, coef[:-1]) + self.C)
            alpha = coef[-1]
            params = np.concatenate(
                (mu.reshape(-1, 1), alpha * np.ones((self.ncells, 1))), axis=1
            )
            u_mu, u_alpha = self.score_statistic(params)
            w_mu, w_alpha = self.hessian_weights(params)
            z_mu = np.matmul(self.X, coef[:-1]) + u_mu / w_mu
            B_mu = np.linalg.inv(np.matmul(self.X.T * w_mu, self.X))
            v_mu = np.matmul(self.X.T * w_mu, z_mu)
            coef_mu = np.matmul(B_mu, v_mu)
            coef_alpha = alpha - u_alpha / w_alpha
            return np.concatenate((coef_mu, [coef_alpha]))

        elif self.distribution == "ZeroPoisson":
            p = np.exp(np.matmul(self.X, coef[:, 0])) / (
                1 + np.exp(np.matmul(self.X, coef[:, 0]))
            )
            p = p.reshape((self.ncells, 1))
            lam = np.exp(np.matmul(self.X, coef[:, 1])).reshape((self.ncells, 1))
            params = np.concatenate((p, lam), axis=1)
            u_p, u_l = self.score_statistic(params)
            w_p, w_l = self.hessian_weights(params)

            z_p = np.matmul(self.X, coef[:, 0]) + u_p / w_p
            B_p = np.linalg.inv(np.matmul(self.X.T * w_p, self.X))
            v_p = np.matmul(self.X.T * w_p, z_p)
            coef_p = np.matmul(B_p, v_p).reshape((-1, 1))

            z_l = np.matmul(self.X, coef[:, 1]) + u_l / w_l
            B_l = np.linalg.inv(np.matmul(self.X.T * w_l, self.X))
            v_l = np.matmul(self.X.T * w_l, z_l)
            coef_l = np.matmul(B_l, v_l).reshape((-1, 1))

            return np.concatenate((coef_p, coef_l), axis=1)

        elif self.distribution == "GeneralDispersion":
            mu = np.exp(np.matmul(X, coef))
            h = 0.1
            l = 2 * np.pi / h
            counts = np.zeros(int(l) + 1)
            counts_est = np.zeros(int(l) + 1)
            var_est = np.zeros(int(l) + 1)

            for i in range(int(l) + 1):
                sub_col = self.response[
                    np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                ]
                sub_est = mu[
                    np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                ]
                counts[i] = np.mean(sub_col)
                counts_est[i] = np.mean(sub_est)
                var_est[i] = np.mean((sub_col - sub_est) ** 2)

            a, b, c = np.polyfit(counts_est, var_est, 2)
            alpha = np.zeros(self.ncells)
            for i in range(int(l) + 1):
                indices = np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                alpha[indices] = (b - 1) / counts_est[i] + a + c / counts_est[i] ** 2

            params = np.concatenate((mu.reshape(-1, 1), alpha.reshape((-1, 1))), axis=1)
            u_mu = self.score_statistic(params)
            w_mu = self.hessian_weights(params)
            z_mu = np.matmul(self.X, coef) + u_mu / w_mu
            B_mu = np.linalg.inv(np.matmul(self.X.T * w_mu, self.X))
            v_mu = np.matmul(self.X.T * w_mu, z_mu)
            coef_mu = np.matmul(B_mu, v_mu).reshape(-1)
            return coef_mu

        else:
            mu = np.exp(np.matmul(X, coef))
            sigma = np.sum((self.response - mu) ** 2) / self.ncells
            params = np.concatenate(
                (mu.reshape((-1, 1)), sigma * np.ones((self.ncells, 1))), axis=1
            )

        z = np.matmul(self.X, coef) + self.score_statistic(
            params
        ) / self.hessian_weights(params)
        try:
            B = np.linalg.inv(
                np.matmul(self.X.T * self.hessian_weights(params), self.X)
            )
            v = np.matmul(self.X.T * self.hessian_weights(params), z)
            new_coef = np.matmul(B, v)
            return new_coef.reshape(-1)
        except:
            print("except")
            return coef

    def fit(self, n_iter, alpha_init=1.0):

        if self.distribution == "ZeroPoisson":
            coef = np.zeros((self.ncovs, 2))
            coef[0, 0] = 1
            coef[0, 1] = 1

        elif self.distribution == "NB":
            coef = np.zeros(self.ncovs + 1)
            coef[0] = np.log(np.mean(self.response))  # log of mean of rawcounts
            coef[-1] = alpha_init  # initilize alpha

        else:
            coef = np.zeros(self.ncovs)
            coef[0] = np.log(np.mean(self.response))

        for i in range(n_iter):
            coef = self.iterate(coef)

        return coef

    def regularized_iteration(self, coef, tau):

        if self.distribution == "Poisson":
            params = np.exp(np.matmul(self.X, coef) + self.C)

        elif self.distribution == "NB":
            mu = np.exp(np.matmul(self.X, coef[:-1]) + self.C)
            alpha = coef[-1]
            params = np.concatenate(
                (mu.reshape(-1, 1), alpha * np.ones((self.ncells, 1))), axis=1
            )
            u_mu, u_alpha = self.score_statistic(params)
            w_mu, w_alpha = self.hessian_weights(params)
            z_mu = np.matmul(self.X, coef[:-1]) + u_mu / w_mu
            B_mu = np.matmul(self.X.T * w_mu, self.X) + tau * np.identity(self.ncovs)
            v = np.matmul(self.X.T * w_mu, z_mu)
            coef_mu = np.matmul(np.linalg.inv(B_mu), v)
            coef_alpha = alpha - u_alpha / w_alpha
            return np.concatenate((coef_mu, [coef_alpha]))

        elif self.distribution == "ZeroPoisson":
            p = np.exp(np.matmul(self.X, coef[:, 0])) / (
                1 + np.exp(np.matmul(self.X, coef[:, 0]))
            )
            p = p.reshape((self.ncells, 1))
            lam = np.exp(np.matmul(self.X, coef[:, 1])).reshape((self.ncells, 1))
            params = np.concatenate((p, lam), axis=1)
            u_p, u_l = self.score_statistic(params)
            w_p, w_l = self.hessian_weights(params)

            z_p = np.matmul(self.X, coef[:, 0]) + u_p / w_p
            B_p = np.matmul(self.X.T * w_p, self.X) + tau * np.identity(self.ncovs)
            v_p = np.matmul(self.X.T * w_p, z_p)
            coef_p = np.matmul(np.linalg.inv(B_p), v_p).reshape((-1, 1))

            z_l = np.matmul(self.X, coef[:, 1]) + u_l / w_l
            B_l = np.matmul(self.X.T * w_l, self.X) + tau * np.identity(self.ncovs)
            v_l = np.matmul(self.X.T * w_l, z_l)
            coef_l = np.matmul(np.linalg.inv(B_l), v_l).reshape((-1, 1))

            return np.concatenate((coef_p, coef_l), axis=1)

        elif self.distribution == "GeneralDispersion":
            mu = np.exp(np.matmul(self.X, coef))
            h = 0.1
            l = 2 * np.pi / h
            counts = np.zeros(int(l) + 1)
            counts_est = np.zeros(int(l) + 1)
            var_est = np.zeros(int(l) + 1)

            for i in range(int(l) + 1):
                sub_col = self.distribution[
                    np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                ]
                sub_est = mu[
                    np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                ]
                counts[i] = np.mean(sub_col)
                counts_est[i] = np.mean(sub_est)
                var_est[i] = np.mean((sub_col - sub_est) ** 2)

            a, b, c = np.polyfit(counts_est, var_est, 2)
            alpha = np.zeros(self.ncells)
            for i in range(int(l) + 1):
                indices = np.logical_and(phase > i * h, phase < (i + 1) * h).squeeze()
                alpha[indices] = (b - 1) / counts_est[i] + a + c / counts_est[i] ** 2

            params = np.concatenate((mu.reshape(-1, 1), alpha.reshape((-1, 1))), axis=1)
            u_mu = self.score_statistic(params)
            w_mu = self.hessian_weights(params)
            z_mu = np.matmul(self.X, coef) + u_mu / w_mu
            B_mu = np.matmul(self.X.T * w_mu, self.X) + tau * np.identity(self.ncovs)
            v = np.matmul(self.X.T * w_mu, z_mu)
            coef_mu = np.matmul(np.linalg.inv(B_mu), v)
            return coef_mu

        else:
            mu = np.exp(np.matmul(self.X, coef))
            sigma = np.sum((self.response - mu) ** 2) / self.ncells
            params = np.concatenate(
                (mu.reshape((-1, 1)), sigma * np.ones((self.ncells, 1))), axis=1
            )

        z = np.matmul(self.X, coef) + self.score_statistic(
            params
        ) / self.hessian_weights(params)
        B = np.matmul(
            self.X.T * self.hessian_weights(params), self.X
        ) + tau * np.identity(self.ncovs)
        v = np.matmul(self.X.T * self.hessian_weights(params), z)
        new_coef = np.matmul(np.linalg.inv(B), v)

        return new_coef.reshape(-1)

    def regularized_fit(self, n_iter, tau):

        if self.distribution == "ZeroPoisson":
            coef = np.zeros((self.ncovs, 2))
            coef[0, 0] = 1
            coef[0, 1] = 1

        elif self.distribution == "NB":
            coef = np.zeros(self.ncovs + 1)
            coef[0] = np.log(np.mean(self.response))
            coef[-1] = 1

        else:
            coef = np.zeros(self.ncovs)
            coef[0] = np.log(np.mean(self.response))

        for i in range(n_iter):
            coef = self.regularized_iteration(coef, self.ncells * tau)

        return coef


def fourier_basis(phi, num_harmonics):
    basis = np.ones((len(phi), 1 + 2 * num_harmonics))
    for n in range(num_harmonics):
        basis[:, 2 * n + 1] = np.cos((n + 1) * phi)
        basis[:, 2 * n + 2] = np.sin((n + 1) * phi)
    return basis
