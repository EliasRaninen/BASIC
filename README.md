# Bias Adjusted SIgn Covariance matrix (BASIC)

This repository contains MATLAB code for the paper E. Raninen and E. Ollila, "Bias Adjusted Sign Covariance Matrix," in IEEE Signal Processing Letters, vol. 29, pp. 339-343, 2022, doi: [10.1109/LSP.2021.3134940](https://dx.doi.org/10.1109/LSP.2021.3134940).

## Abstract

The spatial sign covariance matrix (SSCM), also known as the normalized sample covariance matrix (NSCM), has been widely used in signal processing as a robust alternative to the sample covariance matrix (SCM). It is well-known that the SSCM does not provide consistent estimates of the eigenvalues of the shape matrix (normalized scatter matrix). To alleviate this problem, we propose BASIC (Bias Adjusted SIgn Covariance), which performs an approximate bias correction to the eigenvalues of the SSCM under the assumption that the samples are generated from zero mean unspecified complex elliptically symmetric distributions (the real-valued case is also addressed). We then use the bias correction in order to develop a robust regularized SSCM based estimator, BASIC Shrinkage estimator (BASICS), which is suitable for high dimensional problems, where the dimension can be larger than the sample size. We assess the proposed estimator with several numerical examples as well as in a linear discriminant analysis (LDA) classification problem with real data sets. The simulations show that the proposed estimator compares well to competing robust covariance matrix estimators but has the advantage of being significantly faster to compute.

## Authors

- Elias Raninen, Doctoral Candidate, Department of Signal Processing and Acoustics, Aalto University.
- Esa Ollila, Professor, Department of Signal Processing and Acoustics, Aalto University.
