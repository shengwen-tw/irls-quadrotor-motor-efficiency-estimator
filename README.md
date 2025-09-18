# Data-Driven Estimation of Quadrotor Motor Efficiency via Residual Minimization

This repository contains the MATLAB implementation of the algorithms presented in our paper **"Data-Driven Estimation of Quadrotor Motor Efficiency via Residual Minimization"**.

The method formulates trajectory prediction residuals from quadrotor dynamics and solves a constrained optimization problem to recover motor efficiency parameters. To ensure robustness against outliers, the estimator employs **iteratively reweighted least squares (IRLS)**, combined with a **primal–dual interior-point method** that enforces physical constraints and guarantees numerical stability.  

## Citation

Please feel free to cite the paper:

```bibtex
@inproceedings{IRLS_QuadrotorMotorEfficiencyEstimator2026,
  author    = {Sheng-Wen Cheng and Teng-Hu Cheng},
  title     = {Data-Driven Estimation of Quadrotor Motor Efficiency via Residual Minimization},
  year      = {2026},
}
```

## License
This code is released under the BSD-2-Clause license. See [LICENSE](LICENSE) for details.
