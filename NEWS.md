# bayesSSM 0.5.0

* The `particles` argument in `init_fn`, which is passed to `particle_filter` 
and `pmmh`, is **deprecated**. Please use `num_particles` instead. 
A warning will be issued if `particles` is used.
* Added support for time dependency in functions. You can now use `t` in 
`transition_fn` and `likelihood_fn` when passing them to `particle_filter`
and `pmmh`. This allows for time-varying transition and likelihood functions.
* Fixed a bug in `particle_filter` in likelihood calculation causing it to be
shifted by a constant.
* Improved robustness of `pmmh` when encountering very low log-likelihood 
values.
* Added scaling for the proposal covariance when using `"logit"` in `pmmh`.
* Improved documentation: updated package description, 
clarified text in the README and vignette, and added more unit tests.


# bayesSSM 0.4.7

* Initial CRAN submission.
