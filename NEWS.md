# bayesSSM 0.7.0

* Add `auxiliary_filter` and `resample_move_filter` for auxiliary and 
  resample-move particle filters.  
  See all supported filters in `?particle_filter`.  
* Refactor `pmmh` and internals to support multiple particle filters.  
* Add `pf_wrapper` argument to `pmmh` to specify which particle filter to use.  
* Improve argument checking by switching to the `checkmate` package, 
  providing clearer and more consistent error messages.  
* Make minor edits to README and vignettes for clarity.  

# bayesSSM 0.6.1

* Improved reproducibility of `pmmh`: Setting a seed now ensures consistent 
results regardless of the number of cores used.
* Improvement performance of `pmmh` slightly. 
* Rewrote resampling step in `particle_filter` to C++ for improved 
performance.
* Added new article `detailed-overview` which goes over some theory, 
implementation details, and tips for using the package effectively.
* Improved the documentation of the package, including the README and
articles, to provide clearer explanations and examples.

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
