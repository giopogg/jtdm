

## Resubmission
This is a resubmission. In this new version I modified the parameter definition in the function to sample from 
the regression coefficient matrix B. 

## R CMD check results

0 errors | 0 warnings | 0 note


### Resubmission
This is a resubmission. In this new version I have included additional checks on the input data to prevent errors.

### R CMD check results

0 errors | 0 warnings | 0 note


### Resubmission
This is a resubmission. In this new version I have updated package citation and I have also fixed some minor bugs.

### R CMD check results

0 errors | 0 warnings | 0 note



### Resubmission
This is a resubmission. After the manual check by the CRAN's team, I:

* Removed redundant terms in the description of the DECRIPTION file
* Currently, no reference is available concerning the method described in the package (the paper is in revision)
* Wrote TRUE and FALSE instead of T and F
* Ensured that users environment is not modified in vignettes and examples, nor any package is installed

Please notice that examples and vignette calling ‘joint_trait_prob’ are nevere run with parallel = TRUE. Therefefore, it seems to me that the remark that we use more than two cores in examples, vignettes ecc is a false positive

### R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
