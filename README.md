# penlogit
#### A Stata user-written command for penalized logistic regression

- Current version: `1.1.0` 
- Release date: `17jul2017`

---


### Description of the command
`penlogit`  estimates penalized logistic regression models for a binary response via data augmentation. It allows the user to impose Normal and
    generalized log-F prior distributions on one or more model parameters (log odds-ratios).


### Installation
To install the last version of `penlogit` from GitHub, type
```
net install penlogit, from("https://raw.githubusercontent.com/anddis/penlogit/master/")
```
from within a web-aware Stata.


### Authors
Andrea Discacciati and Nicola Orsini (Karolinska Institutet, Stockholm, Sweden)

Sander Greenland (University of California, Los Angeles, CA)


### References
Andrea Discacciati, Nicola Orsini, and Sander Greenland. _Approximate Bayesian logistic regression via penalized likelihood by data augmentation._ The Stata Journal 2015; 15(3):712–736