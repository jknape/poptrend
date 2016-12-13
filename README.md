# poptrend
poptrend is an R package for estimating simple population trends from count survey data using Generalized Additive Mixed Models (GAMM). The package uses [mgcv](https://cran.r-project.org/package=mgcv) as the model fitting engine and provides methods for plotting trends, and computing indices and estimates of population change.

## Installation

#### CRAN version

The main version of poptrend is available from [CRAN](https://cran.r-project.org/package=poptrend)   and can be installed from R by running

```R 
install.packages("poptrend")
```

#### Development version

Alternatively, the most recent development version can be installed directly from Github. To do so, first install the devtools package:

```R
install.packages("devtools")
```

Then run:

```R
library(devtools)
install_github("jknape/poptrend")
```
