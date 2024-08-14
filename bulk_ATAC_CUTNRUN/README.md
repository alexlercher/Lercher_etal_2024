# LercherEtAl
Peak Files for Lercher et al "Antiviral innate immune memory in alveolar macrophages following SARS-CoV-2 infection"


Filtered peaks were generated using the following approach in R:
```
library(tidyverse)
library(MASS)

myRiP <- read.delim('/path/to/file.txt')
myRowMaxs <- apply(
  myRiP,
  1,
  max
)

myRowMaxs %>%
  log10 %>%
  ecdf %>%
  knots %>%
  fitdistr('normal') -> myNormFit

myOutlierThresh <- myNormFit$estimate[1] - 1.5 * myNormFit$estimate[2]

myFilteredData <- myRiP %>% filter((myRowMaxs %>% log10) >= myOutlierThresh)

```
