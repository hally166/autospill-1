This fork removes mclapply (parallel processing) allowing AutoSpill to work on Windows machines.  This is my test version and it has only been superficially tested against MAC and the [VIB cloud implementation](https://autospill.vib.be/public/#/home). You should almost certainly be looking at [https://github.com/carlosproca/autospill](https://github.com/carlosproca/autospill) instead. 

Check out the using Autospill R file too: https://github.com/hally166/autospill-1/blob/master/using_AutoSpill.r

# autospill

The **`autospill`** package implements the AutoSpill algorithm for calculating 
spillover coefficients, used to compensate or unmix flow cytometry data. 

For more details, please see:  
Roca *et al*: AutoSpill: A method for calculating spillover coefficients to 
compensate or unmix high-parameter flow cytometry data. 
*bioRxiv* 2020.06.29.177196; 
[doi:10.1101/2020.06.29.177196](https://doi.org/10.1101/2020.06.29.177196) 
\(2020\). 


## Installation

To install **`autospill`** from this GitHub repository, 
use the function `install_github` in the 
[devtools](https://cran.r-project.org/package=devtools) package. 

```R
library( devtools )

install_github( "hally166/autospill-1" )
```


## Help

You can use the standard help in R.

```R
library( autospillCH ) #changed the name for my testing

? get.marker.spillover
? refine.spillover
```


## Examples

Please see the example scripts in the `batch` folder after installing the 
package. 

The scripts `calculate_compensation_paper.r` and 
`calculate_compensation_paper.sh` can be used to reproduce the results of 
AutoSpill for single-color controls appearing in the paper above. 
For this, you will need to download the datasets (FCS files and auxiliary 
`fcs_control.csv` files) from FlowRepository: 
[MM1 dataset](https://flowrepository.org/id/FR-FCM-Z2SS), 
[HS1 & HS2 datasets](https://flowrepository.org/id/FR-FCM-Z2ST), and 
[Be1 dataset](https://flowrepository.org/id/FR-FCM-Z2SV). 

The scripts `calculate_compensation_website.r` and 
`calculate_compensation_website.sh` can be used to reproduce the results 
obtained at the [AutoSpill website](https://autospill.vib.be). 

