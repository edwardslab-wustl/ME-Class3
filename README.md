# ME-Class3

## Requirements
 * Python Packages
     *  numpy
     *  pysam
     *  biopython

## Installation
Clone and install with pip:

````
    git clone https://github.com/edwardslab-wustl/ME-Class3.git
    cd ME-Class3
    pip install .
````

## Usage

### Input files

### Example

## License information


## Parameters

```

1. Process samples: 

	This module processes fractional methylaton files and generates differentail methylation files. 

    -ctp	: Cell type Fofn (Please refer to example data). 
    -pto	: Path to Output.
    -expr	: Name of expression file. Samples without expression information will be omitted.

2. Generate proximity list:

	This module generates list of enhancers based on their proximity to respective genes. 

	-pst	:	(+)end start index.
	-ped	:	(+)end end index.
	-nst	:	(-)end start index.
	-ned	:	(-)end end index.
	-glf	:	Gene list file.

3. Merge Gene and Enhancer features togather. 
		
	This module addes gene and enhancer features togather.

	-pst    :   (+)end start index.
    -ped    :   (+)end end index.
    -nst    :   (-)end start index.
    -ned    :   (-)end end index.
    -tif    :   TSS interpolation file.

4. Interpolation.
	
	This module generates interpolation of HRPS and regulatory regions. 

	-sigma	:	Value of sigma. (default=50)
	-nip	:	Number of interp points. (default=500)
	-ibin	:	Size of bin around TSS/GB. (default=5000)
	-ach	:	Anchor window. (default=100000). ( Used for flank normalization)
	-rfl	:	RE flnk length. (default=500)
	-rff	:	RE flnk features. (default=25)
	-rfn	:	Region fofn file.
	-tag	:	Tag for output.
	-sfn	:	Sample file name.
	-fln	:	Flank Normalized interpolation. (default=False)
	-gsl	:	Interpolation  for slected genes?. (default=False)
	-frf	:	Fixed features for RE flank. (default=False)
	-mmg	:	Minimum CpGs assayed. (default=40)
	-mmr	:	Minimum CpGs assayed. (default=2 )

4. Add expression information.

	-intrp	:	Interpolation CSV file.
	-expr	:	Expression file.
	-fef	:	Floor expression value? (default=True)
	-efv	:	Expression floor value. (default=5.0)
	-def	:	Method for differential expression. (default=1)

5. Run classification.

	-dfi	:	Interpolation-expression dataframe.
	-ntr	:	Number of trees for Random Forest Classifier. (default=5001)
	-npr	:	Number of Processors for RF run. (default=8)
	-tag	:	Tag for Output Writing. (test)
	-fsl	:	Feature Selection. 1: TSS; 2: TSS+RE, 3: RE.
	-suf	:	Shuffle true ot false. (default=True)

