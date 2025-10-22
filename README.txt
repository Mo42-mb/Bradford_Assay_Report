Use this script to create an HTML report that transforms absorbance data for Bradford Assay standards into equations that predict the protein concentrations of your samples. 

1. Enter absorbance data according to the template_samples_Bradford.xlsx format. This data can be derived from standard cuvette readings or from plate readers. Due to the regression-prediction method, pathlength corrections are not necessary.

2. Update the working directory (line #20 in .Rmd or #8 in .R) and the file names (lines #21-22 in .Rmd or #9-10 in .R). 

3. Click "Knit" in R Markdown or run the script in R to generate a report comparing predictions from linear and quadratic regressions. 

Take a look at example_template_Bradford_analysis.html in the example folder to see what to expect for an output.
