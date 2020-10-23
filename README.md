CODE OF THE ANALYSIS PERFORMED IN THE T-ALL RELAPSE EVOLUTION IN ADULT PATIENTS PROJECT
---

*The code is organized in the following folders:*

	filters/
		Here are the scripts that are used at the beginning to process the VCF to a MAF 
	modules/
		Here are some functions and data (e.g. such as lists and dictionaries of the colors used) recurrently used along the analysis that can be imported. 
		Must be added to you $PYTHONPATH
	processing/
		Here are some jupyter-notebooks that were used in some intermediate step between the filtered MAFs and the final results as a figures
	notebook_figures/
		Here are jupyter-notebooks generating the figures of the paper
	mut_rate_models/
		Here is the code used for the simulations and figures of the mutation rate increment models
	ext_runs/
		Here are the scripts showing how external computational tools and software were run

*Extra folders and files:*
		
	ext_files/
		Files that are necessary to run certain part of the analysis such as files required to run external tools
	intermediate_files/
		Files generated at some point of the analysis that are the result of some intermediate 

*Workspace*

All patients had al least a tumoral and normal sample and, in some cases, two tumoral (primary and relapse) samples and a normal sample. Therefore, the data was store as:

```
cohort/
	patientID/
		DxTumorID_vs_normalID/
		ReTumorID_vs_normalID/ (sometimes)

```

This working directory structure was used in the code			

**STEPS TO REPRODUCE THE ANALYSIS**

The following intructions are in the right order to obtain the same results

ALIGNMENT AND CALLING

1. Sarek pipeline v2.2.1 (intruccions to install and run --> https://github.com/nf-core/sarek)

	Run preprocessing step and Strelka in variantcalling step

2. FACETS (v0.5.6)

    ```ext_runs/run_FACETS/facets_cnv_wgs_hg19.sh```

	the README file within the run_FACETS/ directory shows how to install it 
 
3. DELLY (v.0.7.9)

    ```ext_runs/run_delly/run_delly.sh```

	there is a yml file to install an enviroment with delly for conda

FILTER FROM STRELKA VCFs

The following list of scripts from filters/ shows the order in which the scripts were executed to obtain the MAF files

1. ```filters/process_strelka_v2.9.3_vcf.py```

2. ```filters/refish_shared.py```

3. ```filters/check_for_missed_MNVs_muts.py```

4. ```ext_runs/run_vep/run_vep_cannonical.sh```

5. ```filters/process_vep.py```

6. ```filters/filter_snps_maf.py```

7. ```filters/clonal_classification_maf.py```


Within ```ext_runs/run_vep``` there is a README file with instructions on how to install vep.

The python scripts can be run within a conda environment of python 3.7. There is a yml file to build de environment (environment.yml)

```conda env create -f environment.yml```

PROCESSING

1. Run IntOGen v20191009

	In the FAQ section of IntOGen are the instructions to install it and run it locally https://www.intogen.org/faq 

2. Processing of IntOGen results. Inspect list of driver genes and filter out possible FP

    ```processing/drivers_intogen.ipynb```

3. Get a list of protein affecting mutations (SNVs and InDels) in driver genes
   
		processing/
			driver_mutations_primary_ALL.ipynb
			driver_mutations_TALL.ipynb	
	
  
4. Processing FACETS results, get private and shared copy number changes

    ```processing/cnv_overlapping_segments.ipynb```

5. Annotate CNV. Add cytobands and check for driver gains and losses

		processing/
			annotate_cnv_TALL.ipynb
			driver_cnv_TALL.ipynb

6. Process SV variants and get the driver known ones
	
		processing/
			SV_parser.ipynb
			driver_SV_TALL.ipynb	

7.  To get exons of NOTCH1 the mutations within this gene were re-annotated with VEP 
    
    ```processing/re-annotate_NOTCH1_muts.ipynb```

8. Run deconstructSigs to fit known signatures

	Inside the folder ```ext_runs/run_deconstructSig``` there is a README.txt file on how it was run 

FIGURES AND RESULTS

The jupyter-notebooks within ```notebook_figures/``` generate the figures of the paper separately.

The final figures were layed out with SVG editable software. 

Figure 1

```
notebook_figures/
	number_mutations_cohort.ipynb
	UMAP_mutational_profile.ipynb
	signature_barplots_by_cohort.ipynb
	clustering_driver_combinations.ipynb
```
	
Figure 2

```
notebook_figures/
	table_driver_alterations_TALL.ipynb (also contains supplementary figure)
	NOTCH1_needle_plot.ipynb
	NOTCH1_pathway_change.ipynb
```
	
Figure 3

```
notebook_figures/
	counts_with_probability_exposure_assignments_phylotree.ipynb
	signature_barplots_evolution_part.ipynb
``` 
(```signature_barplots_evolution_part.ipynb``` also contains supplementary figures)
	
Figure 4

```
notebook_figures/mut_rate_test_&_plot.ipynb
``` 
(also contains supplementary figures)

(missing figures are generated with R code in ```mut_rate_models/mutation_models.R```)
		
Figure 5

```
notebook_figures/relapse_evolution_with_doubling_estimates.ipynb
```
(missing figures are generated with R code in ```mut_rate_models/relapse_fixation.R```)

Additional files with supplementary figures:

_Additional 3_
	
```
notebook_figures/
	barplot_recalling_comparative.ipynb
	barplot_snps_filter.ipynb
	ccf_clonality_scatterplots.ipynb
```

_Additional 1_
	
```
notebook_figures/
	table_driver_mutations_primary_ALL.ipynb
	clinical_plot.ipynb
	plot_total_minor_copy_number.ipynb
	relapse_growth_model.ipynb
```


MUTATION RATE MODELS

The code in this section is written in R and can be found in ```mut_rate_models/mutation_models.R```.

The Rscript performs the simulations and the figures with the summary of the results

TUMOR GROWTH SIMULATIONS

The simulations of the tumor growth were performed with Clonex:  https://github.com/gerstung-lab/clonex.

The Rscript in ```mut_rate_models/relapse_fixation.R``` takes the results of Clonex and outputs the figures of the paper
