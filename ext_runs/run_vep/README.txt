1. Download cache files (corresponding to VEP 92 and genome build GRCh37) 
	to use VEP more efficiently (https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#pre)
	
2. Install vep 
	conda create -n vep92 -c bioconda ensembl-vep=92.4

3. Download file to run Gnomad along with VEP
	ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
