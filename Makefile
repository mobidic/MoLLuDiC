get-vcf-from-muse:

	scp guignardt@muse-login.hpc-lr.univ-montp2.fr:/home/guignardt/scratch/output/19949/FASTEXOME6/MERGED_SAMPLES/*/MERGED_SAMPLES* ./


MPA:

	python3 /softs/MPA/MPA.py -i FASTEXOME6.hg19_multianno.vcf -o FASTEXOME6.hg19_multianno.MPA.vcf

achab:

	perl /softs/Captain-ACHAB/achab.pl --vcf FASTEXOME6.hg19_multianno.MPA.vcf --case CAD180045-CI --dad CSG180964-pere --mum CSG180965-mere --trio YES --candidats DETRESSE_NEONAT_OMIM.genes.txt --phenolyzerFile out.predicted_gene_scores
