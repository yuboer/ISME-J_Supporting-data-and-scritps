# ----Step 1. Quality control of the raw/clean illumina .fq sequences
## 1.1 Reads filtration
fastp -i EBPR420_R1.fq.gz -I EBPR420_R2.fq.gz -o ./fastp_out/clean_EBPR420_R1.fq.gz -O ./fastp_out/clean_EBPR420_R2.fq.gz --detect_adapter_for_pe --dedup --trim_poly_g --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --length_required 50 --thread 2 --html EBPR420_fastp_report.html --json EBPR420_fastp_report.json
fastp -i EBPR617_R1.fq.gz -I EBPR617_R2.fq.gz -o ./fastp_out/clean_EBPR617_R1.fq.gz -O ./fastp_out/clean_EBPR617_R2.fq.gz --detect_adapter_for_pe --dedup --trim_poly_g --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --length_required 50 --thread 2 --html EBPR617_fastp_report.html --json EBPR617_fastp_report.json
fastp -i EBPR723_R1.fq.gz -I EBPR723_R2.fq.gz -o ./fastp_out/clean_EBPR723_R1.fq.gz -O ./fastp_out/clean_EBPR723_R2.fq.gz --detect_adapter_for_pe --dedup --trim_poly_g --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --length_required 50 --thread 2 --html EBPR723_fastp_report.html --json EBPR723_fastp_report.json
## 1.2 Quality check and visualization:
fastqc *.gz -t 6 -o fastaqc_out
multiqc ./fastaqc_out/
# ----Step 2. Assembly of clean reads to contigs
## 2.1 de novo assembly using SPAdes (v4.2.0) in metagenomic mode (metaSPAdes). 
spades.py --pe1-1 clean_EBPR420_R1.fq.gz --pe1-2 clean_EBPR420_R2.fq.gz --meta --only-assembler -k 21,33,55,77,99 -o e4_spades_output
spades.py --pe1-1 clean_EBPR617_R1.fq.gz --pe1-2 clean_EBPR617_R2.fq.gz --meta --only-assembler -k 21,33,55,77,99 -o e6_spades_output
spades.py --pe1-1 clean_EBPR723_R1.fq.gz --pe1-2 clean_EBPR723_R1.fq.gz --meta --only-assembler -k 21,33,55,77,99 -o e7_spades_output
## 2.2 co-assembly using MEGAHIT (v1.2.9)
megahit -1 clean_EBPR617_R1.fq.gz,clean_EBPR723_R1.fq.gz -2 clean_EBPR617_R2.fq.gz,clean_EBPR723_R1.fq.gz --k-list 21,33,55,77,99 -o eco2_output
megahit -1 clean_EBPR420_R1.fq.gz,clean_EBPR617_R1.fq.gz,clean_EBPR723_R1.fq.gz -2 clean_EBPR420_R2.fq.gz,clean_EBPR617_R2.fq.gz,clean_EBPR723_R1.fq.gz --k-list 21,33,55,77,99 -o eco3_output
# ----Step 3. Bowtie mapping + Samtool alignment
## 3.1 build the index file for mapping
bowtie2-build e4_contigs.fasta e4_mapping.contigs > bowtie2_build_e4.log 2>&1
bowtie2-build e6_contigs.fasta e6_mapping.contigs > bowtie2_build_e6.log 2>&1
bowtie2-build e7_contigs.fasta e7_mapping.contigs > bowtie2_build_e7.log 2>&1
bowtie2-build eco2_final.contigs.fa eco2_mapping.contigs > bowtie2_build_eco2.log 2>&1
bowtie2-build eco3_final.contigs.fa eco3_mapping.contigs > bowtie2_build_eco3.log 2>&1
## 3.2 Bowtie mapping
bowtie2 --threads 12 -x e4_mapping.contigs -1 clean_EBPR420_R1.fq.gz -2 clean_EBPR420_R2.fq.gz -S e4_e4.bowtie2.sam > e4_e4_mapping.log 2>&1
bowtie2 --threads 12 -x e4_mapping.contigs -1 clean_EBPR617_R1.fq.gz -2 clean_EBPR617_R2.fq.gz -S e4_e6.bowtie2.sam > e4_e6_mapping.log 2>&1
bowtie2 --threads 12 -x e4_mapping.contigs -1 clean_EBPR723_R1.fq.gz -2 clean_EBPR723_R2.fq.gz -S e4_e7.bowtie2.sam > e4_e7_mapping.log 2>&1

bowtie2 --threads 12 -x e6_mapping.contigs -1 clean_EBPR420_R1.fq.gz -2 clean_EBPR420_R2.fq.gz -S e6_e4.bowtie2.sam > e6_e4_mapping.log 2>&1
bowtie2 --threads 12 -x e6_mapping.contigs -1 clean_EBPR617_R1.fq.gz -2 clean_EBPR617_R2.fq.gz -S e6_e6.bowtie2.sam > e6_e6_mapping.log 2>&1
bowtie2 --threads 12 -x e6_mapping.contigs -1 clean_EBPR723_R1.fq.gz -2 clean_EBPR723_R2.fq.gz -S e6_e7.bowtie2.sam > e6_e7_mapping.log 2>&1

bowtie2 --threads 12 -x e7_mapping.contigs -1 clean_EBPR420_R1.fq.gz -2 clean_EBPR420_R2.fq.gz -S e7_e4.bowtie2.sam > e7_e4_mapping.log 2>&1
bowtie2 --threads 12 -x e7_mapping.contigs -1 clean_EBPR617_R1.fq.gz -2 clean_EBPR617_R2.fq.gz -S e7_e6.bowtie2.sam > e7_e6_mapping.log 2>&1
bowtie2 --threads 12 -x e7_mapping.contigs -1 clean_EBPR723_R1.fq.gz -2 clean_EBPR723_R2.fq.gz -S e7_e7.bowtie2.sam > e7_e7_mapping.log 2>&1

bowtie2 --threads 12 -x eco2_mapping.contigs -1 clean_EBPR420_R1.fq.gz -2 clean_EBPR420_R2.fq.gz -S eco2_e4.bowtie2.sam > eco2_e4_mapping.log 2>&1
bowtie2 --threads 12 -x eco2_mapping.contigs -1 clean_EBPR617_R1.fq.gz -2 clean_EBPR617_R2.fq.gz -S eco2_e6.bowtie2.sam > eco2_e6_mapping.log 2>&1
bowtie2 --threads 12 -x eco2_mapping.contigs -1 clean_EBPR723_R1.fq.gz -2 clean_EBPR723_R2.fq.gz -S eco2_e7.bowtie2.sam > eco2_e7_mapping.log 2>&1

bowtie2 --threads 12 -x eco3_mapping.contigs -1 clean_EBPR420_R1.fq.gz -2 clean_EBPR420_R2.fq.gz -S eco3_e4.bowtie2.sam > eco3_e4_mapping.log 2>&1
bowtie2 --threads 12 -x eco3_mapping.contigs -1 clean_EBPR617_R1.fq.gz -2 clean_EBPR617_R2.fq.gz -S eco3_e6.bowtie2.sam > eco3_e6_mapping.log 2>&1
bowtie2 --threads 12 -x eco3_mapping.contigs -1 clean_EBPR723_R1.fq.gz -2 clean_EBPR723_R2.fq.gz -S eco3_e7.bowtie2.sam > eco3_e7_mapping.log 2>&1
## 3.3 Samtool alignment with the .bowtie2.sam files generated from step 3.2
samtools view -bS $BASENAME.bowtie2.sam > $BASENAME.bam
samtools sort $BASENAME.bam -o $BASENAME.sorted.bam
# ----Step 4. Genome binning, genome quality control, and taxonomy annotation
## 4.1 generate depth file from .sorted.bam generated from step 3.3
jgi_summarize_bam_contig_depths --outputDepth e4_depth.txt ./*.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth e6_depth.txt ./*.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth e7_depth.txt ./*.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth eco2_depth.txt ./*.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth eco3_depth.txt ./*.sorted.bam
## 4.2 metabat2 binning
metabat2 -i e4_scaffolds_2.5k.fasta -a e4_depth.txt -o e4_bins/e4_bin
metabat2 -i e6_scaffolds_2.5k.fasta -a e6_depth.txt -o e6_bins/e6_bin
metabat2 -i e7_scaffolds_2.5k.fasta -a e7_depth.txt -o e7_bins/e7_bin
metabat2 -i eco2_final.contigs_2.5k.fasta -a eco2_depth.txt -o eco2_bins/eco2_bin
metabat2 -i eco3_final.contigs_2.5k.fasta -a eco3_depth.txt -o eco3_bins/eco3_bin
## 4.3 combine the genome bins generated and dereplicate
dRep dereplicate drep_bin_ebpr/ -g ./all_bins/*.fa
## 4.4 genome quality check
checkm lineage_wf ./dereplicated_e_genomes/ ./checkm_e_bins/ -x fa --reduced_tree
## 4.5 genome taxonomy annotation
gtdbtk classify_wf -x fa --cpus 8 --genome_dir ./dereplicated_e_genomes --out_dir e_drep_gtdb

# ----Step 5. Genome Annotation and Consensus Pathway Reconstruction
## 5.1 prokka
prokka $BinID.fa --prefix $BinID --kingdom Bacteria --gcode 11 --force --rfam --addgenes --addmrna --cpus 1 --outdir ../prokka/$BinID

## 5.2 bakta
for genome in /dereplicated_genomes/*.fa; do
    base=$(basename ${genome} .fa)
    echo "=== Running Bakta on ${base} ==="
    bakta \
      --db /bakta_db/db \
      --output /bakta_out/${base} \
      --prefix ${base} \
      --threads 4 \
      --keep-contig-headers \
      ${genome}
done
echo "✅ All genomes processed!"

## 5.3 metabolic annoation of the .faa files yielded from prokka and/or bakta annoation
perl METABOLIC-C.pl \
  -in /folder_faa \
  -o /metabolic_out \
  -t 4 \
  -m-cutoff 0.75 \
  -kofam-db full

echo "🎉 METABOLIC run completed successfully!"

## 5.4 eggnog annotation of the .faa files yielded from prokka and/or bakta annoation
export EGGNOG_DATA_DIR=/scratch/ywang116/software/eggnog_db
INPUT_DIR="/scratch/ywang116/bakta_faa"
OUTPUT_DIR="/scratch/ywang116/emapper_out"
for faa in /folder_faa/*.faa; do
    base=$(basename "${faa}" .faa)
    echo "Processing ${base}..."
    emapper.py \
        -i "${faa}" \
        --itype proteins \
        --cpu 4 \
        --data_dir "${EGGNOG_DATA_DIR}" \
        --output "${base}" \
        --output_dir "${OUTPUT_DIR}" \
        --sensmode more-sensitive \
        --tax_scope auto_broad \
        --excel \
        --override \
        > ${OUTPUT_DIR}/${base}.log 2>&1
done

# ---- Step 6. Phylogenetic tree construction
# 6.1 muscle alignment
muscle -align gtdbtk.bac120_msa.fasta -output ./muscle.fasta
# 6.2 If want to keep the most informative sites without being too aggressive, a common "middle-ground" approach with trimal is using the -automated1 flag, which selects the best procedure based on the alignment characteristics:
trimal -in muscle.fasta -out trimmed.fasta -automated1
# 6.3 fasttree construction
fasttree -lg -boot 1000 trimmed.fasta > tree.nwk

