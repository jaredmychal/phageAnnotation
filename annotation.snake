
# **** Imports ****
import os
import glob
# **** Variables ****

configfile: "config.yaml"
sequence_directory = config['fast_files']
miniconda_env = config['conda_env']
INTERPRO_BIN = config['interproscan']
BTSS_BIN = config['btssfinder']
GFFREAD_BIN = config['gffread']
REFERENCE_GENOME = config['ref_genome']
PROKKA_HMM_DB = config['hmmdb']
NUM_THREADS = 6

# **** Begin script ****

Samples, = glob_wildcards("raw_files/{sample}.fasta")


rule all:
    input:
        expand("results/gff/{sample}.master.gff3", sample=Samples)

# **** Prokka runs annotation ****

rule prokka:
    input:
        "raw_files/{sample}.fasta"
    output:
        r1 = "results/{sample}.gbk",
        r2 = "results/{sample}.faa",
        r3 = "results/{sample}.fna"
    conda:
        "envs/prokka.yaml"
    params:
        sample_name = "{sample}"
    threads: NUM_THREADS
    shell:
        """
        echo Annotating...
        prokka --kingdom Viruses --gcode 11 --prefix {params.sample_name} --locustag {params.sample_name} --addmrna --addgenes --rfam --rnammer --cdsrnaolap --outdir results {input} --force
        """

# **** Search for acquired AMR and virulence genes with Abricate****

rule abricate:
    input:
        rules.prokka.output.r1
    output:
        r1 = "results/{sample}_amr.tbl",
        r2 = "results/{sample}_amr_resfinder.tbl",
        r3 = "results/{sample}_amr_argannot.tbl",
        r4 = "results/{sample}_virulence.tbl",
        r5 = "results/{sample}_summary.tbl"
    conda:
        "envs/prokka.yaml"
    params:
        sample_name = "{sample}"
    threads: NUM_THREADS
    shell:
        """
        echo Annotating Virulence and AMR...
        abricate -db ncbi {input} > {output.r1}
        abricate -db resfinder {input} > {output.r2}
        abricate -db argannot {input} > {output.r3}
        abricate -db vfdb {input} > {output.r4}
        abricate --summary {output.r1} {output.r2} {output.r3} {output.r4} > {output.r5}
        """

# **** NCBI _amrfinder

rule ncbi_amr:
    input:
        r1 = rules.prokka.output.r2,
        r2 = rules.abricate.output
    output:
        "results/{sample}_ncbi_amr.txt"
    params:
        sample_name = "{sample}"
    threads: NUM_THREADS
    shell:
        """
        echo Annotating AMR using NCBI-AMRFINDERPLUS...
        amrfinder -p {input.r1} --plus --threads {NUM_THREADS} -o {output}
        """

# **** Annotate from pVOG database

rule pvog:
    input:
        r1 = rules.ncbi_amr.output,
        r2 = rules.prokka.output.r2
    output:
        r1 = "results/{sample}_alignment.tbl",
        r2 = "results/{sample}_alignment.align"
    conda:
        "envs/hmmer.yaml"
    params:
        sample_name = "{sample}"
    threads: NUM_THREADS
    shell:
        """
        echo Begginning hmmsearch for {input} with pVOG
        hmmsearch --domtblout {output.r1} -A {output.r2} -o results/{params.sample_name} {PROKKA_HMM_DB} {input.r2}
        """

# **** GBK2PTT creates a coordinate file from GBK input for use with TransTermHP ****

rule gbk2ptt:
    input:
        r1 = rules.prokka.output.r1,
        r2 = rules.pvog.output
    output:
        "results/{sample}.ptt"
    conda:
        "envs/regulatory.yaml"
    threads: NUM_THREADS
    shell:
        """
        perl scripts/gbk2ptt.pl < {input.r1} > {output}
        """

# **** Predict Rho-independent terminatrs with TransTermHP ****

rule transterm:
    input:
        r1 = rules.gbk2ptt.output
    output:
        r1 = "results/{sample}_tthp.bag",
        r2 = "results/{sample}_tthp.tt"
    conda:
        "envs/regulatory.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        echo Annotating Terminators...
        transterm -p {miniconda_env}/data/expterm.dat raw_files/{params.sample_name}.fasta {input.r1} --bag-output {output.r1} > {output.r2}
        """

# **** Predict bacterial Sigma promoters with bTSSfinder ****

rule btssfinder:
    input:
        r1 = rules.prokka.output.r1,
        r2 = "raw_files/{sample}.fasta",
        r3 = rules.transterm.output
    output:
        r1 = "results/{sample}_btss.bed"
    conda:
        "envs/regulatory.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        python scripts/get_intergene.py {input.r1}
        mv {params.sample_name}_ign.fasta results
        export bTSSfinder_Data="{BTSS_BIN}/Data"
        echo Annotating promoters...
        {BTSS_BIN}/bTSSfinder -i results/{params.sample_name}_ign.fasta -o results/{params.sample_name}_btss -h 2 -a 1.96 -c 70 -t e
        """

# run INTERPROSCAN

rule interproscan:
    input:
        r1 = rules.prokka.output.r2,
        r2 = rules.btssfinder.output
    output:
        "results/{sample}.faa.gff3"
    conda:
        "envs/interproscan.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        echo Annotating protein functions...
        {INTERPRO_BIN}/interproscan.sh -i {input.r1} --appl SignalP_GRAM_Negative,TMHMM --iprlookup --goterms --pathways -d results
        """

rule eggnog:
    input:
        r1 = rules.prokka.output.r2,
        r2 = rules.interproscan.output
    output:
        "results/{sample}.emapper.annotations"
    conda:
        "envs/interproscan.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        python /home/js/eggnog-mapper/emapper.py -i {input.r1} --output results/{params.sample_name} -m diamond --cpu 6
        """

# **** Clean files for gff3 format

rule gff_format:
    input:
        rules.eggnog.output
    output:
        "results/gff/{sample}.master.gff3"
    conda:
        "envs/interproscan.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        echo Aligning Sequence Data...

        {GFFREAD_BIN} results/{params.sample_name}.gff -g results/{params.sample_name}.faa --tlf > results/edits/{params.sample_name}.tlf
        sed -i 's,\({params.sample_name}\)\(.*ID=\)\(.*\)\(_mRNA.*\),\\3\\2\\3\\4,g' results/edits/{params.sample_name}.tlf
        cp results/{params.sample_name}.gff results/edits/{params.sample_name}_prokka.edit.gff
        sed -i 's/inference.*//g' results/edits/{params.sample_name}_prokka.edit.gff

        echo Adding BLASTP results
        cp results/blast/{params.sample_name}_blastp.txt results/edits/{params.sample_name}_blastp.txt
        tr -d '\n' < results/edits/{params.sample_name}_blastp.txt > results/edits/{params.sample_name}_blastp.2.txt
        sed 's/Query=/\\nQuery=/g' results/edits/{params.sample_name}_blastp.2.txt > results/edits/{params.sample_name}_blastp.txt
        rm results/edits/{params.sample_name}_blastp.2.txt
        sed -i 's/Query= \({params.sample_name}_[[:digit:]]\+\)/\\1\tBLASTP\tCDS\tStart\tEnd\t;ID=\\1;locus_tag=\\1;note=Similar to /g' results/edits/{params.sample_name}_blastp.txt
        sed -i 's/Similar to >YP/Similar to YP/; s/Similar to >NP/Similar to NP/g' results/edits/{params.sample_name}_blastp.txt
        sed -i 's/\]>\(YP.*\)/\] and less similar to \\1;/; s/\]>\(NP.*\)/\] and less similar to \\1;/g' results/edits/{params.sample_name}_blastp.txt
        sed -i 's/[[:alpha:]]>\(YP.*\)/\] and less similar to \\1;/; s/[[:alpha:]]>\(NP.*\)/\] and less similar to \\1;/g' results/edits/{params.sample_name}_blastp.txt
        sed -i 's/\(note=Similar to [[:alpha:]]P_[[:digit:]]\+\.[[:digit:]] \)\(.*\)\( \[.*\)\( and less\)/product=\\2;\\1\\2\\3\\4/g' results/edits/{params.sample_name}_blastp.txt
        sed -i 's/product=.*hypothetical protein.*;note/product=hypothetical protein;note/; s/product=protein of unknown function/product=hypothetical protein/; s/product=unnamed protein product/product=hypothetical protein/; s/putative//g' results/edits/{params.sample_name}_blastp.txt
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}_blastp.txt > results/edits/{params.sample_name}_blastp.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$6; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_blastp.gff results/edits/{params.sample_name}_prokka.edit.gff > results/edits/{params.sample_name}_prokka.blastp.gff

        echo Collecting InterProScan5 Data...

        cp -i results/{params.sample_name}.faa.gff3 results/edits/{params.sample_name}.faa.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}.faa.gff > results/edits/{params.sample_name}_interpro.edit.gff

        sed -i '/##FASTA/,$d' results/edits/{params.sample_name}_interpro.edit.gff
        sed -i 's/{params.sample_name}_[[:digit:]]\+[[:blank:]]\+/{params.sample_name}\t/; s/##sequence-region.*//; s/ID=match$.*;signature/signature/; s/ID=match$.*;Name/Name/; s/Ontology_term=/note=Ontology term:/; s/Target={params.sample_name}_[[:digit:]]\+[[:blank:]]\+[[:digit:]]\+[[:blank:]]\+[[:digit:]]\+//; s/;status=T//; s/Name/note/; s/date\x3d[[:digit:]]\+\x2d[[:digit:]]\+\x2d[[:digit:]]\+\x3b//g' results/edits/{params.sample_name}_interpro.edit.gff
        sed -i '/{params.sample_name}\t\x2e\tpolypeptide/d' results/edits/{params.sample_name}_interpro.edit.gff
        tr -s '\n' < results/edits/{params.sample_name}_interpro.edit.gff > results/edits/{params.sample_name}_interpro.final.gff3
        sed -i 's,\({params.sample_name}\t\)\([[:alpha:]]\+\)\(\t.*\),\\1\\2\\3\x3b,g' results/edits/{params.sample_name}_interpro.final.gff3

        sed -i 's/protein_match/CDS/g' results/edits/{params.sample_name}_interpro.final.gff3
        cat results/edits/{params.sample_name}_interpro.final.gff3 | grep 'SignalP_GRAM_' > results/edits/{params.sample_name}_interpro_signalp.edit.gff3 || true
        sed -i 's/SignalP-noTM/Signal Peptide, no TM, predicted by SignalP4.1;/; s/SignalP-TM/Signal Peptide predicted by SignalP4.1;/g' results/edits/{params.sample_name}_interpro_signalp.edit.gff3
        cat results/edits/{params.sample_name}_interpro.final.gff3 | grep 'TMHMM' > results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 || true
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{$9=\x22;note=Region of a membrane-bound protein predicted to be embedded in the membrane as predicted by TMHMM\x3b\x22; print}}' results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 > results/edits/{params.sample_name}_interpro_tmhmm.edit.2.gff3 && mv results/edits/{params.sample_name}_interpro_tmhmm.edit.2.gff3 results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3
        sort -k 4,4 -n -o results/edits/{params.sample_name}_interpro_signalp.edit.gff3 results/edits/{params.sample_name}_interpro_signalp.edit.gff3
        sort -k 4,4 -n -o results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3
        if [[ -s results/edits/{params.sample_name}_interpro_signalp.edit.gff3 ]]; then awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_interpro_signalp.edit.gff3 results/edits/{params.sample_name}_prokka.blastp.gff > results/edits/{params.sample_name}_prokka.signalp.edit.gff; else cp results/edits/{params.sample_name}_prokka.blastp.gff results/edits/{params.sample_name}_prokka.signalp.edit.gff; fi
        if [[ -s results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 ]]; then awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 results/edits/{params.sample_name}_prokka.signalp.edit.gff > results/edits/{params.sample_name}_prokka.interpro.edit.gff; else cp results/edits/{params.sample_name}_prokka.signalp.edit.gff results/edits/{params.sample_name}_prokka.interpro.edit.gff; fi

        echo Collecting pVOG Data...
        sed 's/[[:blank:]]\+/\t/g' results/{params.sample_name}_alignment.tbl > results/edits/{params.sample_name}_alignment.edit.tbl
        cut -f 1,4,23- results/edits/{params.sample_name}_alignment.edit.tbl > results/edits/{params.sample_name}_alignment.edit.2.tbl
        sort -k 1,1 results/edits/{params.sample_name}_alignment.edit.2.tbl > results/edits/{params.sample_name}_alignment.edit.3.tbl
        sed -i '1,13d' results/edits/{params.sample_name}_alignment.edit.3.tbl
        sed -i 's/\t/ /g' results/edits/{params.sample_name}_alignment.edit.3.tbl
        sed -i 's/\({params.sample_name}_[[:digit:]]\+\)[[:space:]]\(VOG[[:digit:]]\+\)[[:space:]]\(.*\)/\\1\tpVOG\tCDS\tStart\tEnd\t;note=pVOG family \\2: \\3;/g' results/edits/{params.sample_name}_alignment.edit.3.tbl
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}_alignment.edit.3.tbl > results/edits/{params.sample_name}_alignment.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$6; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_alignment.gff results/edits/{params.sample_name}_prokka.interpro.edit.gff > results/edits/{params.sample_name}_prokka.pvog.gff

        echo Collecting Promoter Data...

        cp results/{params.sample_name}_btss.gff results/edits/{params.sample_name}_btss.edit.gff
        sed -i 's/tss\t/Promoter\t/; s/\([[:digit:]]\+\)\x2d\([[:digit:]]\+\)/\t\\1\t\\2\t/; s/promoter=/note=/g' results/edits/{params.sample_name}_btss.edit.gff
        sed -i '1,8d' results/edits/{params.sample_name}_btss.edit.gff
        awk -F'\t' -v OFS='\t' '{{$7=$7+$2-35; $8=$8+$2-10}} {{print}}' results/edits/{params.sample_name}_btss.edit.gff > results/edits/{params.sample_name}_btss.edit.2.gff
        sed -i 's/.*btssfinder\tPromoter/{params.sample_name}\tbtssfinder\tPromoter/g' results/edits/{params.sample_name}_btss.edit.2.gff
        grep 'btssfinder\tPromoter' results/edits/{params.sample_name}_btss.edit.2.gff > results/edits/{params.sample_name}_btss.final.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{$9=\x22;note=Sigma70 promoter predicted by bTSSfinder v1.206\x3b\x22; print}}' results/edits/{params.sample_name}_btss.final.gff > results/edits/{params.sample_name}_btss.final.gff3

        echo Collecting Terminator Data...

        sed 's/.*{params.sample_name}_[0-9]\{{5\}}/{params.sample_name}\tTranstermHP\tTerminator\t/; s/.*NONE.*//g' results/{params.sample_name}_tthp.bag > results/edits/{params.sample_name}_tthp.edit.bag
        tr -s '\n' < results/edits/{params.sample_name}_tthp.edit.bag > results/edits/{params.sample_name}_tthp.edit.gff
        sed -i 's,\({params.sample_name}\tTranstermHP\tTerminator\t\)[[:blank:]]\+\([[:digit:]]\+\)[[:blank:]]\+\x2e\x2e[[:blank:]]\+\([[:digit:]]\+\)[[:blank:]]\+\([\x2b\x2d]\)[[:blank:]]\+[\x2d\x2b].*[\x2d\x2b][[:digit:]]\+\x2e[[:digit:]]\+[[:blank:]]\+\([[:upper:]]\+\)[[:blank:]]\+\(.*[[:blank:]]\+.*[[:blank:]]\+.*\)[[:blank:]]\+\([[:upper:]]\+\)[[:blank:]]\+[[:digit:]]\+[[:blank:]]\+[[:digit:]]\+$,\\1\\2\t\\3\t\x2e\t\\4\t\x2e\tnote=Rho-Independent Terminator predicted by TranstermHP v2.08\x3b,g' results/edits/{params.sample_name}_tthp.edit.gff
        grep '{params.sample_name}' results/edits/{params.sample_name}_tthp.edit.gff > results/edits/{params.sample_name}_tthp.final.gff3

        sed -i 's/\t\x3b/\x3b/; s/.*\tmRNA\t.*//g' results/edits/{params.sample_name}_prokka.pvog.gff
        sed -i 's/signature_desc=/note=/; s/;;/;/g' results/edits/{params.sample_name}_prokka.pvog.gff
        tr -s '\n' < results/edits/{params.sample_name}_prokka.pvog.gff > results/edits/{params.sample_name}_prokka.gff
        tr -s ';' < results/edits/{params.sample_name}_prokka.gff > results/edits/{params.sample_name}_prokka.final.gff

        cp results/edits/{params.sample_name}_tthp.final.gff3 results/gff
        cp results/edits/{params.sample_name}_btss.final.gff3 results/gff
        cp results/edits/{params.sample_name}_prokka.final.gff results/gff


        cat results/gff/{params.sample_name}*.gff3 results/gff/{params.sample_name}_prokka.final.gff > results/gff/{params.sample_name}.master.gff3
        rm results/edits/*

        """
