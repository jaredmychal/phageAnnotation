### -------- phageAnnote - Bacteriophage annotation pipeline using conda and snakemake -------- ###


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
PROKKA_HMM_DB = config['hmmdb']
VIRAL_PROTEINS = config['virus_refseq_protein']
NUM_THREADS = 6


# **** Begin script ****

Samples, = glob_wildcards("raw_files/{sample}.fasta")


rule all:
    input:
        expand("results/{sample}_files/gff/{sample}.master.gff3", sample=Samples)


# **** Prokka annotation ****

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
        mkdir -p results
        mkdir -p results/edits
        prokka --kingdom Viruses --gcode 11 --proteins {VIRAL_PROTEINS} --prefix {params.sample_name} --locustag {params.sample_name} --addmrna --addgenes --rfam --rnammer --cdsrnaolap --outdir results {input} --force
        """


#rule blastp:
#    input:
#        rules.prokka.output.r2
#    output:
#        "results/{sample}_blastp_alignment.tbl"
#    conda:
#        "envs/prokka.yaml"
#    params:
#        sample_name = "{sample}"
#    threads: NUM_THREADS
#    shell:
#        """
#        blastp -db /home/js/phageAnnote/viral_proteins/refseq_viral -query {input} -max_target_seqs 1 -outfmt '6 qaccver saccver stitle pident bitscore score' -out {output}
#        """

# **** Search for acquired AMR and virulence genes with Abricate****

rule abricate:
    input:
        r1 = "results/{sample}.gbk"
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
        abricate -db ncbi {input.r1} > {output.r1}
        abricate -db resfinder {input.r1} > {output.r2}
        abricate -db argannot {input.r1} > {output.r3}
        abricate -db vfdb {input.r1} > {output.r4}
        abricate --summary {output.r1} {output.r2} {output.r3} {output.r4} > {output.r5}
        """

# **** Predict AMR genes using NCBI _amrfinder

rule ncbi_amr:
    input:
        "results/{sample}.faa"
    output:
        "results/{sample}_ncbi_amr.txt"
    params:
        sample_name = "{sample}"
    threads: NUM_THREADS
    shell:
        """
        echo Annotating AMR using NCBI-AMRFINDERPLUS...
        amrfinder -p {input} --plus --threads {NUM_THREADS} -o {output}
        """

# **** Annotate prokaryotic virus orthologous group proteins from pVOG database using Hmmer

rule pvog:
    input:
        "results/{sample}.faa"
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
        hmmsearch --domtblout {output.r1} -A {output.r2} -o results/{params.sample_name}_vog.txt {PROKKA_HMM_DB} {input}
        """

# **** GBK2PTT creates a coordinate file from GBK input for use with TransTermHP ****

rule gbk2ptt:
    input:
        "results/{sample}.gbk"
    output:
        "results/{sample}.ptt"
    conda:
        "envs/regulatory.yaml"
    threads: NUM_THREADS
    shell:
        """
        perl scripts/gbk2ptt.pl < {input} > {output}
        """

# **** Predict Rho-independent terminatrs with TransTermHP ****

rule transterm:
    input:
        "results/{sample}.ptt"
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
        transterm -p {miniconda_env}/data/expterm.dat raw_files/{params.sample_name}.fasta {input} --bag-output {output.r1} > {output.r2}
        """

# **** Predict bacterial Sigma 70 promoters with bTSSfinder ****

rule btssfinder:
    input:
        "raw_files/{sample}.fasta"
    output:
        r1 = "results/{sample}_btss.bed"
    conda:
        "envs/regulatory.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        echo Annotating promoters...
        export bTSSfinder_Data="{BTSS_BIN}/Data"
        {BTSS_BIN}/bTSSfinder -i {input} -o results/{params.sample_name}_btss -h 2 -a 1.96 -c 70 -t e
        """
    #    python scripts/get_intergene.py {input.r1}
    #    mv {params.sample_name}_ign.fasta results
    #    export bTSSfinder_Data="{BTSS_BIN}/Data"
# Predict protein families and additional protein functions using  INTERPROSCAN and Eggnog-Mapper

rule interproscan:
    input:
        "results/{sample}.faa"
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
        {INTERPRO_BIN}/interproscan.sh -i {input} --appl TMHMM,pfam,coils --iprlookup --goterms --pathways -d results
        """

rule eggnog:
    input:
        "results/{sample}.faa"
    output:
        "results/{sample}.emapper.annotations"
    conda:
        "envs/interproscan.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        python /home/js/eggnog-mapper/emapper.py -i {input} --output results/{params.sample_name} --cpu 6 -m diamond --cpu 6
        """

# **** Clean files for gff3 format

rule gff_format:
    input:
        r1 = "results/{sample}.gff",
        r2 = "results/{sample}_summary.tbl",
        r3 = "results/{sample}_ncbi_amr.txt",
        r4 = "results/{sample}_alignment.tbl",
        r5 = "results/{sample}_tthp.bag",
        r6 = "results/{sample}_btss.bed",
        r7 = "results/{sample}.faa.gff3",
        r8 = "results/{sample}.emapper.annotations"
    output:
        "results/{sample}_files/gff/{sample}.master.gff3"
    conda:
        "envs/interproscan.yaml"
    threads: NUM_THREADS
    params:
        sample_name = "{sample}"
    shell:
        """
        mkdir -p results/edits
        echo Aligning Sequence Data...

        {GFFREAD_BIN} results/{params.sample_name}.gff -g results/{params.sample_name}.faa --tlf > results/edits/{params.sample_name}.tlf
        sed -i 's,\({params.sample_name}\)\(.*ID=\)\(.*\)\(_mRNA.*\),\\3\\2\\3\\4,g' results/edits/{params.sample_name}.tlf
        cp results/{params.sample_name}.gff results/edits/{params.sample_name}_prokka.edit.gff


        echo Collecting InterProScan5 Data...

        sed '1,3d' results/{params.sample_name}.faa.gff3 | sed 's/##sequence.*//g' | sed '/##FASTA/,$d' | sed 's/protein_match/CDS/; s/date=[[:digit:]]\+\x2d[[:digit:]]\+\x2d[[:digit:]]\+;//; s/Target={params.sample_name}_[[:digit:]]\+[[:space:]][[:digit:]]\+[[:space:]][[:digit:]]\+//; s/ID=match\$[[:digit:]]\+\x5f[[:digit:]]\+\x5f[[:digit:]]\+//; s/status=T;//; s/;;/;/; s/polypeptide/CDS/g' | tr -s '\n' | sort -k 1,1 > results/edits/{params.sample_name}.faa.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}.faa.gff > results/edits/{params.sample_name}_interpro.edit.gff

        cat results/edits/{params.sample_name}_interpro.edit.gff | grep 'TMHMM' > results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 || true
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{$9=\x22;note=Region of a membrane-bound protein predicted to be embedded in the membrane as predicted by TMHMM\x3b\x22; print}}' results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 > results/edits/{params.sample_name}_interpro_tmhmm.edit.2.gff3 && mv results/edits/{params.sample_name}_interpro_tmhmm.edit.2.gff3 results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3
        cat results/edits/{params.sample_name}_interpro.edit.gff | grep 'Pfam' > results/edits/{params.sample_name}_interpro_pfam.edit.gff3 || true
        sed -i 's/signature_desc=/note=Pfam predicted protein /; s/Name=/;note=/g' results/edits/{params.sample_name}_interpro_pfam.edit.gff3
        cat results/edits/{params.sample_name}_interpro.edit.gff | grep 'Coils' > results/edits/{params.sample_name}_interpro_coils.edit.gff3 || true
        sed -i 's/Name=Coil/note=Coiled-coil domain predicted in protein by InterProScan5/g' results/edits/{params.sample_name}_interpro_coils.edit.gff3
        if [[ -s results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 ]]; then awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_interpro_tmhmm.edit.gff3 results/edits/{params.sample_name}_prokka.edit.gff >  results/edits/{params.sample_name}_prokka.tmhmm.edit.gff; else cp results/edits/{params.sample_name}_prokka.edit.gff results/edits/{params.sample_name}_prokka.tmhmm.edit.gff; fi
        if [[ -s results/edits/{params.sample_name}_interpro_pfam.edit.gff3 ]]; then awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_interpro_pfam.edit.gff3 results/edits/{params.sample_name}_prokka.tmhmm.edit.gff >  results/edits/{params.sample_name}_prokka.pfam.edit.gff; else cp results/edits/{params.sample_name}_prokka.tmhmm.edit.gff results/edits/{params.sample_name}_prokka.pfam.edit.gff; fi
        if [[ -s results/edits/{params.sample_name}_interpro_coils.edit.gff3 ]]; then awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_interpro_coils.edit.gff3 results/edits/{params.sample_name}_prokka.pfam.edit.gff >  results/edits/{params.sample_name}_prokka.interpro.edit.gff; else cp results/edits/{params.sample_name}_prokka.pfam.edit.gff results/edits/{params.sample_name}_prokka.interpro.edit.gff; fi


        echo Collecting pVOG Data...

        sed '1,3d' results/{params.sample_name}_alignment.tbl | head -n -11 | tr -s '\n' | awk -F" " -v OFS='\t' '{{print $1, "pVOG", "CDS", "4", "5", "6", "7", "8", ";note=pVOG family " $4 ";"}}' | sort -k 1,1 > results/edits/{params.sample_name}_alignment.tbl
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}_alignment.tbl > results/edits/{params.sample_name}_alignment.edit.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_alignment.edit.gff results/edits/{params.sample_name}_prokka.interpro.edit.gff > results/edits/{params.sample_name}_prokka.pvog.edit.gff


        echo Collecting Eggnog_mapper Data...

        sed '1,4d' results/{params.sample_name}.emapper.annotations | head -n -3 | awk -v FS='\t' -v OFS='\t' '{{print $1, "EGG", "CDS", "4", "5", "6", "7", "8", ";gene="$6";db_xref=KEGG:"$9";note=Eggnog-mapper predicted functional annotation: " $22 ";"}}' > results/edits/{params.sample_name}.emapper.annotations.edit
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}.emapper.annotations.edit > results/edits/{params.sample_name}.emapper.annotations.edit.2
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}.emapper.annotations.edit.2 results/edits/{params.sample_name}_prokka.pvog.edit.gff > results/edits/{params.sample_name}_prokka.eggnog.edit.gff


        echo Collecting AMR and Virulence Data...

        awk -v FS='\t' -v OFS='\t' '{{print $1, "NCBI", "CDS", "4", "5", "6", "7", "8",  ";note="$5" gene" " similar to " $15 ", " $3 " predicted by NCBI AMRFinder Plus;"}}' results/{params.sample_name}_ncbi_amr.txt > results/edits/{params.sample_name}_ncbi_amr.gff
        awk -v FS='\t' -v OFS='\t' '{{print $2, "ABRICATE","CDS", $3, $4, ".", $5, ".", ";note=Virulence factor similar to " $14 " predicated by Abricate v1.0.0;"}}' results/{params.sample_name}_virulence.tbl > results/edits/{params.sample_name}_virulence.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$1]=$4;b[$1]=$5; next}}{{$4=a[$1]; $5=b[$1]; print}}' results/edits/{params.sample_name}.tlf results/edits/{params.sample_name}_ncbi_amr.gff > results/edits/{params.sample_name}_ncbi_amr.edit.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_ncbi_amr.edit.gff results/edits/{params.sample_name}_prokka.eggnog.edit.gff > results/edits/{params.sample_name}_prokka.amr.edit.gff
        awk -v FS='\t' -v OFS='\t' 'NR==FNR{{a[$3$4]=$9; next}}{{$3$4 in a; $9 = $9a[$3$4]; print}}' results/edits/{params.sample_name}_virulence.gff results/edits/{params.sample_name}_prokka.amr.edit.gff > results/edits/{params.sample_name}_prokka.abricate.edit.gff


        echo Collecting Promoter Data...

        grep '##sequence-region' results/{params.sample_name}_btss.gff | sed 's/\([[:digit:]]\+\)\x2d\([[:digit:]]\+\)/\\1 \\2/g' > results/edits/{params.sample_name}_btss.intergene.edit.gff
        cp results/{params.sample_name}_btss.gff results/edits/{params.sample_name}_btss.edit.gff
        sed -i 's/tss\t/Promoter\t/; s/\([[:digit:]]\+\)\x2d\([[:digit:]]\+\)/\t\\1\t\\2\t/; s/promoter=/note=/g' results/edits/{params.sample_name}_btss.edit.gff
        sed -i 's/\(sigma[[:digit:]]\+\)\(.*\)$/\\1 promoter predicted by bTSSfinder v1.206\\2\x3b/g' results/edits/{params.sample_name}_btss.edit.gff
        awk -F'\t' -v OFS='\t' '{{$4=$4-35; $5=$5-10}} {{print}}' results/edits/{params.sample_name}_btss.edit.gff > results/edits/{params.sample_name}_btss.final.gff
        sed -i 's/box10seq=/box10seq:/; s/box35seq=/box35seq:/; s/box10pos=/box10pos:/; s/box35pos=/box35pos:/; s/;box10pos/,box10pos/; s/;box10seq/,box10seq/; s/;box35seq/,box35seq/; s/;box35pos/,box35pos/g' results/edits/{params.sample_name}_btss.final.gff
        sed -i 's/.*btssfinder\tPromoter/{params.sample_name}\tbtssfinder\tPromoter/g' results/edits/{params.sample_name}_btss.final.gff
        sed -i '1,8d' results/edits/{params.sample_name}_btss.final.gff
        grep 'btssfinder' results/edits/{params.sample_name}_btss.final.gff > results/edits/{params.sample_name}_btss.final.gff3


        echo Collecting Terminator Data...

        sed 's/.*{params.sample_name}_[0-9]\{{5\}}/{params.sample_name}\tTranstermHP\tTerminator\t/; s/.*NONE.*//g' results/{params.sample_name}_tthp.bag > results/edits/{params.sample_name}_tthp.edit.bag
        tr -s '\n' < results/edits/{params.sample_name}_tthp.edit.bag > results/edits/{params.sample_name}_tthp.edit.gff
        sed -i 's,\({params.sample_name}\tTranstermHP\tTerminator\t\)[[:blank:]]\+\([[:digit:]]\+\)[[:blank:]]\+\x2e\x2e[[:blank:]]\+\([[:digit:]]\+\)[[:blank:]]\+\([\x2b\x2d]\)[[:blank:]]\+[\x2d\x2b].*[\x2d\x2b][[:digit:]]\+\x2e[[:digit:]]\+[[:blank:]]\+\([[:upper:]]\+\)[[:blank:]]\+\(.*[[:blank:]]\+.*[[:blank:]]\+.*\)[[:blank:]]\+\([[:upper:]]\+\)[[:blank:]]\+[[:digit:]]\+[[:blank:]]\+[[:digit:]]\+$,\\1\\2\t\\3\t\x2e\t\\4\t\x2e\tnote=Rho-Independent Terminator predicted by TranstermHP v2.08\x3b,g' results/edits/{params.sample_name}_tthp.edit.gff
        grep '{params.sample_name}' results/edits/{params.sample_name}_tthp.edit.gff > results/edits/{params.sample_name}_tthp.final.gff3


        echo Cleaning up empty tags...

        sed 's,\(inference=ab initio prediction:Prodigal:002006.*sequence:.*:\)\(.*\)\(;locus_tag=.*product=\)\(.*\)\s\[\(.*\)\],\\3\\4;note=Similar to \\5 \\4 \\2,g' results/edits/{params.sample_name}_prokka.abricate.edit.gff > results/edits/{params.sample_name}_prokka.edit.2.gff
        sed -i 's/\t\x3b/\x3b/; s/inference=ab initio prediction:Prodigal:002006;//; s/.*\tmRNA\t.*//; s/signature_desc=/note=/; s/[[:alpha:]]\+=;//; s/note=Eggnog-mapper predicted functional annotation: ;//; s/db_xref=KEGG:;//; s/;;/;/g' results/edits/{params.sample_name}_prokka.edit.2.gff
        tr -s '\n' < results/edits/{params.sample_name}_prokka.edit.2.gff > results/edits/{params.sample_name}_prokka.edit.3.gff
        tr -s ';' < results/edits/{params.sample_name}_prokka.edit.3.gff > results/edits/{params.sample_name}_prokka.final.gff


        echo Moving files...

        mkdir -p finished_genomes
        mkdir -p results/{params.sample_name}_files
        mkdir -p results/{params.sample_name}_files/gff
        cp results/edits/{params.sample_name}_tthp.final.gff3 results/{params.sample_name}_files/gff
        cp results/edits/{params.sample_name}_btss.final.gff3 results/{params.sample_name}_files/gff
        cp results/edits/{params.sample_name}_prokka.final.gff results/{params.sample_name}_files/gff

        echo Generating GFF3 master file...
        cat results/{params.sample_name}_files/gff/{params.sample_name}*.gff3 results/{params.sample_name}_files/gff/{params.sample_name}_prokka.final.gff > results/{params.sample_name}_files/gff/{params.sample_name}.master.gff3
        mv results/*{params.sample_name}*.* results/{params.sample_name}_files
        mv raw_files/{params.sample_name}.fasta finished_genomes

        echo Removing temporary files...
        rm -r results/edits
        echo Done!
        """
