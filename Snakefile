import glob

SAMPLES = {}
SAMPLES['ERR1036235'] = {'reference':'ERS764956'}

max_threads = 64

rule all:
    input:
        "results/figures/figure_quast_NG50_vs_misassemblies.pdf",
        "results/figures/figure_quast_completeness_vs_duplication.pdf",
        "results/figures/figure_cpu_time_vs_memory.pdf",
        "results/figures/supplementary/figure_comp_res_gapless.pdf",
        "results/figures/supplementary/figure_coverage_cdf1.pdf",
        "results/figures/supplementary/figure_coverage_cdf2.pdf"

rule figure_comp_res:
    input:
        "results/csv/timelog.csv",
        "input/csv/assemblies.csv"
    output:
        "results/figures/figure_cpu_time_vs_memory.pdf",
        "results/figures/figure_elapsed_time_vs_memory.pdf"
    shell:
        "Rscript bin/figure_comp_res.R {input} {output}"

rule timelog_csv:
    input:
        expand("storage/assemblies/ecoli/00-flye-pacbio_clr_{coverage}x/timelog.txt", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/00-flye-pacbio_clr_{coverage}x/timelog.txt", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/00-flye-hifi_{coverage}x/timelog.txt", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/00-flye-nanopore_{coverage}x/timelog.txt", coverage=[61,30,15,8]),
        expand("storage/assemblies/ecoli/01-gapless-pacbio_clr_{coverage}x-shassembly-SRR3191692/timelog.txt", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/01-gapless-pacbio_clr_{coverage}x-supernova-chromium/timelog.txt", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/01-gapless-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/01-gapless-nanopore_{coverage}x-flye-hifi_33x/timelog.txt", coverage=[61,30,15,8]),
        expand("storage/assemblies/human/01-gapless-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[61,30,15,8]),
        expand("storage/assemblies/ecoli/01-lrscaf-pacbio_clr_{coverage}x-shassembly-SRR3191692/timelog{logs}.txt", coverage=[113,57,28,14,7], logs=['_minimap2','_lrscaf']),
        expand("storage/assemblies/truncatus/01-lrscaf-pacbio_clr_{coverage}x-supernova-chromium/timelog{logs}.txt", coverage=[86,43,21,11,5], logs=['_minimap2','_lrscaf']),
        expand("storage/assemblies/human/01-lrscaf-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/timelog{logs}.txt", coverage=[33,16,8,4], logs=['_minimap2','_lrscaf']),
        expand("storage/assemblies/human/01-lrscaf-nanopore_{coverage}x-flye-hifi_33x/timelog{logs}.txt", coverage=[121,61,30,15,8], logs=['_minimap2','_lrscaf']),
        expand("storage/assemblies/human/01-lrscaf-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/timelog{logs}.txt", coverage=[121,61,30,15,8], logs=['_minimap2','_lrscaf']),
        expand("storage/assemblies/ecoli/02-lr_gapcloser-lrscaf-pacbio_clr_{coverage}x-shassembly-SRR3191692/timelog.txt", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/02-lr_gapcloser-lrscaf-pacbio_clr_{coverage}x-supernova-chromium/timelog.txt", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-nanopore_{coverage}x-flye-hifi_33x/timelog.txt", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/ecoli/02-lr_gapcloser-lrscaf-pacbio_clr_{coverage}x-shassembly-SRR3191692/timelog.txt", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/02-tgs_gapcloser-lrscaf-pacbio_clr_{coverage}x-supernova-chromium/timelog.txt", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-nanopore_{coverage}x-flye-hifi_33x/timelog.txt", coverage=[8]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/timelog.txt", coverage=[30,15,8]),
        expand("storage/assemblies/ecoli/00-shassembly-SRR3191692/timelog{logs}.txt", logs=['_ntcard','_denoise','_contiger','_minia']),
        expand("storage/assemblies/truncatus/00-supernova-chromium/timelog{logs}.txt", logs=['','_output']),
        expand("storage/assemblies/human/00-supernova-T2T_10X_NovaSeq/timelog{logs}.txt", logs=['','_output'])
    output:
        "results/csv/timelog.csv"
    shell:
        """
        echo "folder,cpu_time_s,elapsed_time_s,memory_kb" > {output}
        for f in storage/assemblies/*/*; do
            if [ $(ls -1 ${{f}}/timelog*.txt 2>/dev/null | wc -l) != 0 ]; then
                printf "$f/," >> {output}
                grep ' time (seconds)' ${{f}}/timelog*.txt | awk '{{time += $NF}}END{{printf "%.0f,", time}}' >> {output}
                grep 'Elapsed (wall clock) time (h:mm:ss or m:ss):' ${{f}}/timelog*.txt | awk '{{timestr=$NF;p=index(timestr,":");ctime=0;while(p>0){{ctime=ctime*60+substr(timestr,1,p-1);timestr=substr(timestr,p+1);p=index(timestr,":")}};ctime=ctime*60+timestr;time+=ctime}}END{{printf "%.0f,", time}}' >> {output}
                grep 'Maximum resident set size (kbytes):' ${{f}}/timelog*.txt | awk 'BEGIN{{mem=0}}{{if($NF > mem){{mem = $NF}}}}END{{printf "%d\\n", mem}}' >> {output}
            else
                echo "${{f}} does not contain timelog file"
            fi
        done
        """

rule figure_quast:
    input:
        expand("storage/assemblies/ecoli/00-flye-pacbio_clr_{coverage}x/quast_assembly/transposed_report.tsv", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/00-flye-pacbio_clr_{coverage}x/quast_assembly/transposed_report.tsv", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/00-flye-hifi_{coverage}x/quast_assembly/transposed_report.tsv", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/00-flye-nanopore_{coverage}x/quast_assembly/transposed_report.tsv", coverage=[61,30,15,8]),
        expand("storage/assemblies/ecoli/01-gapless-pacbio_clr_{coverage}x-shassembly-SRR3191692/quast_gapless/transposed_report.tsv", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/01-gapless-pacbio_clr_{coverage}x-supernova-chromium/quast_gapless/transposed_report.tsv", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/01-gapless-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapless/transposed_report.tsv", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/01-gapless-nanopore_{coverage}x-flye-hifi_33x/quast_gapless/transposed_report.tsv", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/human/01-gapless-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapless/transposed_report.tsv", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/ecoli/02-lr_gapcloser-lrscaf-pacbio_clr_{coverage}x-shassembly-SRR3191692/quast_gapclosed/transposed_report.tsv", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/02-lr_gapcloser-lrscaf-pacbio_clr_{coverage}x-supernova-chromium/quast_gapclosed/transposed_report.tsv", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapclosed/transposed_report.tsv", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-nanopore_{coverage}x-flye-hifi_33x/quast_gapclosed/transposed_report.tsv", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/human/02-lr_gapcloser-lrscaf-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapclosed/transposed_report.tsv", coverage=[121,61,30,15,8]),
        expand("storage/assemblies/ecoli/02-tgs_gapcloser-lrscaf-pacbio_clr_{coverage}x-shassembly-SRR3191692/quast_gapclosed/transposed_report.tsv", coverage=[113,57,28,14,7]),
        expand("storage/assemblies/truncatus/02-tgs_gapcloser-lrscaf-pacbio_clr_{coverage}x-supernova-chromium/quast_gapclosed/transposed_report.tsv", coverage=[86,43,21,11,5]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-hifi_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapclosed/transposed_report.tsv", coverage=[33,16,8,4]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-nanopore_{coverage}x-supernova-T2T_10X_NovaSeq/quast_gapclosed/transposed_report.tsv", coverage=[30,15,8]),
        expand("storage/assemblies/human/02-tgs_gapcloser-lrscaf-nanopore_{coverage}x-flye-hifi_33x/quast_gapclosed/transposed_report.tsv", coverage=[8]),
        "storage/assemblies/ecoli/00-shassembly-SRR3191692/quast_minia.contigs/transposed_report.tsv",
        "storage/assemblies/truncatus/00-supernova-chromium/quast_pseudohap.1/transposed_report.tsv",
        "storage/assemblies/human/00-supernova-T2T_10X_NovaSeq/quast_pseudohap.1/transposed_report.tsv",
        list="input/csv/assemblies.csv"
    output:
        "results/figures/figure_quast_NG50_vs_misassemblies.pdf",
        "results/figures/figure_quast_completeness_vs_duplication.pdf"
    shell:
        "Rscript bin/figure_quast.R . {input.list} {output}"

rule figure_coverage_cdf:
    input:
        lambda wildcards: "work/assemblies/human/01-gapless-nanopore_61x-flye-hifi_33x/pass1/gapless_stats.pdf" if "1" == wildcards.i else "work/assemblies/human/01-gapless-nanopore_61x-supernova-T2T_10X_NovaSeq/pass1/gapless_stats.pdf"
    output:
        "results/figures/supplementary/figure_coverage_cdf{i}.pdf"
    shell:
        "pdftk {input} cat 5 output {output}"

rule figure_comp_res_gapless:
    input:
        "storage/assemblies/truncatus/01-gapless-pacbio_clr_86x-supernova-chromium/timelog.pdf",
        "storage/assemblies/human/01-gapless-hifi_33x-supernova-T2T_10X_NovaSeq/timelog.pdf",
        "storage/assemblies/human/01-gapless-nanopore_61x-flye-hifi_33x/timelog.pdf",
        "storage/assemblies/human/01-gapless-nanopore_61x-supernova-T2T_10X_NovaSeq/timelog.pdf"
    output:
        "results/figures/supplementary/figure_comp_res_gapless.pdf"
    params:
        workdir="work/figures/figure_comp_res_gapless"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        bash bin/figure_4panels.sh {params.workdir}/
        pdfcrop {params.workdir}/scfinal.pdf {output} 1>/dev/null
        """

rule figure_2panels:
    input:
        "results/figures/subfigures/figure_graph_traversal.pdf",
        "results/figures/subfigures/figure_graph_traversal_example.pdf"
    output:
        "results/figures/figure_graph_traversal_full.pdf"
    params:
        workdir="work/figures/figure_graph_traversal_full"
    shell:
        """
        mkdir -p {params.workdir}
        i=0
        for fi in {input}; do
            i=$(($i+1))
            pdfcrop $fi {params.workdir}/figure${{i}}.pdf 1>/dev/null
        done
        bash bin/figure_2panels.sh {params.workdir}/
        pdfcrop {params.workdir}/scfinal.pdf {output} 1>/dev/null
        """

rule figure_gapless_comp_res:
    input:
        "storage/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/timelog.csv"
    output:
        "storage/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/timelog.pdf"
    shell:
        "Rscript bin/figure_gapless_comp_res.R {input} {output}"

rule gapless_timelog_csv:
    input:
        expand("work/assemblies/{{species}}/01-gapless-{{data}}_{{cov}}x-{{assembler}}-{{data2}}/pass{it}/timing/{call}.txt", it=[1,2,3], call=["gapless_split","minimap2_repeats","minimap2_reads","gapless_scaffold","minimap2_extension","gapless_extend","gapless_finish","minimap2_consensus","racon"])
    output:
        "storage/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/timelog.csv"
    params:
        folder="work/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}"
    shell:
        """
        echo "pass,call,cpu_time_s,elapsed_time_s,memory_kb" > {output}
        for p in {params.folder}/pass*; do
            for f in $p/timing/*.txt; do
                printf "${{p##*/}}," >> {output}
                ff=${{f##*/}};printf "${{ff%.txt}}," >> {output}
                grep ' time (seconds)' ${{f}} | awk '{{time += $NF}}END{{printf "%.0f,", time}}' >> {output}
                grep 'Elapsed (wall clock) time (h:mm:ss or m:ss):' ${{f}} | awk '{{timestr=$NF;p=index(timestr,":");ctime=0;while(p>0){{ctime=ctime*60+substr(timestr,1,p-1);timestr=substr(timestr,p+1);p=index(timestr,":")}};ctime=ctime*60+timestr;time+=ctime}}END{{printf "%.0f,", time}}' >> {output}
                grep 'Maximum resident set size (kbytes):' ${{f}} | awk 'BEGIN{{mem=0}}{{if($NF > mem){{mem = $NF}}}}END{{printf "%d\\n", mem}}' >> {output}
            done
        done
        """

ref_dict = {}
ref_dict['human'] = "chm13.draft_v1.1.fasta.gz"
ref_dict['truncatus'] = "GCF_011762595.1_mTurTru1.mat.Y_genomic.fna"
ref_dict['ecoli'] = "../../../storage/ecoli/reference/ERS764956.fa"
rule quast:
    input:
        asm="storage/assemblies/{species}/{asmmeth}-{data}/{asmfile}.fasta.gz",
        ref=lambda wildcards: "input/{}/reference/{}".format(wildcards.species, ref_dict[wildcards.species])
    output:
        "storage/assemblies/{species}/{asmmeth}-{data}/quast_{asmfile}/report.pdf",
        "storage/assemblies/{species}/{asmmeth}-{data}/quast_{asmfile}/transposed_report.tsv"
    params:
        outdir="storage/assemblies/{species}/{asmmeth}-{data}/quast_{asmfile}",
        fragmented=lambda wildcards: "" if wildcards.species in ["human","ecoli"] else "--fragmented"
    threads: 2
    log:
        "logs/assemblies/{species}/quast/{asmmeth}-{data}-{asmfile}.txt"
    shell:
        "quast-lg.py -t {threads} -o {params.outdir} -r {input.ref} {params.fragmented} -s -l {wildcards.asmmeth} --no-icarus {input.asm} >{log} 2>&1"

rule samba_finish:
    input:
        "work/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}/input_asm.fa.scaffolds.fa"
    output:
        "storage/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}/scaffolds.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

assembler_output_dict = {}
assembler_output_dict['supernova'] = "pseudohap.1.fasta.gz"
assembler_output_dict['shassembly'] = "minia.contigs.fasta.gz"
assembler_output_dict['flye'] = "assembly.fasta.gz"
rule samba:
    input:
        asm=lambda wildcards: "storage/assemblies/{}/00-{}-{}/{}".format(wildcards.species, wildcards.assembler, wildcards.data2, assembler_output_dict[wildcards.assembler]),
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "work/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}/input_asm.fa.scaffolds.fa",
        timelog="storage/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}/timelog_scaf.txt"
    params:
        datadir="work/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}_data/",
        outdir="work/assemblies/{species}/01-samba-{data}_{cov}x-{assembler}-{data2}/",
        minlen=lambda wildcards: "-m 2500" if wildcards.species == "ecoli" else ""
    log:
        "logs/assemblies/{species}/samba/01-samba-{data}_{cov}x-{assembler}-{data2}_scaf.log"
    threads: max_threads
    shell:
        """
        export org_path=$(pwd)
        mkdir -p {params.datadir}
        zcat {input.reads} > {params.datadir}/reads.fq
        zcat {input.asm} > {params.datadir}/input_asm.fa
        cd {params.outdir}
        env time -v -o ${{org_path}}/{output.timelog} bash ${{org_path}}/bin/masurca/bin/samba.sh -t {threads} {params.minlen} -r ${{org_path}}/{params.datadir}/input_asm.fa -q ${{org_path}}/{params.datadir}/reads.fq >${{org_path}}/{log} 2>&1
        cd -
        rm -f {params.datadir}/{{reads.fa,input_asm.fa}}
        """

rule TGS_GapCloser_finish:
    input:
        "work/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/gapclosed.scaff_seqs"
    output:
        "storage/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/gapclosed.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

tgs_gapcloser_data = {}
tgs_gapcloser_data['pacbio_clr'] = 'pb'
tgs_gapcloser_data['hifi'] = "pb --minmap_arg '-x asm20'"
tgs_gapcloser_data['nanopore'] = 'ont'
rule TGS_GapCloser:
    input:
        asm="storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}/scaffolds.fasta.gz",
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "work/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/gapclosed.scaff_seqs",
        timelog="storage/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/timelog.txt",
    params:
        datadir="work/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}_data/",
        outdir="work/assemblies/{species}/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/",
        type=lambda wildcards: tgs_gapcloser_data[wildcards.data]
    log:
        "logs/assemblies/{species}/tgs_gapcloser/02-tgs_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}.log"
    threads: max_threads
    shell:
        """
        export org_path=$(pwd)
        mkdir -p {params.datadir}
        seqtk seq -A {input.reads} > {params.datadir}/reads.fa
        zcat {input.asm} > {params.datadir}/input_asm.fa
        cd {params.outdir}
        env time -v -o ${{org_path}}/{output.timelog} bash ${{org_path}}/bin/TGS-GapCloser/TGS-GapCloser.sh --thread {threads} --scaff ${{org_path}}/{params.datadir}/input_asm.fa --reads ${{org_path}}/{params.datadir}/reads.fa --output gapclosed --racon ${{org_path}}/bin/racon/racon --tgstype {params.type} >${{org_path}}/{log} 2>&1
        cd -
        rm -f {params.datadir}/{{reads.fa,input_asm.fa}}
        """

rule LR_Gapcloser_finish:
    input:
        lambda wildcards: "work/assemblies/{}/02-lr_gapcloser-lrscaf-{}_{}x-{}-{}/iteration-{}/gapclosed.fasta".format(wildcards.species, wildcards.data, wildcards.cov, wildcards.assembler, wildcards.data2, 2 if wildcards.species == "ecoli" and wildcards.cov == "113" else 3)
    output:
        "storage/assemblies/{species}/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/gapclosed.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

lr_gapcloser_data = {}
lr_gapcloser_data['pacbio_clr'] = 'p'
lr_gapcloser_data['hifi'] = 'p'
lr_gapcloser_data['nanopore'] = 'n'
rule LR_Gapcloser:
    input:
        asm="storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}/scaffolds.fasta.gz",
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "work/assemblies/{species}/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/iteration-3/gapclosed.fasta",
        timelog="storage/assemblies/{species}/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/timelog.txt",
    params:
        datadir="work/assemblies/{species}/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}_data/",
        outdir="work/assemblies/{species}/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}/",
        type=lambda wildcards: lr_gapcloser_data[wildcards.data]
    log:
        "logs/assemblies/{species}/lr_gapcloser/02-lr_gapcloser-lrscaf-{data}_{cov}x-{assembler}-{data2}.log"
    threads: max_threads
    shell:
        """
        export org_path=$(pwd)
        mkdir -p {params.datadir}
        seqtk seq -A {input.reads} > {params.datadir}/reads.fa
        zcat {input.asm} > {params.datadir}/input_asm.fa
        rm -rf {params.outdir}
        cd bin/LR_Gapcloser/src
        env time -v -o ${{org_path}}/{output.timelog} bash LR_Gapcloser.sh -t {threads} -i ${{org_path}}/{params.datadir}/input_asm.fa -l ${{org_path}}/{params.datadir}/reads.fa -s {params.type} -o ${{org_path}}/{params.outdir} >${{org_path}}/{log} 2>&1
        cd -
        rm -f {params.datadir}/{{reads.fa,input_asm.fa}}
        """

minimap_data = {}
minimap_data['pacbio_clr'] = 'map-pb'
minimap_data['hifi'] = 'asm20'
minimap_data['nanopore'] = 'map-ont'
rule lrscaf:
    input:
        asm=lambda wildcards: "storage/assemblies/{}/00-{}-{}/{}".format(wildcards.species, wildcards.assembler, wildcards.data2, assembler_output_dict[wildcards.assembler]),
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}/scaffolds.fasta.gz",
        timelog="storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}/timelog_minimap2.txt",
        timelog2="storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}/timelog_lrscaf.txt"
    params:
        workdir="work/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}",
        outdir="storage/assemblies/{species}/01-lrscaf-{data}_{cov}x-{assembler}-{data2}",
        mapx=lambda wildcards: minimap_data[wildcards.data]
    log:
        minimap2="logs/assemblies/{species}/lrscaf/01-lrscaf-{data}_{cov}x-{assembler}-{data2}_minimap2.log",
        lrscaf="logs/assemblies/{species}/lrscaf/01-lrscaf-{data}_{cov}x-{assembler}-{data2}_lrscaf.log"
    threads: max_threads
    shell:
        """
        mkdir -p {params.workdir}
        env time -v -o {output.timelog} minimap2 -t {threads} -x {params.mapx} {input.asm} {input.reads} >{params.workdir}/aln.mm 2>{log.minimap2}
        env time -v -o {output.timelog2} java -Xms100g -Xmx100g -jar bin/LRScaf-1.1.11.jar -p {threads} -c <(zcat {input.asm}) -a {params.workdir}/aln.mm -t mm -o {params.outdir}/ >{log.lrscaf} 2>&1
        gzip {params.outdir}/scaffolds.fasta
        """

rule gapless_finish:
    input:
        "work/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/gapless.fa"
    output:
        "storage/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/gapless.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

datatype = {}
datatype['pacbio_clr'] = 'pb_clr'
datatype['hifi'] = 'pb_hifi'
datatype['nanopore'] = 'nanopore'
rule gapless:
    input:
        asm=lambda wildcards: "storage/assemblies/{}/00-{}-{}/{}".format(wildcards.species, wildcards.assembler, wildcards.data2, assembler_output_dict[wildcards.assembler]),
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "work/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/gapless.fa",
        timelog="storage/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/timelog.txt",
    params:
        outdir="work/assemblies/{species}/01-gapless-{data}_{cov}x-{assembler}-{data2}/",
        datatype=lambda wildcards: datatype[wildcards.data]
    threads: max_threads
    shell:
        """
        env time -v -o {output.timelog} gapless.sh -j {threads} -o {params.outdir} -i {input.asm} -t {params.datatype} {input.reads}
        """

rule flye_finish:
    input:
        "work/assemblies/{species}/00-flye-{data}/assembly.fasta"
    output:
        "storage/assemblies/{species}/00-flye-{data}/assembly.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

flye_par_dict = {}
flye_par_dict['nanopore'] = "--nano-raw"
flye_par_dict['hifi'] = "--pacbio-hifi"
flye_par_dict['pacbio_clr'] = "--pacbio-raw"
rule flye:
    input:
        reads="work/data/{species}/{data}_{cov}x.fastq.gz"
    output:
        "work/assemblies/{species}/00-flye-{data}_{cov}x/assembly.fasta",
        timelog="storage/assemblies/{species}/00-flye-{data}_{cov}x/timelog.txt"
    params:
        outdir="work/assemblies/{species}/00-flye-{data}_{cov}x/",
        data=lambda wildcards: flye_par_dict[wildcards.data]
    log:
        "logs/assemblies/{species}/flye/00-flye-{data}_{cov}x.txt"
    threads: max_threads
    shell:
        "env time -v -o {output.timelog} flye {params.data} {input} --out-dir {params.outdir} --threads {threads} >{log} 2>&1"

supsample_full_cov_dict = {}
supsample_full_cov_dict[('human','hifi')] = 33
supsample_full_cov_dict[('human','nanopore')] = 121
supsample_full_cov_dict[('truncatus','pacbio_clr')] = 86
supsample_full_cov_dict[('ecoli','pacbio_clr')] = 113
supsample_divisor_dict = {}
supsample_divisor_dict[('human','hifi','16')] = 2
supsample_divisor_dict[('human','hifi','8')] = 4
supsample_divisor_dict[('human','hifi','4')] = 8
supsample_divisor_dict[('human','nanopore','61')] = 2
supsample_divisor_dict[('human','nanopore','30')] = 4
supsample_divisor_dict[('human','nanopore','15')] = 8
supsample_divisor_dict[('human','nanopore','8')] = 16
supsample_divisor_dict[('truncatus','pacbio_clr','43')] = 2
supsample_divisor_dict[('truncatus','pacbio_clr','21')] = 4
supsample_divisor_dict[('truncatus','pacbio_clr','11')] = 8
supsample_divisor_dict[('truncatus','pacbio_clr','5')] = 16
supsample_divisor_dict[('ecoli','pacbio_clr','57')] = 2
supsample_divisor_dict[('ecoli','pacbio_clr','28')] = 4
supsample_divisor_dict[('ecoli','pacbio_clr','14')] = 8
supsample_divisor_dict[('ecoli','pacbio_clr','7')] = 16
rule subsampled_data:
    input:
        lambda wildcards: "work/data/{}/{}_{}x.fastq.gz".format(wildcards.species, wildcards.type, supsample_full_cov_dict[(wildcards.species, wildcards.type)])
    output:
        "work/data/{species}/{type}_{cov}x.fastq.gz"
    params:
        divisor=lambda wildcards: supsample_divisor_dict[(wildcards.species, wildcards.type, wildcards.cov)],
        frac=lambda wildcards: 1/supsample_divisor_dict[(wildcards.species, wildcards.type, wildcards.cov)]
    shell:
        "seqtk sample -s20210610{params.divisor} {input} {params.frac} | gzip > {output}"

rule clr_full_data_ecoli:
    input:
        "storage/ecoli/ERR1036235/real/ERR1036235.fq.gz"
    output:
        "work/data/ecoli/pacbio_clr_113x.fastq.gz"
    shell:
        "ln -s ../../../{input} {output}"

rule clr_full_data_truncatus:
    input:
        glob.glob('input/truncatus/PacBio_CLS/*.subreads.bam')
    output:
        "work/data/truncatus/pacbio_clr_86x.fastq.gz"
    params:
        out="work/data/truncatus/pacbio_clr_86x.fastq"
    shell:
        """
        rm -f {output} {params.out}
        for fq in {input}; do
          samtools fastq $fq >> {params.out}
        done
        gzip {params.out}
        """

rule hifi_full_data_human:
    input:
        "input/human/PRJNA530776_SequelII_HiFi/SRR11292120_subreads.fastq.gz",
        "input/human/PRJNA530776_SequelII_HiFi/SRR11292121_subreads.fastq.gz",
        "input/human/PRJNA530776_SequelII_HiFi/SRR11292122_subreads.fastq.gz",
        "input/human/PRJNA530776_SequelII_HiFi/SRR11292123_subreads.fastq.gz"
    output:
        "work/data/human/hifi_33x.fastq.gz"
    shell:
        "zcat {input} | gzip > {output}"

rule nanopore_full_data_human:
    input:
        "input/human/PRJNA559484_Nanopore/rel7.fastq.gz"
    output:
        "work/data/human/nanopore_121x.fastq.gz"
    shell:
        "ln -s ../../../{input} {output}"

rule gapless_test:
    input:
        asm=lambda wildcards: "storage/{0}/reference/test_sets/{1}/{2}/{2}.fa".format(wildcards.species,wildcards.reffile,wildcards.testset),
        reads=lambda wildcards: "storage/{}/{}/{}/{}.fq.gz".format(wildcards.species, wildcards.sample, wildcards.type, wildcards.sample if "real" == wildcards.type else "PaSS")
    output:
        "storage/{species}/reference/test_sets/{reffile}/{testset}/{sample}/{type}/gapless/pass{iter,[0-9]*}/gapless.fa"
    params:
        outdir="storage/{species}/reference/test_sets/{reffile}/{testset}/{sample}/{type}/gapless"
    threads: max_threads
    shell:
        """
        gapless.sh -j {threads} -n {wildcards.iter} -o {params.outdir} -i {input.asm} -t pb_clr {input.reads}
        """

supernova_maxreads = {}
supernova_maxreads[('truncatus','chromium')] = "1001737466"
supernova_maxreads[('human','T2T_10X_NovaSeq')] = "all"
chromium_data_files = {}
chromium_data_files[('truncatus','chromium')] = ["mTurTru1_S1_L001_R1_001.fastq.gz", "mTurTru1_S1_L001_R2_001.fastq.gz", "mTurTru1_S2_L001_R1_001.fastq.gz", "mTurTru1_S2_L001_R2_001.fastq.gz", "mTurTru1_S3_L001_R1_001.fastq.gz", "mTurTru1_S3_L001_R2_001.fastq.gz", "mTurTru1_S4_L001_R1_001.fastq.gz", "mTurTru1_S4_L001_R2_001.fastq.gz"]
chromium_data_files[('human','T2T_10X_NovaSeq')] = ["CHM13_prep5_S13_L002_R1_001.fastq.gz", "CHM13_prep5_S13_L002_R2_001.fastq.gz", "CHM13_prep5_S14_L002_R1_001.fastq.gz", "CHM13_prep5_S14_L002_R2_001.fastq.gz", "CHM13_prep5_S15_L002_R1_001.fastq.gz", "CHM13_prep5_S15_L002_R2_001.fastq.gz", "CHM13_prep5_S16_L002_R1_001.fastq.gz", "CHM13_prep5_S16_L002_R2_001.fastq.gz"]
rule supernova_run:
    input:
        lambda wildcards: expand("input/{}/{}/{{files}}".format(wildcards.species, wildcards.data), files=chromium_data_files[(wildcards.species, wildcards.data)])
    output:
        "storage/assemblies/{species}/00-supernova-{data}/raw.fasta.gz",
        "storage/assemblies/{species}/00-supernova-{data}/megabubbles.fasta.gz",
        "storage/assemblies/{species}/00-supernova-{data}/pseudohap.1.fasta.gz",
        "storage/assemblies/{species}/00-supernova-{data}/pseudohap.2.fasta.gz",
        timelog="storage/assemblies/{species}/00-supernova-{data}/timelog.txt",
        timelog2="storage/assemblies/{species}/00-supernova-{data}/timelog_output.txt"
    params:
        inputdir="input/{species}/{data}/",
        dir="assemblies/{species}/00-supernova-{data}/",
        maxreads=lambda wildcards: supernova_maxreads[(wildcards.species, wildcards.data)]
    log:
        run="logs/assemblies/{species}/supernova/{data}_run.log",
        mkoutput="logs/assemblies/{species}/supernova/{data}_mkoutput.log"
    threads: max_threads
    shell:
        """
        export org_path=$(pwd)
        mkdir -p work/{params.dir}
        cd work/{params.dir}
        env time -v -o ${{org_path}}/{output.timelog} supernova run --localcores={threads} --id=assembly --fastqs=${{org_path}}/{params.inputdir} --description="Truncatus assembly full coverage" --maxreads={params.maxreads} 1>${{org_path}}/{log.run} 2>&1 &&\
        supernova mkoutput --asmdir=assembly/outs/assembly --outprefix=raw --style=raw 1>${{org_path}}/{log.mkoutput} 2>&1 &&\
        supernova mkoutput --asmdir=assembly/outs/assembly --outprefix=megabubbles --style=megabubbles --minsize=0 1>>${{org_path}}/{log.mkoutput} 2>&1 &&\
        env time -v -o ${{org_path}}/{output.timelog2} supernova mkoutput --asmdir=assembly/outs/assembly --outprefix=pseudohap --style=pseudohap2 --minsize=0 1>>${{org_path}}/{log.mkoutput} 2>&1 &&\
        cd - &&\
        mv work/{params.dir}/*.fasta.gz work/{params.dir}/assembly/_log work/{params.dir}/assembly/outs/report.txt work/{params.dir}/assembly/outs/summary.csv work/{params.dir}/assembly/outs/assembly/stats storage/{params.dir} &&\
        rm -rf work/{params.dir}
        """

rule sh_assembly_finish:
    input:
        "work/assemblies/{species}/00-shassembly-{data}/minia.contigs.fa"
    output:
        "storage/assemblies/{species}/00-shassembly-{data}/minia.contigs.fasta.gz"
    shell:
        "cat {input} | gzip > {output}"

illumina_data_dict = {}
illumina_data_dict[('ecoli','SRR3191692','fq1')] = ["SRR3191692_1.fq.gz"]
illumina_data_dict[('ecoli','SRR3191692','fq2')] = ["SRR3191692_2.fq.gz"]
illumina_data_dict[('human','PRJNA559484_NovaSeq','fq1')] = ["SRR10035390_1.fastq.gz","SRR10035391_1.fastq.gz","SRR10035392_1.fastq.gz","SRR10035393_1.fastq.gz"]
illumina_data_dict[('human','PRJNA559484_NovaSeq','fq2')] = ["SRR10035390_2.fastq.gz","SRR10035391_2.fastq.gz","SRR10035392_2.fastq.gz","SRR10035393_2.fastq.gz"]
rule sh_assembly:
    input:
        fq1=lambda wildcards: expand("input/{}/{}/{{data}}".format(wildcards.species, wildcards.data), data=illumina_data_dict[(wildcards.species, wildcards.data, 'fq1')]),
        fq2=lambda wildcards: expand("input/{}/{}/{{data}}".format(wildcards.species, wildcards.data), data=illumina_data_dict[(wildcards.species, wildcards.data, 'fq2')])
    output:
        "work/assemblies/{species}/00-shassembly-{data}/minia.contigs.fa",
        timelog="storage/assemblies/{species}/00-shassembly-{data}/timelog_ntcard.txt",
        timelog2="storage/assemblies/{species}/00-shassembly-{data}/timelog_denoise.txt",
        timelog3="storage/assemblies/{species}/00-shassembly-{data}/timelog_contiger.txt",
        timelog4="storage/assemblies/{species}/00-shassembly-{data}/timelog_minia.txt"
    params:
        outdir="work/assemblies/{species}/00-shassembly-{data}"
    log:
        denoise="logs/assemblies/{species}/shassembly/{data}/denoise.log",
        contiger="logs/assemblies/{species}/shassembly/{data}/contiger.log",
        minia="logs/assemblies/{species}/shassembly/{data}/minia.log"
    threads: max_threads
    shell:
        """
        zcat {input.fq1} > {params.outdir}/reads1.fq
        zcat {input.fq2} > {params.outdir}/reads2.fq
        export org_path=$(pwd)
        cd {params.outdir}
        env time -v -o ${{org_path}}/{output.timelog} ntcard -t{threads} -k47 -o freq_k47.hist reads1.fq reads2.fq 1>freq_k47.stats 2>&1
        echo "reads1.fq" > input_list.txt
        echo "reads2.fq" >> input_list.txt
        env time -v -o ${{org_path}}/{output.timelog2} CQF-deNoise -t {threads} -i input_list.txt -o k47.cqf -k 47 $(awk '($2=="F1"){{F1 = $3}}($2=="F0"){{F0 = $3}}(($1 == 47) && ($2 == 1)){{f1 = $3}}(($1 == 47) && ($2 == 2)){{f2 = $3}}END{{print "-N", F1, "-n", F0-f1-f2, "-e", 1-((F1-f1-2*f2)/F1)^(1/28) }}' freq_k47.stats freq_k47.hist) -f g 1>${{org_path}}/{log.denoise} 2>&1
        sleep 30
        env time -v -o ${{org_path}}/{output.timelog3} Contiger -t {threads} -k 47 -i input_list.txt -c k47.cqf -o unitigs.fa -f g 1>${{org_path}}/{log.contiger} 2>&1
        sleep 30
        env time -v -o ${{org_path}}/{output.timelog4} minia -nb-cores {threads} -kmer-size 47 -unitig -in unitigs.fa -out minia 1>${{org_path}}/{log.minia} 2>&1
        rm reads{{1,2}}.fq
        cd -
        """

#---Tests-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
rule create_testsets:
    input:
        "storage/{species}/reference/{reffile}.fa"
    output:
        "storage/{species}/reference/test_sets/{reffile}/{testset,gaps|inversions}/{testset}.fa"
    shell:
        "./bin/create_testset.py {wildcards.testset} {input} {output}"
        

rule create_test_unchanged:
    input:
        "storage/{species}/reference/{reffile}.fa"
    output:
        "storage/{species}/reference/test_sets/{reffile}/unchanged/unchanged.fa"
    shell:
        "ln -srf {input} {output}"

rule PaSS_simulation:
    input:
        index="storage/{species}/{sample}/PaSS/index",
        per="storage/{species}/{sample}/PaSS/percentage.txt",
        config="storage/{species}/{sample}/PaSS/sim.config",
        reads="storage/{species}/{sample}/real/{sample}.fq.gz"
    output:
        "storage/{species}/{sample}/PaSS/PaSS.fq.gz"
    params:
        prefix="storage/{species}/{sample}/PaSS/PaSS"
    threads: max_threads
    shell:
        "PaSS -list {input.per} -index {input.index} -m pacbio_RS -c {input.config} -r $(($(seqtk seq {input.reads} | wc -l)/4)) -t {threads} -o {params.prefix} && gzip {params.prefix}.fq"
        
rule PaSS_index:
    input:
        lambda wildcards: "storage/{}/reference/{}.fa".format(wildcards.species, SAMPLES[wildcards.sample]['reference'])
    output:
        "storage/{species}/{sample}/PaSS/index",
        "storage/{species}/{sample}/PaSS/percentage.txt"
    params:
        outdir="storage/{species}/{sample}/PaSS/"
    shell:
        "bin/PaSS/pacbio_mkindex.pl {input} {params.outdir}"

rule PaSS_profile:
    input:
        reads="storage/{species}/{sample}/real/{sample}.fq.gz",
        map="storage/{species}/{sample}/real/{sample}.blasr.gz"
    output:
        "storage/{species}/{sample}/PaSS/sim.config"
    params:
        workdir="work/{species}/{sample}/PaSS/"
    shell:
        """
        mkdir -p {params.workdir}
        cp bin/PaSS/run.pl {params.workdir}
        cp -r bin/PaSS/script {params.workdir}
        zcat {input.reads} > {params.workdir}/{wildcards.sample}.fq
        zcat {input.map} > {params.workdir}/{wildcards.sample}.blasr
        cd {params.workdir}
        ./run.pl {wildcards.sample}.fq {wildcards.sample}.blasr RS
        cd -
        mv {params.workdir}/sim.config {output}
        rm -rf {params.workdir}
        """

rule blasr:
    input:
        reads="storage/{species}/{sample}/real/{sample}.fq.gz",
        ref=lambda wildcards: "storage/{}/reference/{}.fa".format(wildcards.species, SAMPLES[wildcards.sample]['reference'])
    output:
        "storage/{species}/{sample}/real/{sample}.blasr.gz"
    params:
        out="storage/{species}/{sample}/real/{sample}.blasr",
        tmpdir="work/{species}/{sample}/real/"
    shell:
# Must be run with only one thread, otherwise split file and recombine afterwards, so that alignments keep original order
        """
        mkdir -p {params.tmpdir}
        zcat {input.reads} > {params.tmpdir}/{wildcards.sample}.fq
        blasr --nproc 1 {params.tmpdir}/{wildcards.sample}.fq {input.ref} --allowAdjacentIndels --hitPolicy randombest --out {params.out} -m 0
        rm -f {params.tmpdir}/{wildcards.sample}.fq
        gzip {params.out}
        """

rule pbh5tools:
    input:
        lambda wildcards: glob.glob('input/{}/{}/*.bas.h5'.format(wildcards.species, wildcards.sample))
    output:
        "storage/{species}/{sample}/real/{sample}.fq.gz"
    params:
        prefix="storage/{species}/{sample}/real/{sample}"
    shell:
        """
        bash5tools.py --outFilePrefix {params.prefix} --readType subreads --outType fastq {input}
        mv {params.prefix}.fastq {params.prefix}.fq
        gzip {params.prefix}.fq
        """
        
rule extract_fa_from_gff:
    input:
        "input/{species}/reference/{reference}.gff"
    output:
        "storage/{species}/reference/{reference,[!/]*}.fa"
    shell:
        """
        awk '{{if(store){{print $0}}else{{if($0 == "##FASTA"){{store = 1}}}}}}' {input} > {output}
        """
