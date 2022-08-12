# DASTools and Eggnog mapper both use diamond so I may be able to stack the analyses. I don't know where this is but I should be able to find it. 
# try to figure out a way to organize the cat files. 

# <!> NOTE some of the issues may be trying to re-use the diamond files above (as they may use a different version of diamond)

# "../" -> goes straight to scripts, CAT environment does not currently add CAT_pack to scripts. 

# think about splitting things up to minimize re-work. 

# <!> add in CAT Bin taxonomy

DBDIR = config["database_dir"]
database = config["abricatedb"]

rule CAT_contigs:
    input: 
        contigs = "{sample}/{sample}_contigs.fasta",
    output:
        out = "{sample}/CAT/{sample}.contig2classification.txt",
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
       "../envs/CAT.yaml"
    log:
        "log/CAT/{sample}.contigs.log", 
    benchmark:
        "log/benchmarks/CAT/{sample}.contigs.tsv"
    shell:
       "pwd -P"
       #" mkdir {wildcards.sample}/CAT/ 2> {log} " # there is a problem here....
       " ; "
        " CAT_pack/CAT contigs -c {input.contigs} --force -n {threads} -d {DBDIR}/CAT_database -t {DBDIR}/CAT_taxonomy --out_prefix {wildcards.sample}/CAT/{wildcards.sample}  2>> {log}" # --tmpdir /dev/shm --out_prefix {wildcards.sample}/CAT/{wildcards.sample} "
        " ; "

rule CAT_summarize:
    input: 
        contigs = "{sample}/{sample}_contigs.fasta",
        info = "{sample}/CAT/{sample}.contig2classification.official_names.txt"
    output:
        out = "genomes/taxonomy/CAT/{sample}.summary.txt"
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/CAT.yaml"
    log:
        "log/CAT/{sample}.contigs_names.log", 
    benchmark:
        "log/benchmarks/CAT/{sample}.contigs.tsv"
    shell:
        " CAT_pack/CAT summarise -c {input.contigs} -i {input.info} -o genomes/taxonomy/CAT/{wildcards.sample}.summary.txt --force  2> {log} "
        " ; "
        " cat  genomes/taxonomy/CAT/{wildcards.sample}.summary.txt >> genomes/taxonomy/CAT/CAT.summary.txt  2>> {log} "

rule CAT_add_names_official:
    input: 
        contig2class = "{sample}/CAT/{sample}.contig2classification.txt",
    output:
        out = "{sample}/CAT/{sample}.contig2classification.official_names.txt"
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/CAT.yaml"
    log:
        "log/CAT/{sample}.contigs_names.log", 
    benchmark:
        "log/benchmarks/CAT/{sample}.contigs.tsv"
    shell:
        " CAT_pack/CAT add_names -i {input.contig2class} -o {wildcards.sample}/CAT/{wildcards.sample}.contig2classification.official_names.txt -t {DBDIR}/CAT_taxonomy --only_official --force  2> {log} "

def get_all_CAT(wildcards): 
    
    all_genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome
     
    return expand(rules.CAT_summarize.output.out, genome=all_genomes, sample = SAMPLES)


# the below expand rules are incomplete as they do not cover all output files this should change when database wildcards added. 
rule all_CAT:
    input:
        get_all_CAT,
    output:
        touch(f"/genomes/taxononomy/CAT/finished_CAT"),
    conda:
        "../envs/CAT.yaml"
    shell:
        " mkdir -p genomes/annotations/abricate/temp "
        " ; "
        " mkdir -p genomes/annotations/abricate/taxonomy "


####### BAT ##############

# Binner location hard coded

BINNER = "DASTool"


def get_all_CATabra(wildcards): 
    
    all_genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome
     
    return expand(rules.join.output.out, genome=all_genomes, sample = SAMPLES, database = database)

def get_all_BAT(wildcards): 
    
    all_genomes = glob_wildcards(os.path.join(genome_dir, "{genome}.fasta")).genome
     
    return expand(rules.BAT_summarize.output.out, genome=all_genomes, sample = SAMPLES)


# the below expand rules are incomplete as they do not cover all output files this should change when database wildcards added. 

        
rule BAT_bins:
    input: 
        bins = "{sample}/binning/DASTool/bins",
    output:
        out = "{sample}/BAT/{sample}.bin2classification.txt",
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
       "../envs/CAT.yaml"
    log:
        "log/BAT/{sample}.bins.log", 
    benchmark:
        "log/benchmarks/BAT/{sample}.bins.tsv"
    shell:
    #   " mkdir {wildcards.sample}/BAT/ 2> {log} " # there is a problem here....
    #  " ; "
        " CAT_pack/CAT bins -f 0.5 -b {input.bins} --force -n {threads} -d {DBDIR}/CAT_database -t {DBDIR}/CAT_taxonomy --out_prefix {wildcards.sample}/BAT/{wildcards.sample} -s .fasta  2> {log}" # --tmpdir /dev/shm --out_prefix {wildcards.sample}/BAT/{wildcards.sample} "
        " ; "

rule BAT_summarize:
    input: 
        bins = "{sample}/binning/DASTool/bins",
        info = "{sample}/BAT/{sample}.bin2classification.official_names.txt"
    output:
        out = "genomes/taxonomy/BAT/{sample}.summary.txt"
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/CAT.yaml"
    log:
        "log/BAT/{sample}.bins_names.log", 
    benchmark:
        "log/benchmarks/BAT/{sample}.bins.tsv"
    shell:
        " CAT_pack/CAT summarise -i {input.info} -o genomes/taxonomy/BAT/{wildcards.sample}.summary.txt --force  2> {log} "
        " ; "
        " cat  genomes/taxonomy/BAT/{wildcards.sample}.summary.txt >> genomes/taxonomy/BAT/BAT.summary.txt  2>> {log} "

rule BAT_add_names_official:
    input: 
        bin2class = "{sample}/BAT/{sample}.bin2classification.txt",
    output:
        out = "{sample}/BAT/{sample}.bin2classification.official_names.txt"
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/CAT.yaml"
    log:
        "log/BAT/{sample}.bins_names.log", 
    benchmark:
        "log/benchmarks/BAT/{sample}.bins.tsv"
    shell:
        " CAT_pack/CAT add_names -i {input.bin2class} -o {wildcards.sample}/BAT/{wildcards.sample}.bin2classification.official_names.txt -t {DBDIR}/CAT_taxonomy --only_official --force  2> {log} "

###############   Joining Abricate and CAT/BAT and MAGs  ###################################

# rule processing: # can really be combined into join but here for clarity.
#     input:
#         BATcomplete = "/genomes/taxononomy/CAT/finished_CAT",
#         CATcomplete = "/genomes/taxononomy/CAT/finished_BAT",
#     output:
#         touch(f"genomes/annotations/abricate/temp/cleaned"), 
#     conda:
#         "../envs/CAT.yaml"
#     shell:
#         " sed 's/bins\//\t/g' genomes/annotations/abricate/summary_bin_megares.tab > genomes/annotations/abricate/temp/processed_summary_bin_megares.tab "
#         " ; "
#         " sed 's/genomes\/genomes\///g' genomes/annotations/abricate/summary_mag_megares.tab | sed 's/.fasta//g' > genomes/annotations/abricate/temp/processed_summary_mag_megares.tab "
    

rule processing: # can really be combined into join but here for clarity.
    input:
        BATcomplete = "/genomes/taxononomy/CAT/finished_CAT",
        CATcomplete = "/genomes/taxononomy/CAT/finished_BAT",
    output:
        touch(f"genomes/annotations/abricate/temp/cleaned"), 
    conda:
        "../envs/CAT.yaml"
    shell:
        " sed 's/bins\//\t/g' genomes/annotations/abricate/summary_bin_{wildcards.database}.tab > genomes/annotations/abricate/temp/processed_summary_bin_{wildcards.database}.tab "
        " ; "
        " sed 's/genomes\/genomes\///g' genomes/annotations/abricate/summary_mag_{wildcards.database}.tab | sed 's/.fasta//g' > genomes/annotations/abricate/temp/processed_summary_mag_{wildcards.database}.tab "




# This may require some other control flow.....
rule join: # might not be sorted on the correct row? See uncommitted sort. 
    input: 
        genomesComplete = "finished_genomes",
        abricatecomplete = "finished_abricate",
        BATcomplete = "/genomes/taxononomy/CAT/finished_CAT",
        CATcomplete = "/genomes/taxononomy/CAT/finished_BAT",
        MAGtaxonomy = "genomes/taxonomy/gtdb_taxonomy.tsv",
        cleaned = "genomes/annotations/abricate/temp/cleaned",
        # abricateMAG = "genomes/annotations/abricate/summary_mag_{database}.tab",
        # abricateContig = "genomes/annotations/abricate/summary_contig_{database}.tab",
        # abricateBin = "genomes/annotations/abricate/summary_bin_{database}.tab",
        contig2class = "{sample}/CAT/{sample}.contig2classification.official_names.txt",
        bin2class = "{sample}/BAT/{sample}.bin2classification.official_names.txt",
    output:
        out = "{sample}/CAT/temp/CATabraMerged.{database}",
    threads: config["threads"]
    resources:
        mem=config["simplejob_mem"],
        time=config["runtime"]["default"],
    conda:
        "../envs/CAT.yaml"
    log:
        "log/BAT/{sample}.{database}.join.log", 
    benchmark:
        "log/benchmarks/BAT/{sample}.{database}.join.tsv"
    shell:
        " join  -1 2 -2 1 -t $'\t' -e NA <(sort -f -k2 genomes/annotations/abricate/summary_contig_{wildcards.database}.tab) <(sort -f -k1 {input.contig2class}) >> genomes/annotations/abricate/taxonomy/contigs_{wildcards.database}_tax.txt 2>> {log} "
        " ; "
        " join  -1 2 -2 1 -t $'\t' -e NA  <(sort -f -k2 genomes/annotations/abricate/temp/processed_summary_bin_{wildcards.database}.tab) <(sort -f -k1 {input.bin2class}) >> genomes/annotations/abricate/taxonomy/bin_{wildcards.database}_tax.txt 2>> {log} "
        " ; "
        " join  -1 1 -2 1 -t $'\t' -e NA <(sort -f -k1 genomes/annotations/abricate/temp/processed_summary_mag_{wildcards.database}.tab) <(sort -f -k1 {input.MAGtaxonomy}) >> genomes/annotations/abricate/taxonomy/mag_{wildcards.database}_tax.txt 2>> {log} "
        " ; "
        " touch {output.out} "


rule all_BAT:
    input:
        get_all_BAT,
    output:
        touch(f"/genomes/taxononomy/CAT/finished_BAT"),
    conda:
        "../envs/CAT.yaml"
    shell:
        "CAT_pack/CAT --help"


rule all_CATabra:
    input:
        get_all_CATabra,
    output:
        touch(f"/genomes/taxononomy/CAT/finished_CATabra"),
    conda:
        "../envs/CAT.yaml"
    shell:
        "CAT_pack/CAT --help"