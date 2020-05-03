
#you do not need this if you supply your email address directly to rule "download_seqs"
configfile: "snakemake_config.json"

#recommend to run the make_database commands first, so you can ensure this all goes well, before trying to do
#the nextstrain run steps. Just for sanity.

rule all:
    input:
        ["auspice/CoV_229E.json", #"auspice/cov.json",
         "auspice/CoV_NL63.json",
         "auspice/CoV_OC43.json",
         "auspice/CoV_HKU1.json"]

rule make_database_all:
    input:
        expand("genbank/CoV-{id}_genbank_meta.tsv", id=["229E", "NL63", "OC43","HKU1"])

#calls for the augur runs
rule build_HKU1:
    input: "auspice/CoV_HKU1_spike.json"
rule build_229E:
    input: "auspice/CoV_229E_spike.json"
rule build_NL63:
    input: "auspice/CoV_NL63_spike.json"
rule build_OC43:
    input: "auspice/CoV_OC43_spike.json"

#calls for individual database runs
rule make_database_HKU1:
    input:
        "genbank/CoV-HKU1_genbank_meta.tsv"
rule make_database_229E:
    input:
        "genbank/CoV-229E_genbank_meta.tsv"
rule make_database_NL63:
    input:
        "genbank/CoV-NL63_genbank_meta.tsv"
rule make_database_OC43:
    input:
        "genbank/CoV-OC43_genbank_meta.tsv"

wildcard_constraints:
    gene="|-last7|-replicase",
    id="229E|NL63|OC43|HKU1"

rule files:
    input:
        #raw vipr tab-delimited download files go here. Change the file name to match your latest download version
        file_229E = "data/vipr_download_229E_01May2020.tsv", #"data/CoV-229E-23Jan20.tsv",
        file_NL63 = "data/vipr_download_NL63_01May2020.tsv", #"data/CoV-NL63-23Jan20.tsv",
        file_beta = "data/vipr_download_OC43_01May2020.tsv", #"data/CoV-Betacoronavirus1-23Jan20.tsv",
        file_HKU1 = "data/vipr_download_HKU1_01May2020.tsv", #"data/CoV-HKU1-19Apr20.tsv",

        #currently manual samples are not used.
        #samples added manually - not from ViPR (ensure not duplicates of anything downloaded)
        #manual_seqs = "data/sequences.fasta",
        #manual_meta = "data/metadata.tsv",

        #config
        regions = "scripts/geo_regions.tsv",
        #have separate align and translate references
        #so that modifying the translate (gene name etc) doesn't
        #cause everything to re-align!!! (takes ages for some viruses!!)
        align_reference = "{id}/config/CoV-{id}-align_reference.gb", #<- try not to modify this unless want re-align!
        translate_reference = "{id}/config/CoV-{id}-translate_reference.gb", #<- modify this for gene names etc!
        
        extra_meta = "{id}/data/extra_meta_{id}.tsv",
        include = "{id}/config/include_{id}.txt",
        exclude = "{id}/config/exclude_{id}.txt", 

        description = "{id}/config/description_{id}.md",
        auspice_config = "config/auspice_config.json",

        colors = "config/colors.tsv",
        lat_longs = "config/lat_longs.tsv",

        blast_ref = "{id}/config/CoV-{id}-translate_reference.fasta",


files = rules.files.input

def is_rerun(wildcards):
    #these lines are just debug...
    #print("genbank/CoV-{}_current_download.tsv".format(wildcards))
    #print(os.path.isfile("genbank/CoV-{}_current_download.tsv".format(wildcards)))
    return os.path.isfile("genbank/CoV-{}_current_download.tsv".format(wildcards))

VIPR_FILES = {"229E": files.file_229E,
              "NL63": files.file_NL63,
              "OC43": files.file_beta,
              "HKU1": files.file_HKU1}

##############################
# Parse metadata from ViPR
# some strain names are identical. So adds accession number like this:
#   strain__ACCESSION
# This way all names are unique. Original strain name is stored in new column
# 'orig_strain'
###############################

rule parse_vipr_meta:
    input:
        #meta = lambda wildcards: "data/CoV-{}.tsv".format(vipr_files[wildcards.id]),
        meta = lambda wildcards: VIPR_FILES[wildcards.id],
        regions = ancient(files.regions)
    output:
        out = "temp/CoV-{id}_current_download.tsv"
    params:
        rerun = lambda wildcards: is_rerun(wildcards.id),
        leng = "vp1" #sorry, but this is needed to trigger whether to add the accession
    shell:
        """
        # Figure out message to show user.
        rrun={params.rerun}
        if [ $rrun == "True" ]; then
            echo "This {wildcards.id} rerun will use existing GenBank files! Only new accession numbers will be downloaded"
        else
            echo "Starting new {wildcards.id} run from scratch. All VIPR samples will be downloaded."
        fi

        python scripts/vipr_parse.py --input {input.meta} --output {output.out} \
            --regions {input.regions} \
            --length {params.leng}
        """

##############################
# FIND NEW
# find only new seqs to download - can be adapated to exclude manual sequences.
# Send possible link to current download file - if empty (newrun) will be ignored.
# If not empty (rerun), these will also be ignored
###############################
rule find_new:
    input:
        new_meta = rules.parse_vipr_meta.output.out
    params:
        old_meta = ancient("genbank/CoV-{id}_current_download.tsv")
    output:
        "temp/CoV-{id}_meta_to_download.tsv" 
    shell:
        """
        python scripts/find_new.py --input-new {input.new_meta} \
            --exclude {params.old_meta} \
            --output {output}
        """
        #{input.man_meta}


##############################
# Download from Genbank only new and non-duplicate sequences
###############################
rule download_seqs:
    input:
        meta = rules.find_new.output[0] 
    output:
        sequences = "temp/CoV-{id}_download_seqs.fasta", 
        meta = "temp/CoV-{id}_download_meta.tsv" 
    run:
        import pandas as pd
        from Bio import Entrez, SeqIO
        from augur.parse import forbidden_characters
        from datetime import datetime
        Entrez.email = config['email'] #REPLACE WITH YOUR OWN EMAIL ADDRESS 

        print(input.meta)
        meta = pd.read_csv(input.meta, sep='\t')
        originalMetaLen = len(meta)
        additional_meta = {}
        #len_cutoff = 6400 if wildcards.length=="genome" else 300
        len_cutoff = 800
        print("Downloading only sequences with length >= {}".format(len_cutoff))

        tooShort = []
        didntWork = []
        with open(output.sequences, 'w') as fh:
            for ri, row in meta.iterrows():
                try:
                    handle = Entrez.efetch(db="nucleotide", id=row.accession, rettype="gb", retmode="text")
                except:
                    print(row.accession, "did not work")
                    didntWork.append("{}\t{}".format(row.strain, row.accession))
                    meta.drop(ri, inplace=True)
                    continue
                print(row.strain, row.accession)
                rec = SeqIO.read(handle, 'genbank')
                if len(rec.seq) - rec.seq.count("N") < len_cutoff:
                    print(row.strain, row.accession, "is too short when Ns removed!")
                    tooShort.append("{}\t{}".format(row.strain, row.accession))
                    meta.drop(ri, inplace=True)
                    continue
                try:
                    authors = rec.annotations['references'][0].authors
                    title = rec.annotations['references'][0].title
                except:
                    authors = ''
                    title = ''

                url = 'https://www.ncbi.nlm.nih.gov/nuccore/'+row.accession
                add_date = datetime.today().strftime('%Y-%m-%d')
                additional_meta[ri] = {'url':url, 'authors':authors, 'title':title, 'date_added':add_date}
                tmp = row.strain
                for c,r in forbidden_characters:
                    tmp=tmp.replace(c,r)
                rec.id = tmp
                rec.name = tmp
                rec.description = ''
                SeqIO.write(rec, fh, 'fasta')

        print("\n")
        print(len(tooShort), "sequences were too short after Ns were removed, and were excluded.")
        shortFile = "temp/too_short_{}.txt".format(wildcards.id)
        if tooShort:
            with open(shortFile, 'w') as f:
                for item in tooShort:
                    f.write("%s\n" % item)
            print("You can see those excluded as too short in '{}'".format(shortFile))

        print(len(didntWork), "sequences weren't able to be downloaded and were excluded.")
        didntFile = "temp/didnt_work_{}.txt".format(wildcards.id)
        if didntWork:
            with open(didntFile, 'w') as f:
                for item in didntWork:
                    f.write("%s\n" % item)
            print("You can see those that failed to download in '{}'".format(didntFile))

        print("\nOf {} files we tried to download, {} were downloaded.".format(originalMetaLen,len(meta)))
        add_meta = pd.DataFrame(additional_meta).transpose()
        all_meta = pd.concat((meta, add_meta), axis=1)
        all_meta.to_csv(output.meta, sep='\t', index=False)


#####################################################################################################
#    BLAST - or not
#       If VP1, need to BLAST. if genome, do not, but must rename files.
#####################################################################################################

##############################
# If VP1 - BLAST
###############################
rule blast:        
    input:
        blast_db_file = files.blast_ref,
        #seqs_to_blast =  "vp1/temp/downloaded_seqs.fasta" #from rule download.seqs output
        seqs_to_blast = "temp/CoV-{id}_download_seqs.fasta"
    output:
        blast_out = "temp/blast_out_{id}.csv"
    run:
        from subprocess import call
        import os
        comm = ["makeblastdb -in", input.blast_db_file, "-out vp1/temp/entero_db_{} -dbtype nucl".format(wildcards.id)]
        cmd = " ".join(comm)
        os.system(cmd)
        comm2 = ["blastn -task blastn -query", input.seqs_to_blast, "-db vp1/temp/entero_db_{} -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -out".format(wildcards.id), output.blast_out, "-evalue 0.0005"]
        cmd2 = " ".join(comm2)
        os.system(cmd2)

        #blast output is specified by the call - here is what the columns are:
        #(query is the Genbank seq we are looking for a match in, subject is the reference VP1 alignment)
        ##First 7 columns:
        #10     qseqid  sseqid  pident  length      mismatch    gapopen 
        #csv    queryID subjID  %match  alignLen    #mismatch   #gap    

        ##Next 4 columns
        #qstart                     qend                    sstart                  send 
        #start of align in query    end of align in query   start of align in subj  end of align in subj

        ##Last 3 columns
        #evalue         bitscore    qcovs
        #Expect value   Bit score   Querty coverage per subject

##############################
# If VP1 - after BLAST, take only those that match
###############################
rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out,
        #input_meta = "vp1/temp/downloaded_meta.tsv", #from rule download.seqs output
        #input_seqs = "vp1/temp/downloaded_seqs.fasta" #from rule download.seqs output
        input_seqs = "temp/CoV-{id}_download_seqs.fasta", 
        input_meta = "temp/CoV-{id}_download_meta.tsv" 
    output:
        out_seqs = "temp/CoV-{id}_blast_sequences.fasta",
        out_meta = "temp/CoV-{id}_blast_metadata.tsv"
    params:
        matchLen = 300,
        dups = "vp1/temp/duplicates.tsv"
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --meta {input.input_meta} \
            --seqs {input.input_seqs} \
            --out_seqs {output.out_seqs} \
            --out_meta {output.out_meta} \
            --match_length {params.matchLen} \
            --dup_file {params.dups}
        """


################################
# Align whatever's new to reduce time spent aligning whole download again...
################################

rule align_download:
    message:
        """
        Aligning newly downloaded sequences to {input.ref}
          - filling gaps with N
        """
    input:
        sequences = rules.blast_sort.output.out_seqs, #rules.download_seqs.output.sequences,
        ref = files.align_reference
    output:
        out = "temp/CoV-{id}_download_aligned_seqs.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.ref} \
            --output {output.out} \
            --fill-gaps \
            --remove-reference

        # remove debug guff
        rm "{output.out}.log"
        #rm "{output.out}.post_aligner.fasta"
        #rm "{output.out}.pre_aligner.fasta"
        #rm "{input.sequences}.ref.fasta"
        """


#####################################################################################################
#    Bring together old and new files & create database
#####################################################################################################

##############################
# Concat meta and sequences to existing Genbank
# If it's a first run it concats to nothing....
###############################
rule add_meta:
    input:
        metadata = rules.blast_sort.output.out_meta, #rules.download_seqs.output.meta
    output:
        metadata = "temp/CoV-{id}_genbank_meta.tsv"
    params:
        rerun = lambda wildcards: is_rerun(wildcards.id)
    run:
        import pandas as pd
        from augur.parse import fix_dates, forbidden_characters
        md = []
        input_files = [input.metadata]
        if params.rerun:
            input_files.append("genbank/CoV-{}_genbank_meta.tsv".format(wildcards.id))
        for fname in input_files: #input.metadata:
            tmp = pd.read_csv(fname, sep='\t' if fname.endswith('tsv') else ',')
            tmp_name = []
            for x in tmp.strain:
                f = x
                for c,r in forbidden_characters:
                    f=f.replace(c,r)
                tmp_name.append(f)
            tmp.strain = tmp_name
            md.append(tmp)
        all_meta = pd.concat(md)
        all_meta.to_csv(output.metadata, sep='\t', index=False)
#concat sequences - if first run , just copies with different name...
rule add_sequences:
    input:
        "temp/CoV-{id}_download_aligned_seqs.fasta"
    output:
        "temp/CoV-{id}_genbank_seqs.fasta"
    params:
        rerun = lambda wildcards: is_rerun(wildcards.id)
    shell:
        '''
        rrun={params.rerun}
        if [ $rrun == "True" ]; then
            cat {input} "genbank/CoV-{wildcards.id}_genbank_seqs.fasta" > {output}
        else
            cp {input} {output}
        fi
        '''

##############################
# If all has gone well - make a new Database!
###############################
rule make_database:
    input:
        gen_seqs = rules.add_sequences.output, #"temp/CoV-{id}_genbank_seqs.fasta",
        gen_meta = rules.add_meta.output.metadata, #"temp/CoV-{id}_genbank_meta.tsv",
        download = rules.parse_vipr_meta.output.out #"temp/CoV-{id}_current_download.tsv"
    output:
        gen_seqs = "genbank/CoV-{id}_genbank_seqs.fasta",
        gen_meta = "genbank/CoV-{id}_genbank_meta.tsv",
        download = "genbank/CoV-{id}_current_download.tsv"
    params:
        rerun = lambda wildcards: is_rerun(wildcards.id)
    shell:
        '''
        cp {input.gen_seqs} genbank
        cp {input.gen_meta} genbank
        cp {input.download} genbank

        #important to rename stuff to prevent messing up of future snakemake runs bc 'files already there'
        mv {input.gen_meta} "temp/CoV-{wildcards.id}_genbank_meta_old.tsv"
        mv {input.gen_seqs} "temp/CoV-{wildcards.id}_genbank_seqs_old.fasta"
        
        mv temp/CoV-{wildcards.id}_download_meta.tsv  temp/CoV-{wildcards.id}_download_meta_old.tsv 
        mv temp/CoV-{wildcards.id}_download_seqs.fasta  temp/CoV-{wildcards.id}_download_seqs_old.fasta 
        mv temp/CoV-{wildcards.id}_download_aligned_seqs.fasta  temp/CoV-{wildcards.id}_download_aligned_seqs_old.fasta 

        #count number of genbank sequences (minus 1 for header)

        totalSeqs=$(($(cat {output.gen_meta} | wc -l)-1))
        # Figure out message to show user.
        rrun={params.rerun}
        if [ $rrun == "True" ]; then
            #count number of new sequences added (minus 1 for header)
            newSeq=$(($(cat temp/CoV-{wildcards.id}_download_meta_old.tsv | wc -l)-1))
            echo "Existing Genbank files updated with $newSeq new sequences!"
            echo "There is now a total of $totalSeqs Genbank sequences"
        else
            echo "$totalSeqs Genbank sequences stored in database. Reruns will only download new accession numbers."
        fi
        '''
        # totalSeqs=$(($(cat {wildcards.id}/genbank/genbank_meta.tsv | wc -l)-1))
        # newSeq=$(($(cat {wildcards.id}/temp/add_meta.tsv | wc -l)-1))

        #mv {input.download} "temp/CoV-{wildcards.id}_current_download_old.tsv"


#####################################################################################################
#####################################################################################################
#    Bring together ViPR/Genbank (and any manual samples) and do Nextstrain run
#####################################################################################################
#####################################################################################################


##############################
# If want to Concatenate genbank data with Manual at any point - put this here.
# We want this to run if changes to Manual, even if not new Genbank!
###############################

#just copy over metadata so ensure it isn't touched again, and all in one place
rule copy_meta:
    input:
        metadata = rules.make_database.output.gen_meta
    output:
        meta = "{id}/results/metadata.tsv"
    shell:
        '''
        cp {input.metadata} {output.meta}
        '''

# Add extra metadata manually curated, if have it
rule extra_meta:
    input:
        metadata = rules.copy_meta.output.meta,
        extra_meta = files.extra_meta
    output:
        meta = "{id}/results/metadata_extra.tsv"
    shell:
        """
        python scripts/add_meta.py --meta-in {input.metadata} \
            --extra-meta-in {input.extra_meta} \
            --meta-out {output.meta}
        """

#Filter things - different filtering according to different strains
rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - minimum genome length of {params.min_length}
        """
    input:
        sequences = rules.make_database.output.gen_seqs, 
        metadata = rules.extra_meta.output.meta, 
        include = files.include,
        exclude = files.exclude
    output:
        sequences = "{id}/results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 500, #100,
        min_length = 800,
    shell:
        """
        #should filter by host?
        if [ "{wildcards.id}" == "229E" ]; then
            echo "Including only human samples (run {wildcards.id})"
            exclude_where="--exclude-where host!=Human"
            include_where=""
        elif [ "{wildcards.id}" == "SARS" ]; then
            echo "Excluding all bat & mouse samples (run {wildcards.id})"
            exclude_where="--exclude-where host=Bat host=Mouse"
            include_where=""
        elif [ "{wildcards.id}" == "bOC43" ]; then
            echo "Including only human & chimp samples (run {wildcards.id})"
            exclude_where="--exclude-where host!=Human"
            include_where="--include-where host=Chimpanzee"
        else
            echo "Not filtering samples by host"
            exclude_where=""
            include_where=""
        fi

        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            $exclude_where \
            $include_where \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-length {params.min_length}
        """

#Do not align again - takes loads of time, and it's already done. 

##Mask spike gene if wanted.... only done for NL63 currently
rule mask:
    input:
        alignment = rules.filter.output.sequences
    output:
        sequences = "{id}/results/masked.fasta"
    shell:
        """
        if [ "{wildcards.id}" == "NL63" ]; then
            echo "Masking bases from 20500 to 22300"
            echo "(run {wildcards.id})"
            from_to="--mask-from-X-to-X 20500 22300"
        else
            echo "Not masking"
            from_to=""
        fi

        python scripts/mask-alignment.py \
            --alignment {input.alignment} \
            $from_to \
            --output {output.sequences}
        """

rule sub_alignments:
    input:
        alignment = rules.mask.output.sequences
    output:
        sequences = "{id}/results/masked{gene}.fasta"
    run:
        real_gene = wildcards.gene.replace("-","",1)
        boundaries = {
            'replicase':{'Betacoronavirus1':(210,21496), '229E':(293,20568)},
            'last7':{'Betacoronavirus1':(21506,30425), '229E':(20570,26855)}}
        b = boundaries[real_gene][wildcards.id]
        from Bio import SeqIO
        alignment = SeqIO.parse(input[0], "fasta")
        with open(output.sequences, "w") as oh:
            for record in alignment:
                sequence = record.seq.tomutable()
                gene_keep = sequence[b[0]:b[1]]
                sequence[0:len(sequence)] = len(sequence)*"N"
                sequence[b[0]:b[1]] = gene_keep
                record.seq = sequence
                SeqIO.write(record, oh, "fasta")

#Build initial tree
rule tree:
    message: "Building initial tree"
    input:
        alignment = rules.sub_alignments.output.sequences #rules.mask.output.sequences
    output:
        tree = "{id}/results/tree_raw{gene}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

## Rooting & pruning outgroup could go here - excluded for moment.

# Refine tree
rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.sub_alignments.output.sequences, #rules.mask.output.sequences,
        metadata = rules.extra_meta.output.meta
    output:
        tree = "{id}/results/tree{gene}.nwk",
        node_data = "{id}/results/branch_lengths{gene}.json"
    params:
        clock_rate = 0.000459, # estimate taken from MERS via Dudas et al. 2018. eLife.
        clock_std_dev = 0.0003,
        coalescent = "skyline",
        date_inference = "marginal"
    shell:
        """
        #which ones do we allow to estimate rate?
        if [ "{wildcards.id}" == "229E" ] || [ "{wildcards.id}" == "SARS" ] || [ "{wildcards.id}" == "Betacoronavirus1" || [ "{wildcards.id}" == "OC43"]; then
            echo "Estimating clock rate (run {wildcards.id})"
            clock_rate=""
            clock_std_dev=""
        else
            echo "Setting clock rate at {params.clock_rate} with std dev {params.clock_std_dev}"
            echo "(run {wildcards.id})"
            clock_rate="--clock-rate {params.clock_rate}"
            clock_std_dev="--clock-std-dev {params.clock_std_dev}"
        fi

        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            $clock_rate \
            $clock_std_dev \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --date-confidence \
            --no-covariance
        """
        #    --keep-root \
        #--clock-rate {params.clock_rate} \
        #--clock-std-dev {params.clock_std_dev} \

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.sub_alignments.output.sequences #rules.mask.output.sequences
    output:
        node_data = "{id}/results/nt_muts{gene}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.translate_reference
    output:
        node_data = "{id}/results/aa_muts{gene}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.extra_meta.output.meta,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = files.auspice_config,
        colors = files.colors,
        lat_longs = files.lat_longs,
        description = files.description
    output:
        auspice_json = "auspice/CoV_{id}{gene}_spike.json" #rules.all.input.auspice_json
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} \
            --auspice-config {input.auspice_config} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --output {output.auspice_json} \
            --title "Genomic epidemiology of Coronavirus (CoV) {wildcards.id}{wildcards.gene} using data from ViPR & GenBank"
        """
