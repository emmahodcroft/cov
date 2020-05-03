import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from augur.parse import forbidden_characters

if __name__ == '__main__':

    import argparse
    parser = parser = argparse.ArgumentParser(description='parse blast results',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--blast', help="blast result file")
    parser.add_argument('--meta', help="input meta from genbank")
    parser.add_argument('--seqs', help="input sequences from genbank")
    parser.add_argument('--out_seqs', help="output seqs")
    parser.add_argument('--out_meta', help="output meta")
    parser.add_argument('--match_length', type=int, help="min number of bases to match")
    parser.add_argument('--dup_file', help="file to write duplicate info to")
    args = parser.parse_args()

    #read in blast results and give them a header
    blast = pd.read_csv(args.blast, index_col=False,
        names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"])

    #read in sequences
    fasta_sequences = SeqIO.parse(open(args.seqs), 'fasta')

    # OLD: pick the ones we want - matches of at least 700 bp. 
    # OLD: In initial tests, 600bp added only 47 sequences more, 800bp lost 289 sequences.

    #Now allow users to specify match length...
    matchLen = args.match_length

    keep_blast = list(blast.loc[blast.send - blast.sstart >= matchLen, "qseqid"])

    blast['matchlength'] = blast.apply(lambda row: row.send - row.sstart, axis=1)

    print("{} of {} sequences had some kind of BLAST hit in VP1.".format(len(set(blast.qseqid)), len(list(fasta_sequences))))
    print("{} of the {} sequences have at least {}bp in VP1 and will be included".format(len(keep_blast), 
        len(set(blast.qseqid)),matchLen))

    # These meta strain names are still uncorrected for forbidden chars, but the sequence names have been corrected...
    # so we need to make them comparable by removing the forbidden chars before we can find them in the meta...
    meta = pd.read_csv(args.meta, sep='\t', index_col=False)
    keep_meta_forbid_names = []
    for f in meta.strain:
        name = f
        for c,r in forbidden_characters:
            name=name.replace(c,r)
        meta.loc[meta.strain == f, "qseqid"] = name
        if name in keep_blast:
            keep_meta_forbid_names.append(f)
    keep_meta = meta[meta['strain'].isin(keep_meta_forbid_names)]

    #add match lengths
    meta_matchLen = pd.merge(keep_meta, blast[['qseqid','matchlength']], on=['qseqid'], how='left')

    #WARNING - If there are two matches in VP1 then the TOTAL match length will be wrong.
    #the 'matchlength' will only be for the first match...
    #look for two-match-hits & remove
    dups = meta_matchLen.duplicated(subset='strain', keep='first')
    if sum(dups) != 0:
        print("\nRemoving {} two-match-hits. Note that two-match 'matchlength's will be only for the FIRST match.".format(sum(dups)))
        print("Here are the accessions:")
        print(meta_matchLen[dups].accession)
        meta_matchLen = meta_matchLen[~dups]

    #Need to remove duplicates - sequences with same strain that were sequenced twice.
    #Some debug if want to look into duplicates:
        #   dups = meta_matchLen.duplicated(subset="orig_strain", keep=False)
        #   print(meta_matchLen.loc[dups, ['accession','orig_strain','matchlength']].to_string())
    #First get all of them, to write out, for debugging/checking
    dups2 = meta_matchLen.duplicated(subset="orig_strain", keep=False)
    if not dups2.empty:
        meta_matchLen[dups2].to_csv(args.dup_file, sep='\t', index=False)

    dups2 = meta_matchLen.duplicated(subset="orig_strain", keep="first")
    dup_strains = list(meta_matchLen.loc[dups2,'orig_strain'])
    for dup_s in dup_strains:
        #sort by matchlength (longest at top), take first
        ind = meta_matchLen.loc[meta_matchLen.orig_strain == dup_s].sort_values(by="matchlength", ascending=False).index[0]
        meta_matchLen.drop(ind, inplace=True)
    if dup_strains:
        print("\n{} sequences had different accession numbers but the same strain name.".format(len(dup_strains)))
        print("This is the same sample, sequenced multiple times and submitted under different accessions.")
        print("The longest match will be taken and others discarded.")
        print("Duplicates have been written out to {}".format(args.dup_file))

    over700 = len(meta_matchLen.loc[meta_matchLen.matchlength >= 700, "strain"])
    under700 = len(meta_matchLen.loc[meta_matchLen.matchlength < 700, "strain"])
    print("\nThe final number of sequences to be included is {}".format(len(meta_matchLen)))
    print("{} sequences have >= 700bp, and {} sequences have < 700bp.".format(over700, under700))

    #create new keep sequences list from the meta file.
    keep_sequences = list(meta_matchLen["qseqid"])

    # OLD -- Now find them in the fasta, but only take the matching area of the sequence
    #Now take the whole sequence and let 'augur align' figure this out. Prevent over-snipping at the ends of
    #matching regions. 
    fasta_sequences = SeqIO.parse(open(args.seqs), 'fasta')
    keepSeqs = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name in keep_sequences:#keep_blast: 
            #these dont work without removing dups from blast via: blast = blast[~blast.duplicated(subset='qseqid', keep='first')]
            #strt = int(blast.loc[blast.qseqid == name].qstart-1)
            #en = int(blast.loc[blast.qseqid == name].qend)
            vp1_seq = sequence #[strt:en]
            keepSeqs.append(SeqRecord(Seq(vp1_seq), id=name, name=name, description=""))
    SeqIO.write(keepSeqs, args.out_seqs, "fasta")

    #remove the added qseqid column
    meta_matchLen = meta_matchLen.drop("qseqid", axis=1)
    #write out only those we are taking
    meta_matchLen.to_csv(args.out_meta, sep='\t', index=False)


