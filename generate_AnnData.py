"""
Build a cell_line by kmer-count matrix for all cell lines.
Save as AnnData object.
Generate PCA and UMAP projection of cell lines

Requires bedtools, kmer-counter

Usage:
python generate_AnnData.py K out_dir
where K is the size of your desired kmer.
"""

import subprocess
import os
import re
import sys

import anndata as ad
import pandas as pd
import scanpy as sc

def main():
    kmer_size = sys.argv[1]
    out_dir = sys.argv[2]

    # Constants:
    data_file = "AA_CCLE_hg38_aggregated_050323/aggregated_results.csv"
    bed_file = "AA_CCLE_hg38_aggregated_050323.bed"
    fasta_file = os.path.join(out_dir,"aggreaged_ecDNA.fa")
    counts_file = os.path.join(out_dir, f"aggregated_ec_{kmer_size}-mer.counts")
    adata_file = os.path.join(out_dir, f"aggreaged_sample_by_{kmer_size}-mer.h5ad")
    pca_plot = f"_aggreaged_sample_by_{kmer_size}-mer.pca.pdf"
    umap_plot = f"_aggreaged_sample_by_{kmer_size}-mer.umap.pdf"


    if os.path.isfile(bed_file): 
        skip_bed = True
    else: 
        skip_bed = False

    # Load in initial matrix
    if not skip_bed:
        df = loadDF(data_file)

    # Generate bed file for all regions (include length and cell line of origin)
    if not skip_bed:
        generateBed(df, bed_file)

    # Generate fasta file for all beds:
    generateFasta(bed_file, fasta_file)

    # Perform kmer counting from bed file with kmer-counter
    generateKmerCounts(kmer_size, fasta_file, counts_file)

    # Generate sample by kmer dataframe
    sample_by_kmer_df = kmerCountsDataFrame(
        os.path.join(counts_file,"count.txt"), 
        bed_file
    )

    # Build anndata and save file
    adata = buildAnnData(sample_by_kmer_df, adata_file)

    # Generate and save plots
    savePlots(adata, pca_plot, umap_plot, kmer_size)

def loadDF(data_file):
    keep_cols = ['Tissue of origin', "Sample name", 
        "Classification", "Location", 'Captured interval length']
    df = pd.read_csv(
        data_file,
        usecols = keep_cols,
        dtype = {"length":"int", 'Tissue of origin':"str"}
        ).rename(
            columns={
                'Tissue of origin':"tissue_of_origin", 
                "Sample name":"sample_name",
                'Captured interval length':"length"
            }
        ).query(
            "Classification == 'ecDNA'"
        )
    print("Read in df")
    return df

def buildAnnData(sample_by_kmer_df, adata_file):
    adata = ad.AnnData(sample_by_kmer_df)

    tissue_types = []
    for string in adata.obs.index:
        tissue_types.append(string.split("_",1)[1])
    adata.obs["tissue_type"] = tissue_types

    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.pca(adata)

    adata.write(adata_file)
    print(f"Saved {adata_file}")
    return adata

def savePlots(adata, pca_plot, umap_plot, kmer_size):
    # PCA
    sc.pl.pca(
        adata,
        color=["tissue_type"],
        save = pca_plot,
        title = f"PCA from {kmer_size}-mer, Tissue Type",
        show = False
    )

    # UMAP
    sc.pl.umap(
        adata,
        color="tissue_type",
        save = umap_plot,
        title = f"UMAP from {kmer_size}-mer, Tissue Type",
        show = False
    )

def generateKmerCounts(kmer_size, fasta_file, counts_file):
    count_command = ["kmer-counter/kmer-counter", f"--k={kmer_size}", "--fasta", 
        fasta_file ,"--rc", f"--results-dir={counts_file}"]
    subprocess.run(count_command)
    print(f"Built Counts {counts_file}")

def generateFasta(bed_file, fasta_file):
    fasta_command = ["bedtools", "getfasta", "-fo", fasta_file, "-fi", 
        "data/hg38.fa.gz", "-bed", bed_file]
    subprocess.run(fasta_command)
    print(f"Built Fasta {fasta_file}")

def generateBed(df, bed_file):
    regions = []
    for l in zip(df.Location, df.sample_name):
        regions += toBed(l)
    print("Generated regions")

    beregid_ons = "\n".join(regions)
    with open(bed_file, "w") as f:
        f.write(bed_regions)
    print(f"Wrote {bed_file}")

def kmerCountsDataFrame(counts_file, bed_file):
    # Read the input file with tab as the only separator
    df = pd.read_csv(counts_file, 
        header=None, 
        sep="\t"
        )

    # Create a dictionary to store k-mer counts
    kmer_columns = {}

    # Loop through each row to process the k-mer counts
    for idx, row in df.iterrows():
        try:
            kmer_counts = row[1].split(" ")
        except:
            pass
            # Some columns are shorter than the min kmer count. We skip those.
        for kmer_count in kmer_counts:
            try:
                kmer, count = kmer_count.split(':')
                count = int(count)
                if kmer not in kmer_columns:
                    kmer_columns[kmer] = [0] * len(df)
                kmer_columns[kmer][idx] = count
            except:
                pass

    # Convert the kmer_columns dictionary to a DataFrame
    kmer_df = pd.DataFrame(kmer_columns)

    # Add region information to the DataFrame
    kmer_df.insert(0, 'Region', df[0])
    kmer_df["Region"] = kmer_df["Region"].str.replace('>', '')
    kmer_df = kmer_df.set_index("Region")

    # Read the BED file to map regions to sample names
    bed_df = pd.read_csv(bed_file, 
        sep="\t",
        header=None,
        names=["chr", "start", "end", "sample_name"],
        dtype={"start": "str", "end": "str"}
        )

    bed_df["region"] = bed_df["chr"] + ":" +  bed_df["start"] + "-" + bed_df["end"]
    region_sample_map = bed_df.set_index('region')['sample_name'].to_dict()

    # Map regions to sample names and sum k-mer counts by sample
    kmer_df['sample_name'] = kmer_df.index.map(region_sample_map)
    sample_by_kmer = kmer_df.groupby('sample_name').sum()

    return sample_by_kmer

def toBed(inp):
    string_ranges = inp[0]
    sample = inp[1]
    output = []
    string_ranges = re.sub(r"[\[\] \' ]+", "", string_ranges)
    bed_ranges = string_ranges.split(",")
    for br in bed_ranges:
        bed_format_range =  "\t".join(re.split("[-|:]", br))
        bed_line = bed_format_range + "\t" + sample
        output.append(bed_line)
    return output

if __name__ == "__main__":
    main()