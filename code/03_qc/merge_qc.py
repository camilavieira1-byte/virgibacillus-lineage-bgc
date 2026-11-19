import pandas as pd
base = "/PATH/TO/virgibacillus_v2"
meta = pd.read_csv(f"{base}/00_metadata/metadata_raw_flagged.tsv", sep="\t")
cm2 = pd.read_csv(f"{base}/03_qc/checkm2/quality_report.tsv", sep="\t")[
    ["Name","Completeness","Contamination","Coding_Density","Genome_Size","Total_Contigs","Contig_N50","GC_Content"]]
cm2 = cm2.rename(columns={"Name":"accession"})
q = pd.read_csv(f"{base}/03_qc/quast/transposed_report.tsv", sep="\t")
q = q.rename(columns={"Assembly":"accession","# contigs":"quast_contigs","N50":"quast_N50","Largest contig":"largest_contig"})
q = q[["accession","quast_contigs","quast_N50","largest_contig"]]
df = meta.merge(cm2, on="accession", how="inner").merge(q, on="accession", how="left")

df["hq_mimag"]   = (df.Completeness >= 90) & (df.Contamination < 5)
df["mq_or_better"] = (df.Completeness >= 50) & (df.Contamination < 10)
mu, sd = df.Genome_Size.mean(), df.Genome_Size.std()
df["size_z"] = (df.Genome_Size - mu) / sd
df["size_outlier"] = df.size_z.abs() > 3
df.to_csv(f"{base}/03_qc/qc_merged.tsv", sep="\t", index=False)

print("Genomes:", len(df))
print("High quality (>=90% / <5%):", int(df.hq_mimag.sum()))
print("Medium or better:", int(df.mq_or_better.sum()))
print("Size outliers (|z|>3):", int(df.size_outlier.sum()))
print(df[["Completeness","Contamination","Genome_Size","Total_Contigs","quast_N50"]].describe().round(2))
print("\nLowest completeness / highest contamination:")
print(df.sort_values("Completeness").head(5)[["accession","Completeness","Contamination","Genome_Size"]])
print(df.sort_values("Contamination", ascending=False).head(5)[["accession","Completeness","Contamination","Genome_Size"]])
