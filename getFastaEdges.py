
def main():
    input_fasta = "Lung_ecDNA.fa"
    output_fasta = "Lung_ecDNA.trimmed.fa"
    with open(input_fasta,"r") as f:
        with open(output_fasta, "w") as of:
            for line in f.readlines():
                out_line = line
                if not line.startswith(">"):
                    out_line = out_line[:20] + "$" + out_line[-21:]
                of.write(out_line)

if __name__ == "__main__":
    main()