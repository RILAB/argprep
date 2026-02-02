import sys
import gzip
filename = sys.argv[1]

if filename.endswith(".gz"):
    newname = filename.replace(".gz", ".format.vcf")
else:
    newname = f'{filename}.format.vcf'

def open_maybe_gzip(filename, mode="rt"):
    """
    Open plain text or gzipped files transparently.
    """
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    return open(filename, mode, encoding="utf-8")

def vcf4singer(filename):
    """
           Convert the .clean file sample genotype field into 0 or 1
    (assuming homozygous samples only).

           Parameters:
               sample_field (str): e.g. ".:1,0,0:1"
               format_field (str): e.g. "GT:AD:DP"
               we are using the AD fild (the reads mapped to ref (:1,0,0:))
               to infer the genotyping, so need check whether you have this AD filed

           Returns:
               str: inferred genotype (0 or 1)
           """
    with open_maybe_gzip(filename,"rt") as f:
        with open(newname,'w',encoding="utf-8") as nf:
            for line in f:
                if line.startswith('#'):
                    nf.write(f'{line}')
                else:
                    line = line.strip('\n').split('\t')
     		        format_field = line[8]
                    format_keys = format_field.split(':')
                    try:
                        ad_idx = format_keys.index("AD")
                    except ValueError:
                        ad_idx = None
                    allele_str = ''
                    for i in range(9,len(line)):
                        sample = line[i]
                        allele = "."
                        if sample not in (".", "./.", ".|."):
                            sample_parts = sample.split(':')
                            if ad_idx is not None and ad_idx < len(sample_parts):
                                ad_field = sample_parts[ad_idx]
                                ref_read_depth = ad_field.split(',')[0] if ad_field not in (".", "") else "."
                                if ref_read_depth != ".":
                                    allele = '1' if ref_read_depth == '0' else '0'
                        allele_str = allele_str + '\t' + allele

                    newline = '\t'.join(line[:2])
		            newline1 = '\t'.join(line[5:8]) + '\tGT'
                    #newline1 = '\t'.join(line[5:9])
                    ID = f'snp_{line[0]}_{line[1]}'
                    nf.write(f'{newline}\t{ID}\t{line[3]}\t{line[4].replace(",<NON_REF>","")}\t{newline1}{allele_str}\n')

vcf4singer(filename)
