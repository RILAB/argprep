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
                    allele_str = ''
                    for i in range(9,len(line)):
                        ref_read_depth = line[i].split(':')[1].split(',')[0]
                        if ref_read_depth == '0':
                            allele = '1'
                            allele_str = allele_str + '\t' + allele
                        else:
                            allele = "0"
                            allele_str = allele_str + '\t' + allele
                    newline = '\t'.join(line[:2])
                    newline1 = '\t'.join(line[5:9])
                    #name = 'GT'
                    ID = f'snp_{line[0]}_{line[1]}'
                    nf.write(f'{newline}\t{ID}\t{line[3]}\t{line[4].replace(",<NON_REF>","")}\t{newline1}{allele_str}\n')

vcf4singer(filename)
