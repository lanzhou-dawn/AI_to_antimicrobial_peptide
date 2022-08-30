import pandas as pd
from Bio import SeqIO
from Bio import SeqFeature



def gen_bank(path :str):
    data = []
    data_1 = []
    record = list(SeqIO.parse(path, "genbank"))[0]
    id = record.id

    for gene in record.features[1:] :
        if not gene.type == 'CDS' :
            continue


        
        gen_star = gene.location.start
        gen_end = gene.location.end
        gen_strand = gene.location.strand
        gen_type = gene.type
        try:
            protein_seq = gene.qualifiers["translation"][0]
        except KeyError:
            protein_seq = None
            
        try:
            gen_protein_id = gene.qualifiers['protein_id'][0]
        except KeyError:
            gen_protein_id = None
        try:
            gen_id = gene.qualifiers['locus_tag'][0]
        except KeyError:
            gen_id = None   
        

        try:
            gen_product = gene.qualifiers['product'][0]
        except KeyError:
            gen_product = None
            
        try:
            gene_name = gene.qualifiers['gene'][0]
        except KeyError:
            gene_name = None
              
        data.append({'gen_id':gen_id, 
                     'gen_protein_id':gen_protein_id,
                     'gen_star': gen_star,
                     'gen_end': gen_end,
                     'gen_strand': gen_strand,
                     'gen_type': gen_type,
                     
                     'gen_product': gen_product,
                     'gene_name': gene_name,
                     'protein_seq':protein_seq
                     })
        data_1.append({'left': gen_star, 'right': gen_end, 'strand': gen_strand, 'gene_id': gen_id})

        
    return data, data_1, id  

if __name__ == '__main__':
    gen_bank('genomic\GCA_002157365.2_ASM215736v2_genomic.gbff')
'''
    data = gen_bank(path_1)[0]
    df = pd.DataFrame(data)
    df.to_csv('genomeInfo.csv', sep= '\t')
    '''