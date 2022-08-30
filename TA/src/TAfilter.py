import pandas as pd
import numpy as np
def TA_filter(genome_path:str ,opern_path:str):
    #0. 导入数据
    metadf = pd.read_csv(opern_path, sep = '\t' ,header=0,index_col=0)
    metagenome = pd.read_csv(genome_path, sep = '\t',header=0,index_col=0)
    #print(metadf)
    #print(metagenome)
    #1. 整理数据
    #1.1 蛋白大小、注释信息字符化

    metagenome['gen_product'] = metagenome['gen_product'].map(lambda x : str(x))
    metagenome['size'] = metagenome['gen_end'].astype(int) - metagenome['gen_star'].astype(int)
    metagenome['size'] = (metagenome['size'].abs())/3 

    #print(metagenome)

    #1.2 合并数据
    metadf2 = metadf[['Gene2','dist']]
    metadf.drop(columns=['dist'],inplace=True)


    df1 = pd.merge(metadf, metagenome ,how= 'inner' ,left_on= 'Gene1', right_on= 'gen_id' )
    df2 = pd.merge(metadf2, metagenome ,how= 'inner' ,left_on= 'Gene2', right_on= 'gen_id' )
    df = pd.merge(df1, df2 ,how= 'inner', on= 'Gene2', suffixes = ("_g1", "_g2"))


    #2. 筛选数据
    #(i) the two genes, each consisting of less than 200 amino acid (aa) residues, either overlap or are less than 30 bases apart
    df = df[(df['pred'] == True )]
    #使用集合
    set_gene1_2 = set(df['Gene1']).intersection(set(df['Gene2']))
    dict_Gene1_2 = df[["Gene2", "Gene1"]].set_index("Gene2").to_dict(orient='dict')["Gene1"]
    list_gene1 = list(set_gene1_2)

    for i in list_gene1 :
        set_gene1_2.add(dict_Gene1_2[i])




    df = df.set_index('Gene1',drop=True)
    df = df.drop(list(set_gene1_2))
    df.reset_index(inplace = True) 

    df = df[(df['dist'] <= 30)     & 
            (df['size_g1'] < 200 ) &
            (df['size_g2'] < 200)
            ]

    #(ii) the genes consist of an operon encoding two proteins

    #df = df[(df['pOp'] >= 0.98 )]
    #数据转化为字符
    #df=df.astype(str)
    #print(df.dtypes)
    #(iii) both genes encode a hypothetical protein
    #desc
    df = df[
        (df['gen_product_g1'].str.contains('hypothetical protein', case = False) ) |
        (df['gen_product_g2'].str.contains('hypothetical protein', case = False) ) 
    ]
    return df
