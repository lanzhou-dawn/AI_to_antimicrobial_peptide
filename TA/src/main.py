import pandas as pd
import sys 
from genbank import gen_bank
from opern import opern_pred
from TAfilter import TA_filter
path = sys.argv[1]
def main():
    #提取基因组信息
    genomeInfo = gen_bank(path)
    id = genomeInfo[2]
    data = genomeInfo[0]
    genome_path = r'../tmp/{}_genomeInfo.csv'.format(id)
    opern_path = r'../tmp/{}_opern.csv'.format(id)
    TA_path = r'../tmp/{}_TA.csv'.format(id)
    TA_complex = r'../tmp/{}_TA_complex.csv'.format(id)
    TA_fasta = r'../out/{}_TA_complex.fasta'.format(id)
    df = pd.DataFrame(data)
    df.to_csv(genome_path, sep= '\t')
    
    #预测

    predictions = opern_pred(path)


    df1 = pd.DataFrame(predictions)
    df1.to_csv(opern_path, sep= '\t')
    
    #TA筛选
    df2 = TA_filter(genome_path,opern_path)
    df2.to_csv(TA_path, sep= '\t',index=0)
    
    #TA序列合并
    df3 = pd.read_csv(TA_path, usecols=[0,1,12,23], sep='\t')
    df3.to_csv(TA_complex, sep='\t',index=0)

    f1 = open(TA_complex,'r').readlines()#需要整理的文件
    f2 = open(TA_fasta, 'w')#整理之后的文件
    for i in f1[1:]:#有表头，如果无表头则将1：去掉
        pro_id = i.split('\t')[0] + '_' + i.split('\t')[1]
        fa = i.strip('\n').split('\t')[2] + 'NNN' + i.strip('\n').split('\t')[3]
        f2.write('>'+pro_id+' '+id+'\n'+fa+'\n')
    f2.close()
#主函数
if __name__ == '__main__':
    main()
    