#!/usr/bin/env python
#coding=utf-8

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import re
import sys
import argparse
import subprocess as sp
import shutil

markerscanner = '/home/hjz1/hwz_biosoftware/test_metagenome/AMPHORA2-master/Scripts/MarkerScanner.pl'
mafft = '/home/hjz1/miniconda2/envs/hu_metagenome/bin/mafft'
trimal = '/home/hjz1/miniconda2/envs/hu_metagenome/bin/trimal'
raxmlHPC='/home/hjz1/miniconda2/envs/hu_metagenome/bin/raxmlHPC-PTHREADS-SSE3'
fasttree = '/home/hjz1/miniconda2/envs/hu_metagenome/bin/fasttree'

################# functions ###############
#获取文件名字
def fname(f):
    file_name = os.path.splitext(f)[0]
    return file_name
#创建文件夹，已存在则跳过
def mxmkdir(path):
    '''
    According to the provided path, a new dir including its subdir will be created unless the dir already exists. 
    '''
    path=path.strip()
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False
#串联序列
def concatenate_seqs(indir,outfiles):
    flst = os.listdir(indir)
    #获取所有基因组名称
    genes_lst = []
    for f in flst:
        fpath = os.path.join(indir,f)
        fseqs = SeqIO.parse(fpath,"fasta")
        for fseq in fseqs:
            if fseq.id not in genes_lst:
                genes_lst.append(fseq.id)
    temp_path = os.path.join(indir,"temp")
    
    # 生成每个基因组的保守蛋白串联序列文件
    mxmkdir(temp_path)
    for gene in genes_lst:
        newseq = Seq("")
        for i in range(len(flst)):
            fpath = os.path.join(indir,flst[i])
            fseq_dict = SeqIO.to_dict(SeqIO.parse(fpath,"fasta"))
            if gene in fseq_dict.keys():
                newseq = newseq+fseq_dict[gene].seq
            else:
                added_seq = "?"*len(fseq_dict[fseq_dict.keys()[0]].seq)
                newseq = newseq+Seq(added_seq)
        newseq_rec = SeqRecord(seq=newseq,id = gene,description=gene,name = gene)
        gene_path = os.path.join(temp_path,gene+".allpep.fasta")
        SeqIO.write(newseq_rec,gene_path,"fasta")
        print("{}的保守蛋白串联序列文件{}构建完成！".format(gene,gene_path))
        
    # 合并所有基因组到同一个文件中
    genome_flst = os.listdir(temp_path)
    with open(outfiles,'ab') as out_handle:
        for f in genome_flst:
            fpath = os.path.join(temp_path,f)
            subreads = SeqIO.index(fpath,"fasta")
            for rec in subreads:
                out_handle.write(subreads.get_raw(rec))
    print("所有串联基因组已合并至{}".format(outfiles))
    
    
################prameters ######################

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, help='input pep folder')
parser.add_argument('-o', required=True, help='output folder')
args = vars(parser.parse_args())
inputdir = args['i']
outputdir = args['o']

################# Scripts #######################
curDir = os.getcwd()
inputdir = os.path.join(curDir,inputdir)
#输出文件夹
out_peps=os.path.join(curDir,outputdir,"peps")
out_merged_peps=os.path.join(curDir,outputdir,"Merged_peps")
out_aligned_peps=os.path.join(curDir,outputdir,"aligned_peps")
out_trimed_peps=os.path.join(curDir,outputdir,"trimed_peps")
out_trees=os.path.join(curDir,outputdir,"trees")

####Part1, 抽取保守蛋白序列，工具markerscanner
print("""
###########################################################################################
                  @@@@@@@ Step 1. 抽取31个细菌高保守蛋白序列....@@@@@@@@
工具：AMPHORA2 software
高保守蛋白包括：
danG, frr, infC, nusA, pgk, pyrG, rplA, rplB, rplC, rplD, rplE, rplF, rplK, rplL, rplM, rplN,
rplP, rplS, rplT, rpmA, rpoB, rpsB, rpsC, rpsE, rpsI, rpsJ, rpsK, rpsM, rpsS, smpB and tsf
输出文件夹：{}
###########################################################################################
""".format(out_peps))
fasta_list = os.listdir(inputdir)

for fasta_f in fasta_list:
    print("现在抽取基因组{}中的细菌高保守蛋白序列.....请耐心....".format(fasta_f))
    #后缀名检查
    if os.path.splitext(fasta_f)[1].lower() not in ['.fasta','.fas','.fna']:
        print("{} 不是fasta格式的序列文件，将被删除，避免后续过程文件冲突!".format(fasta_f))
        os.remove(os.path.join(inputdir,fasta_f))
        continue
    else:
        # 获取序列文件的名字作为后续合并时序列的名字
        fastaname = os.path.splitext(fasta_f)[0]
        #进入工作文件夹后，调用AMPHORA2程序
        os.chdir(inputdir)
        sp.call(['perl',markerscanner,'-Bacteria','-DNA',fasta_f])
        #对程序生成的文件进行处理
        peplist = os.listdir(os.getcwd())
        for pep in peplist:
            pepname = os.path.splitext(pep)[0]
            pepext = os.path.splitext(pep)[1]
            # 非pep文件跳过
            if not pepext == '.pep':
                continue
            else:
                #向文件名中增加序列来源识别，并且使用基因组名称替换pep文件中的序列名称
                new_pep_name = pepname+'.'+fastaname+'.pep'
                os.rename(pep,new_pep_name)
                with open(new_pep_name,'r') as f_r:
                    lines = f_r.readlines()
                with open(new_pep_name,'w') as f_w:
                    for line in lines:
                        if line.startswith('>'):
                            line = re.sub('>.+\n','>'+fastaname+'\n',line)
                        f_w.write(line)
                #检查拷贝数，如果为多拷贝，则串联拷贝。
                fseq = SeqIO.parse(new_pep_name,"fasta")
                fseq = list(fseq)
                if len(fseq) > 1:
                    print("{}的高保守蛋白{}存在{}个拷贝，下面将进行多拷贝的拼接...".format(fasta_f,pepname,len(fseq)))
                    ids = fseq[0].id
                    names = fseq[0].name
                    descriptions = fseq[0].description
                    newseq = SeqRecord(seq=Seq(""),id=ids,name=names,description=descriptions)
                    for rec in fseq:
                        newseq = newseq+rec
                    SeqIO.write(newseq,new_pep_name,"fasta")
                #移动文件至相应的文件夹
                new_pep_name_old_path = os.path.join(inputdir,new_pep_name)
                mxmkdir(os.path.join(out_peps,pepname))
                new_pep_name_new_path = os.path.join(out_peps,pepname,new_pep_name)
                shutil.move(new_pep_name_old_path,new_pep_name_new_path)
        os.chdir(curDir)
print("""
###########################################################################################
               @@@@@@@ Step 2. 合并各个保守蛋白并进行比对齐平....@@@@@@@@                  
工具：MAFFT比对，TrimAI齐平
输出文件夹:
合并序列文件夹：{}
比对后序列文件夹：{}
齐平后序列文件夹：{}
###########################################################################################
""".format(out_merged_peps,out_aligned_peps,out_trimed_peps))
os.chdir(curDir)
#####Part 2 按基因合并序列并分别比对
peplst = os.listdir(out_peps)
mxmkdir(out_merged_peps)
#按照基因合并序列文件。
for pepdir in peplst:
    peppath = os.path.join(out_peps,pepdir)
    flst = os.listdir(peppath)
    merged_f = os.path.join(out_merged_peps,pepdir+".fasta")
    with open(merged_f,'ab') as out_handle:
        for f in flst:
            fpath = os.path.join(peppath,f)
            subreads = SeqIO.index(fpath,"fasta")
            for rec in subreads:
                num=out_handle.write(subreads.get_raw(rec))
    print('{}合并完成，合并后文件 {} !'.format(pepdir,merged_f))

mxmkdir(out_aligned_peps)
mxmkdir(out_trimed_peps)
peps_for_align = os.listdir(out_merged_peps)
os.chdir(out_merged_peps)
#比对切平，并且移动文件
print("现在对每个高保守蛋白序列进行分别比对切平")
for peps in peps_for_align:
    pepname = fname(peps)
    peppath = os.path.join(os.getcwd(),peps)
    pep_align = pepname+'.align.pep'
    pep_trim = pepname+'.trimed.fasta'
    print("""
&&&&&&&&&&&&&&&&&&&&&&&&&&比对齐平{},输入文件{}&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    """.format(pepname,peppath))
    sp.call(['mafft --maxiterate 1000 --globalpair {} > {}'.format(peps,pep_align)],shell = True)
    sp.call(['trimal -in {} -out {} -fasta -automated1'.format(pep_align,pep_trim)],shell = True)
    shutil.move(os.path.join(out_merged_peps,pep_align),os.path.join(out_aligned_peps,pep_align))
    shutil.move(os.path.join(out_merged_peps,pep_trim),os.path.join(out_trimed_peps,pep_trim))

os.chdir(curDir)

print("""
###########################################################################################
               @@@@@@@ Step 3. 串联保守蛋白，构建进化树....@@@@@@@@                  
工具：fasttree
###########################################################################################
""")
mxmkdir(out_trees)
tree_input_fasta = os.path.join(out_trees,"tree_input.fasta")
concatenate_seqs(out_trimed_peps,tree_input_fasta)
os.chdir(out_trees)
print("""
@@@@@@@@@@@@@@ 建树软件选择 @@@@@@@@@@@@@@
[fasttree] 使用fasttree 建树，参数为 '-lg -gamma'
[raxml] 使用raxmlHPC-PTHREADS-SSE3建树, 参数为 '-f a -n tre -m PROTGAMMALGX -x 1000 -# 1000 -p 1000 -T 40'
后续版本会提供更多参数选择
[q] 退出程序，自定义建树
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
""")
tree_software = ''
while tree_software not in ['fasttree','raxml','q']:
    tree_software = input("请输入 fasttree/raxml/q 进行选择---> ")

if tree_software == 'fasttree':
    fasttreeCMD = ['fasttree -lg -gamma {} > lg_gamma.tree'.format(tree_input_fasta)]
    sp.call(fasttreeCMD,shell=True)
elif tree_software == 'raxml':
    from Bio import AlignIO
    AlignIO.convert("tree_input.fasta","fasta","tree_input.phy","phylip-relaxed")
    raxmlCMD = ['raxmlHPC-PTHREADS-SSE3 -f a -s tree_input.phy -n tre -m PROTGAMMALGX -x 1000 -# 1000 -p 1000 -T 40']
    sp.call(raxmlCMD,shell=True)
elif tree_software == 'q':
    print('运行结束，祝愉快！')
    sys.exit()
