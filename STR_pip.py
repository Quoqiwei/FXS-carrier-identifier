import subprocess
import os 
import re
import argparse
import sys
import pandas as pd
from pandarallel import pandarallel
from Bio import SeqIO
from joblib import Parallel, delayed
parser = argparse.ArgumentParser(description='Data Analysis Pipline')
parser.add_argument('--fastq','-fq',type=str,default=False,help='input fastq file')
parser.add_argument('--fast5','-f5',type=str,default=False,help='input fast5 file')
parser.add_argument('--list','-fl',type=str,default=False,help=' list of fastq or fast5 file')
parser.add_argument('--tarlen','-l',type=int,default=40,help='length of prefix and suffix sequences')
parser.add_argument('--qvalue','-q',type=int,default=4,help='minimum mapping quality value for prefix and suffix sequences')
parser.add_argument('--cov','-c',type=int,default=0,help='coverage of the target area for read')
parser.add_argument('--tarfile','-tf',type=str,default='/home/gqw/STR_soft/STR_test/ref/tar.bed',help='position of amplified sequence')
parser.add_argument('--ref','-r',type=str,default='/home/gqw/STR_soft/STR_test/ref/genome.fa',help='reference genome')
parser.add_argument('--mode','-m',type=str,default='CGG',help='STR mode[default:CGG]')
parser.add_argument('--threads','-t',type=int,default=10,help='threads number')
parser.add_argument('--output','-o',type=str,default='STR_analysisfile',help='output file')
parser.add_argument('--samplename','-n',type=str,default='sample_pass',help='sample name')
parser.add_argument('--barcodelen','-bl',type=int,default=6,help='the length of identifier sequence')
parser.add_argument('--barseq','-bs',type=str,help='prefix and suffix sequence file')
parser.add_argument('--ifbar','-ib',type=str,default='bar', help='whether the sample has identifier sequence')
parser.add_argument('--barcodemis','-bm',type=int,default=1, help='maximum number of mismatched base for identifier sequence')
parser.add_argument('--barcodefa','-bf',type=str,default='/home/gqw/STR_soft/STR_test/script/test_barcode.fa',help='identifier sequence file')
parser.add_argument('--STRmotif','-sf',type=str,default='nomotif', help='AGG interruption')
parser.add_argument('--STRMut','-sm',type=str,default='noMut')
parser.add_argument('--lowper','-lp',type=float,default=0, help='threshold of premutation and full mutation allele')
parser.add_argument('--lowper1','-lp1',type=float,default=0, help=' threshold of AGG interruption')
parser.add_argument('--lowper53','-l53',type= int,default=10000, help='threshold of wildtype allele')
parser.add_argument('--shift_num','-snm',type=int,default=0, help='the base number of intervald between prefix and suffix sequences and CGG repeat sequence')
parser.add_argument('--gc_content','-gc',type=float,default=90, help='the GC content between prefix and suffix sequences')
parser.add_argument('--SE','-se',type=str,default='/home/gqw/STR_soft/STR_test/script/SE.tsv', help='peak range')
argv = parser.parse_args()
#-----Pass in parameters
sn = argv.samplename
if argv.fastq:
    dfile = os.path.abspath(argv.fastq)
    fmode='fq'
elif argv.fast5:
    dfile = os.path.abspath(argv.fast5)
    fmode='f5'
else:
    dfile = os.path.abspath(argv.list)
    fmode='list'
def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

sumdict = {}
tarfile = argv.tarfile
SE = argv.SE
lowper = argv.lowper
lowper1 = argv.lowper1
lowper53= argv.lowper53
qvalue = argv.qvalue
tarlen = argv.tarlen
gc_content = argv.gc_content
cov = argv.cov/100
ref = argv.ref
mode = argv.mode
STRmotif = argv.STRmotif
Mutmotif = argv.STRMut
threads = argv.threads
ofn = argv.output
Ifbar = argv.ifbar 
barmis = argv.barcodemis
barcodefa = argv.barcodefa
shift_num = argv.shift_num
if os.path.isdir(ofn):
    print('The output file already exists.')
    sys.exit()
else:
    os.mkdir(ofn)
wkdir = os.path.abspath(ofn)
Blen = argv.barcodelen
tarbed = '/home/gqw/STR_soft/STR_test/ref/tar.bed'
def get_bar(tmpdf):
    tbardb = tmpdf
    tmpdb = pd.DataFrame(columns=['readsID','barpos','barcodenum'])
    if tbardb.shape[0] == 2:
        if len(set(tbardb['strand'].tolist())) == 2:
            frontpos = int(tbardb[tbardb['strand']=='+']['start'].tolist()[0])
            frontend = int(tbardb[tbardb['strand']=='+']['end'].tolist()[0])
            pregion = f'{frontpos}-{frontend}'
            rearpos = int(tbardb[tbardb['strand']=='-']['start'].tolist()[0])
            rearend = int(tbardb[tbardb['strand']=='-']['end'].tolist()[0])
            rregion = f'{rearpos}-{rearend}'
            return pd.DataFrame({'readsID':tbardb['seqID'].tolist()[0],'barcodenum':2,'class':'pairedbar','Plus_region':pregion,'Rev_region':rregion},index=[0])
        else:
            pregionlist = []
            rregionlist = []
            if tbardb['strand'].tolist()[0] == '+':           
                for i in tbardb[tbardb['strand']=='+'].index.tolist():
                    ttbardb = tbardb.loc[tbardb.index==i,:]
                    frontpos = int(ttbardb[ttbardb['strand']=='+']['start'].tolist()[0])
                    frontend = int(ttbardb[ttbardb['strand']=='+']['end'].tolist()[0])
                    pregion = f'{frontpos}-{frontend}'
                    pregionlist.append(pregion)
            else:       
                for i in tbardb[tbardb['strand']=='-'].index.tolist():
                    ttbardb = tbardb.loc[tbardb.index==i,:]
                    rearpos = int(ttbardb[ttbardb['strand']=='-']['start'].tolist()[0])
                    rearend = int(ttbardb[ttbardb['strand']=='-']['end'].tolist()[0])
                    rregion = f'{rearpos}-{rearend}'
                    rregionlist.append(rregion)
            mpregion = ';'.join(pregionlist)
            mrregion = ';'.join(rregionlist)
            return  pd.DataFrame({'readsID':tbardb['seqID'].tolist()[0],'barcodenum':2,'class':'errorbar','Plus_region':mpregion,'Rev_region':mrregion},index=[0])
    elif tbardb.shape[0] > 2:
        rregionlist = []
        pregionlist = []
        for i in tbardb.index.tolist():
            ttbardb = tbardb.loc[tbardb.index==i,:]
            if ttbardb['strand'].tolist()[0] == '+':
                frontpos = int(ttbardb[ttbardb['strand']=='+']['start'].tolist()[0])
                frontend = int(ttbardb[ttbardb['strand']=='+']['end'].tolist()[0])
                pregion = f'{frontpos}-{frontend}'
                pregionlist.append(pregion)
            else:
                rearpos = int(ttbardb[ttbardb['strand']=='-']['start'].tolist()[0])
                rearend = int(ttbardb[ttbardb['strand']=='-']['end'].tolist()[0])
                rregion = f'{rearpos}-{rearend}'
                rregionlist.append(rregion)
        mpregion = ';'.join(pregionlist)
        mrregion = ';'.join(rregionlist)
        return  pd.DataFrame({'readsID':tbardb['seqID'].tolist()[0],'barcodenum':tbardb.shape[0],'class':'multibar','Plus_region':mpregion,'Rev_region':mrregion},index=[0])
    else:
        if tbardb['strand'].tolist()[0] == '+':
            frontpos = int(tbardb[tbardb['strand']=='+']['start'].tolist()[0])
            frontend = int(tbardb[tbardb['strand']=='+']['end'].tolist()[0])
            pregion = f'{frontpos}-{frontend}'
            return  pd.DataFrame({'readsID':tbardb['seqID'].tolist()[0],'barcodenum':1,'class':'singlebar','Plus_region':pregion,'Rev_region':'-'},index=[0])

        else:
            rearpos = int(tbardb[tbardb['strand']=='-']['start'].tolist()[0])
            rearend = int(tbardb[tbardb['strand']=='-']['end'].tolist()[0])
            rregion = f'{rearpos}-{rearend}'
            return  pd.DataFrame({'readsID':tbardb['seqID'].tolist()[0],'barcodenum':1,'class':'singlebar','Plus_region':'-','Rev_region':rregion},index=[0])

def applyParallel(dfGrouped, func):
    with  Parallel(n_jobs=10) as parallel:
        retLst = parallel(delayed(func)(group) for name, group in dfGrouped)
    return pd.concat(retLst)
def applyParallel1(dfGrouped, func):
    with Parallel(n_jobs=10) as parallel1:
        retLst = parallel1(delayed(func)(group) for name, group in dfGrouped)
    return retLst

def custom_barcoding(inf,ofn,Pre,barfile,mism,ifbar):
    if not os.path.isdir(ofn):
        os.makedirs(ofn)
    subprocess.run(f'seqkit fx2tab -n -i {inf} > {ofn}/AllID.txt',shell=True)
    subprocess.run(f'seqkit fx2tab -n -l -i {inf} > {ofn}/All_IDlength.tsv',shell=True)
    IDldb = pd.read_table(f'{ofn}/All_IDlength.tsv',header=None,names=['seqID','seqLength'])
    if ifbar == 'bar':
        
        subprocess.run(f'seqkit fx2tab {barfile} > {ofn}/bar.tsv',shell=True)
        with open(f'{ofn}/bar.tsv') as f:
            for line in f:
                barID = line.split('\t')[0].strip()
                barfa = line.split('\t')[1].strip()
                if mism!=0:
                    subprocess.run(f'seqkit locate -m {mism} -p {barfa} {inf} -j 10 > {ofn}/{barID}_seqgrep.tsv',shell=True)
                else:
                    subprocess.run(f'seqkit locate -p {barfa} {inf} > {ofn}/{barID}_seqgrep.tsv',shell=True)
                bardb = pd.read_table(f'{ofn}/{barID}_seqgrep.tsv')
                bardb = bardb.merge(IDldb,on='seqID')
                bardb.to_csv(f'{ofn}/{barID}_seqgrepleng.tsv',sep='\t',index=False)
               
                if bardb.shape[0] >= 1:
                    
                    barlocdb = applyParallel(bardb.groupby('seqID'),get_bar)
                    barlocdb.to_csv(f'{ofn}/seqID_bar_{barID}.tsv',sep='\t',index=False)
                    print(barID)
                    print(barlocdb.value_counts(['class']))
                    barlocdb[barlocdb['class'].isin(['singlebar','pairedbar'])]['readsID'].to_csv(f'{ofn}/{Pre}_{barID}_barcoding_ID.txt',sep='\t',header=None,index=None)
                    barlocdb[barlocdb['class']=='pairedbar']['readsID'].to_csv(f'{ofn}/{Pre}_{barID}_bothID.txt',sep='\t',index=None,header=None)
                    barlocdb[barlocdb['class']=='singlebar']['readsID'].to_csv(f'{ofn}/{Pre}_{barID}_nobothID.txt',sep='\t',index=None,header=None)
                else:
                    print(f'{barID} does not detected.')


    else:
        subprocess.run(f'cp {ofn}/AllID.txt {ofn}/{Pre}_nobar_barcoding_ID.txt',shell=True)
#------Extract the amplified sequence------
def tar_select(tf,ref,ofn,tl,nt):
    with open(f'{wkdir}/task.log','a') as f1:
        os.mkdir(f'{ofn}/barcode')
        if not is_fasta(tf):
            with open(tf) as f:
                for line in f:
                    line = line.strip().split('\t')
            seqname = f'{line[0]}:{line[1]}-{line[2]}'
            try:
                subprocess.run(f'bedtools getfasta -fi {ref} -bed {tf} -fo {ofn}/barcode/Target.fa',stdout=f1,stderr=f1,shell=True)
            except LookupError:
                print('There is no common sequence between the amplified sequence and reference genome.')
        else:
            with open(tf) as f:
                for line in f:
                    if line.startswith('>'):
                        seqname = line.strip().replace('>','')
            subprocess.run(f'cp {tf} {ofn}/barcode/Target.fa',shell=True)

        if not argv.barseq:
            rear_bed = f'{seqname}\t191\t191+{tl}-1'
            subprocess.run(f'seqkit subseq -r {130-tl+1}:130 {ofn}/barcode/Target.fa -j {nt} -w0 > {ofn}/barcode/Tar_{tl}.fa',stdout=f1,stderr=f1,shell=True)
            subprocess.run(f'''sed -i 's/{seqname}/front/g' {ofn}/barcode/Tar_{tl}.fa''',stdout=f1,stderr=f1,shell=True)
            subprocess.run(f'seqkit subseq -r 191:{191+tl-1} {ofn}/barcode/Target.fa -j {nt} -w0 >> {ofn}/barcode/Tar_{tl}.fa',stdout=f1,stderr=f1,shell=True)
            subprocess.run(f'''sed -i 's/{seqname}/rear/g' {ofn}/barcode/Tar_{tl}.fa''',stdout=f1,stderr=f1,shell=True)
        else:
            def comparelist(list1,list2):
                if list1[0] > max(list2):
                    return list2,list1
                elif list2[0] > max(list1):
                    return list1,list2
                else:
                    print('Stop analysis because there is an overlap on both sides of the amplified sequence.')
                    sys.exit()

            with open(argv.barseq) as bf:
                tmpl=[]
                tmpdict = {}
                n=1
                for line in bf:
                    if not line.startswith('>'):
                        seq = line.strip()
                        subprocess.run(f'seqkit locate -p {seq} {ofn}/barcode/Target.fa > {ofn}/tmp.bed',shell=True)
                        afile = pd.read_table(f'{ofn}/tmp.bed')
                        if afile.shape[0]>1:
                            print(f'{seq} has poor specificity in the target area, with {afile.shape[0]} similar regions were detected.')
                            sys.exit()
                        else:
                            tmpdict[n] = [int(afile['start']),int(afile['end'])]
                            n+=1
            
                front,rear = comparelist(tmpdict[1],tmpdict[2])
                rear_bed = f'{seqname}\t{rear[0]}\t{rear[1]}'
                tl = rear[1] - rear[0] + 1
                subprocess.run(f'seqkit subseq -r {front[0]}:{front[1]} {ofn}/barcode/Target.fa -j {nt} -w0 > {ofn}/barcode/Tar_{tl}.fa',shell=True)
                subprocess.run(f'''sed -i 's/{seqname}/front/g' {ofn}/barcode/Tar_{tl}.fa''',shell=True)
                subprocess.run(f'seqkit subseq -r {rear[0]}:{rear[1]} {ofn}/barcode/Target.fa -j {nt} -w0 >> {ofn}/barcode/Tar_{tl}.fa',shell=True)
                subprocess.run(f'''sed -i 's/{seqname}/rear/g' {ofn}/barcode/Tar_{tl}.fa''',shell=True)
                subprocess.run(f'makeblastdb -dbtype nucl -in {ofn}/barcode/Target.fa -input_type fasta -out {ofn}/barcode/Target',shell=True)
                subprocess.run(f'''blastn -db  {ofn}/barcode/Target -query {ofn}/barcode/Tar_{tl}.fa -out {ofn}/barcode/Tar_{tl}_blast -num_threads 10 -evalue 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -task blastn-short -word_size 4''',shell=True)
                blastdb = pd.read_table(f'{ofn}/barcode/Tar_{tl}_blast',names=['sequence','Reference genome name','Similarity(%)','Length','Difference quantity','Gap number','Start position','End position','Reference start','Reference end','Evalue','Score'])
                
    return tl
#-----Identify the identifier sequence---------
def finding_barcode(bfile_out,sn):
    os.chdir(bfile_out)
    bar_dict = {6:['GAGTCT','AGACTC'],8:['GAGTCTAG','CTAGACTC'],10:['GAGTCTAGAT','ATCTAGACTC']}
    afile = pd.read_table(f'{bfile_out}/barcoding_summary.txt')
    barnum = len(afile['barcode_front_refseq'].tolist()[0])
    passrdnum = os.popen(f'cat {bfile_out}/barcoding_summary.txt|wc -l').read().strip()
    barcode = bar_dict.get(barnum)
    frotnum = sum(afile['barcode_front_foundseq']==barcode[0])
    rearnum = sum(afile['barcode_rear_foundseq']==barcode[1])
    bothID = afile[(afile['barcode_front_foundseq']==barcode[0]) & (afile['barcode_rear_foundseq']==barcode[1])]['read_id']
    bothnum = sum((afile['barcode_front_foundseq']==barcode[0]) & (afile['barcode_rear_foundseq']==barcode[1]))
    orID = afile[(afile['barcode_front_foundseq']==barcode[0]) | (afile['barcode_rear_foundseq']==barcode[1])]['read_id']
    ornum = sum((afile['barcode_front_foundseq']==barcode[0]) | (afile['barcode_rear_foundseq']==barcode[1]))
    bothID.to_csv(f'{sn}_B{barnum}_bothID.txt',index=0,header=0)
    orID.to_csv(f'{sn}_B{barnum}_barcoding_ID.txt',index=0,header=0)
    subprocess.run(f'grep -vf {sn}_B{barnum}_bothID.txt {sn}_B{barnum}_barcoding_ID.txt > {sn}_B{barnum}_nobothID.txt ',shell=True)
    nobothnum = os.popen(f'cat {sn}_B{barnum}_nobothID.txt|wc -l ').read().strip()
    print(f'Barcode_len:{barnum}\nbarcode:{barcode}\npassrdnum:{passrdnum}\nfrotnum:{frotnum}\nrearnum:{rearnum}\nbothnum:{bothnum}\nornum:{ornum}\nnobothnum:{nobothnum}\n')
#----basecalling-------
def basecalling_pro(f5f,sn,of):
    print(f'Basecalling on {sn} in progress...')
    with open(f'{wkdir}/task.log','a') as f:
        subprocess.run(f'guppy_basecaller -i {f5f} -s {of}/{sn}_output -x auto -c  dna_r9.4.1_450bps_sup.cfg',stdout=f,stderr=f,shell=True)
        
def barcoding(fqf,qscore,of,sn):
    print(f'Barcoding of {sn} in progress...')
    bk_dict = {6:'MY-CUSTOM-BARCODES',8:'MY-CUSTOM-BARCODES08',10:'MY-CUSTOM-BARCODES10'}
    bk = bk_dict.get(Blen)
    cmd=f'''/home/gqw/STR_soft/ont-guppy/bin/guppy_barcoder \\
    --input_path {fqf} \\
    --save_path  {of}/barcoding_output \\
    --data_path /home/gqw/STR_soft/ont-guppy/my_data/barcoding \\
    --barcode_kits {bk} \\
    --min_score {qscore} \\
    --min_score_rear_override {qscore} \\
    --min_score_mask {qscore} \\
    '''
    with open(f'{wkdir}/task.log','a') as f:
        subprocess.run(cmd,stdout=f,stderr=f,shell=True)
    finding_barcode(f'{of}/barcoding_output',sn)
    
print(dfile)    
def precess_rawdata(rf,ofn,nt,md):
    os.makedirs(f'{ofn}/rawdata')
    os.chdir(f'{ofn}/rawdata')
    if md == 'fq':
        sam = argv.samplename
        os.makedirs(f'{ofn}/rawdata/{sn}_output/pass')
        if os.path.isfile(dfile):
            subprocess.run(f'ln -s {dfile} {ofn}/rawdata/{sn}_output/pass/All.fastq',shell=True)
        elif os.path.isdir(dfile):
            subprocess.run(f'cat {dfile}/* > {ofn}/rawdata/{sn}_output/pass/All.fastq',shell=True)
        custom_barcoding(f'{ofn}/rawdata/{sam}_output/pass/All.fastq',f'{ofn}/barcoding_output',sam,barcodefa,barmis,Ifbar)
    elif md == 'f5':
        sam = argv.samplename
        basecalling_pro(dfile,sam,f'{ofn}/rawdata/')
        barcoding(f'{ofn}/rawdata/{sn}_output/pass',60,ofn,sam)

    else:
        with open(dfile) as f:
            for line in f:
                line = line.strip().split('\t')
                sam = line[0]
                nfile = line[1]
                basecalling_pro(nfile,sam,f'{ofn}/rawdata/')
                barcoding(f'{ofn}/rawdata/{sam}_output/pass',60,ofn,sam)


def check_group(group):
    front_refstart = group[group['Query'] == 'front']['RefEnd'].max()
    rear_refstart = group[group['Query'] == 'rear']['RefStart'].min()
    if rear_refstart < front_refstart:
        return group['Ref'].iloc[0]


def iden_tar(nt,qv,tarl,sam,bl):
    os.makedirs(f'{wkdir}/1.Mapping')
    os.chdir(f'{wkdir}/1.Mapping')
    if not os.path.isfile('All.fastq'):
        subprocess.run(f'ln -s {wkdir}/rawdata/{sn}_output/pass/All.fastq ./All.fastq',shell=True)
    subprocess.run(f'minimap2 -ax map-ont {ref} All.fastq -t {nt}|samtools sort --threads {nt} -o tmp.sort.bam -t 10 ',shell=True)
    with open('mapping.log','w') as f:
        subprocess.run(f'samtools view tmp.sort.bam |awk \'{{if($2==0)print $1}}\' > none_senseID.txt',shell=True,stdout=f,stderr=f)
    subprocess.run(f'samtools view -F 276 tmp.sort.bam |cut -f1 > t_forward_ID.txt',shell=True)
    subprocess.run(f'samtools view -F 260 -f 16 tmp.sort.bam |cut -f1 > t_reverse_ID.txt',shell=True)
    subprocess.run(f'samtools view -f 4 tmp.sort.bam|cut -f1 > no_mapping_ID.txt',shell=True)
    subprocess.run(f'cp t_reverse_ID.txt reverse_ID.txt',shell=True)
    subprocess.run(f'grep -vf reverse_ID.txt t_forward_ID.txt > forward_ID.txt',shell=True)
    subprocess.run(f'seqkit grep -f forward_ID.txt All.fastq -j {nt} > forward.fastq',shell=True)
    subprocess.run(f'seqkit grep -f reverse_ID.txt All.fastq -j {nt} > t_reverse.fastq',shell=True)
    subprocess.run(f'seqkit grep -f no_mapping_ID.txt All.fastq > no_mapping.fastq',shell=True)
    subprocess.run(f'seqkit seq -p -r t_reverse.fastq -j {nt} > reverse.fastq',shell=True)
    subprocess.run(f'cat forward.fastq reverse.fastq no_mapping.fastq > t_new.fastq',shell=True)
    subprocess.run(f'seqkit rmdup t_new.fastq -j {nt} > new.fastq',shell=True)
    with open('task.log','a') as f:
        subprocess.run(f'seqkit fish -f ../barcode/Tar_{tarl}.fa new.fastq -b {sam}.tar{tarl}.bam -q {qv} -j {nt} > result.tsv 2>&1',shell=True,stdout=f,stderr=f)
    fishdb = pd.read_table('result.tsv')
    fgtrIDlist = [i for i in applyParallel1(fishdb.groupby('Ref'),check_group) if i != None]
    fgtrID = '\n'.join(fgtrIDlist)
    open('frontgtrearID.txt','w').write(f'{fgtrID}')
    subprocess.run(f'samtools view {sam}.tar{tarl}.bam -F 16 > {sam}.tar{tarl}.q{qv}.sam',shell=True)
    #1.Remove read with double prefix and suffix sequences 2. Remove reads with double identifier sequences
    subprocess.run(f'''grep rear  {sam}.tar{tarl}.q{qv}.sam|cut -f3|sort|uniq -c |awk '{{if($1>1)print $2}}' >  {sam}.tar{tarl}_{qv}_dbrear_ID.txt''',shell=True)
    subprocess.run(f'''grep front  {sam}.tar{tarl}.q{qv}.sam|cut -f3|sort|uniq -c |awk '{{if($1>1)print $2}}' >  {sam}.tar{tarl}_{qv}_dbfront_ID.txt''',shell=True)
    subprocess.run(f'''cat {wkdir}/barcoding_output/**_barcoding_ID.txt |sort|uniq -c |awk '{{if($1>1)print $2}}' > upperone_barcodeID.txt ''',shell=True)
    subprocess.run(f'cat {sam}.tar{tarl}_{qv}_dbrear_ID.txt {sam}.tar{tarl}_{qv}_dbfront_ID.txt upperone_barcodeID.txt frontgtrearID.txt|sort -u > {sam}.tar{tarl}_{qv}_db_ID.txt',shell=True)
    subprocess.run(f'''cut -f3 {sam}.tar{tarl}.q{qv}.sam |sort|uniq -c |awk '{{if($1==2)print $2}}' > t_{sam}.tar{tarl}_{qv}_ID''',shell=True)
    subprocess.run(f'grep -vf  {sam}.tar{tarl}_{qv}_db_ID.txt  t_{sam}.tar{tarl}_{qv}_ID >  {sam}.tar{tarl}_{qv}_ID',shell=True)

    
    if Ifbar == 'bar':
        with open(barcodefa) as f:
            for line in f:
                if line.startswith('>'):
                    nbarID = line.strip().replace('>','')
                    if os.path.isfile(f'{wkdir}/barcoding_output/{sam}_{nbarID}_barcoding_ID.txt'):
                        subprocess.run(f'grep -f {wkdir}/barcoding_output/{sam}_{nbarID}_barcoding_ID.txt {sam}.tar{tarl}_{qv}_ID > {sam}.tar{tarl}_{qv}_{nbarID}.filt_ID',shell=True)
                        subprocess.run(f'grep -vf {wkdir}/barcoding_output/{sam}_{nbarID}_barcoding_ID.txt {sam}.tar{tarl}_{qv}_ID > {sam}.tar{tarl}_{qv}_{nbarID}.filt_ID_failed',shell=True)
                        subprocess.run(f'grep -f {sam}.tar{tarl}_{qv}_{nbarID}.filt_ID {sam}.tar{tarl}.q{qv}.sam > {sam}.tar{tarl}.CGG.q{qv}_{nbarID}.sam',shell=True)
                        subprocess.run(f'grep -vf {sam}.tar{tarl}_{qv}_{nbarID}.filt_ID {sam}.tar{tarl}.q{qv}.sam > {sam}.tar{tarl}.CGG.q{qv}_{nbarID}_failed.sam',shell=True)
    else:
        subprocess.run(f'grep -f {wkdir}/barcoding_output/{sam}_nobar_barcoding_ID.txt {sam}.tar{tarl}_{qv}_ID > {sam}.tar{tarl}_{qv}_nobar.filt_ID',shell=True)
        subprocess.run(f'grep -vf {wkdir}/barcoding_output/{sam}_nobar_barcoding_ID.txt {sam}.tar{tarl}_{qv}_ID > {sam}.tar{tarl}_{qv}_nobar.filt_ID_failed',shell=True)
        subprocess.run(f'grep -f {sam}.tar{tarl}_{qv}_nobar.filt_ID {sam}.tar{tarl}.q{qv}.sam > {sam}.tar{tarl}.CGG.q{qv}_nobar.sam',shell=True)
        subprocess.run(f'grep -vf {sam}.tar{tarl}_{qv}_nobar.filt_ID {sam}.tar{tarl}.q{qv}.sam > {sam}.tar{tarl}.CGG.q{qv}_nobar_failed.sam',shell=True)


def cal_STR(samfile,bcID,md,fqfile,cv,nt):
    if not os.path.isdir(f'{wkdir}/2.STR_result'):
        os.makedirs(f'{wkdir}/2.STR_result')
    os.chdir(f'{wkdir}/2.STR_result')
    subprocess.run(f'cp {wkdir}/1.Mapping/tmp.sort.bam ./',shell=True)
    subprocess.run(f'bamToBed -i tmp.sort.bam |awk -v OFS=\'\t\' \'{{print $1"_"$4,$2,$3}}\'  |sort -k1,1V -k2,2n -k3,3n > bam.bed',shell=True)
    subprocess.run(f'bedtools merge -i bam.bed |awk -v OFS=\'\t\' -F \'[_\t]\' \'{{print $1,$3,$4,$2}}\' > bam.merge.bed',shell=True)
    subprocess.run(f'bedtools coverage -a bam.merge.bed -b {tarbed} > Read_genome_cov.bed',shell=True)
    subprocess.run(f'awk -v OFS=\'\t\' \'{{if($6/310>={cv})print $0,$6/310}}\' Read_genome_cov.bed > Read_genome_high_cov.bed',shell=True)
    subprocess.run(f'cut -f4 Read_genome_high_cov.bed|sort -u > high_cov_ID',shell=True)
    len_dict={}
    lenmode = len(md)
    with open(samfile) as f:
        for line in f:
            line = line.strip().split('\t')
            SeqN = line[0]
            if int(line[1]) > 16:
                if int(line[1])-256==16:
                    SeqStrand = "-"
                else:
                    SeqStrand = "+"
            else:
                if int(line[1]) == 16:
                        SeqStrand = "-"
                else:
                        SeqStrand = "+"
            SeqID = line[2]
            SeqPos = line[3]
            SeqQ = line[5]
            if SeqID not in len_dict.keys():
                len_dict[SeqID] = {}
                len_dict[SeqID][SeqN] = [SeqStrand,SeqPos,SeqQ]
            else:
                len_dict[SeqID][SeqN] = [SeqStrand,SeqPos,SeqQ]

    for keys in len_dict.keys():
        try:
            if len_dict[keys]['front'][0]=='-' and len_dict[keys]['rear'][0]=='-':
                SeqQ = len_dict[keys]['rear'][2]
                strlen = abs(int(len_dict[keys]['front'][1])-int(len_dict[keys]['rear'][1]))-sum([int(i) for i in re.findall('(\d+)[M|D]',SeqQ)]) + shift_num
                STR_count = f'{strlen/lenmode:.1f}'
                newline = f'{keys}\t{strlen}\t{STR_count}\t{qvalue}\t{tarlen}\t{cov}\t{sn}\n'
                startP = int(len_dict[keys]['front'][1])
                endP = int(len_dict[keys]['rear'][1])
                newline1 = f'{keys}\t{endP}\t{startP}\t-\t{strlen}\t{STR_count}\t{qvalue}\t{tarlen}\t{cov}\t{sn}\n'
                if tarlen > 0:
                    open(f'STR_summary_{bcID}.tsv','a').write(newline)
                    open(f'STR_summary1_{bcID}.tsv','a').write(newline1)
            elif len_dict[keys]['front'][0]=='+' and len_dict[keys]['rear'][0]=='+':
                SeqQ = len_dict[keys]['front'][2]
                strlen = abs(int(len_dict[keys]['front'][1])-int(len_dict[keys]['rear'][1]))-sum([int(i) for i in re.findall('(\d+)[M|D]',SeqQ)])+shift_num
                STR_count = f'{strlen/lenmode:.1f}'
                newline = f'{keys}\t{strlen}\t{STR_count}\t{qvalue}\t{tarlen}\t{cov}\t{sn}\n'
                startP = int(len_dict[keys]['front'][1])+sum([int(i) for i in re.findall('(\d+)[M|D]',SeqQ)])
                endP = int(len_dict[keys]['rear'][1])
                newline2 = f'{keys}\t{startP}\t{endP}\t+\t{strlen}\t{STR_count}\t{qvalue}\t{tarlen}\t{cov}\t{sn}\n'
                if tarlen > 0:
                    open(f'STR_summary_{bcID}.tsv','a').write(newline)
                    open(f'STR_summary1_{bcID}.tsv','a').write(newline2)

            else:
                pass
        except:
            open(f'STR_failedID_{bcID}.tsv','a').write(f'{keys}\n')
    if os.path.isfile(f'STR_summary_{bcID}.tsv'):
        subprocess.run(f'seqkit fq2fa {wkdir}/1.Mapping/new.fastq > new.fa',shell=True)
        subprocess.run(f'bedtools getfasta -bed STR_summary1_{bcID}.tsv -fi new.fa|seqkit fx2tab --gc > {bcID}_gcSum.tsv',shell=True)
        gcdb = pd.read_table(f'{bcID}_gcSum.tsv',header=None,names=['TarSeq','tmp1','gc'])
        gcdb.insert(0,'seqID',gcdb.index.str.split(':').str[0])
        gcdb = gcdb[['seqID','TarSeq','gc']]
        gcdb.to_csv(f'{bcID}_gc.tsv',sep='\t',index=False)
        gcdb.loc[gcdb['gc']>gc_content,'seqID'].to_csv(f'{bcID}_highgcID.txt',index=False,header=False)
        subprocess.run(f'''awk -v OFS='\t' '{{if($5>0)print $0}}' STR_summary_{bcID}.tsv > STR_summary_{bcID}_filt.tsv''',shell=True)
        subprocess.run(f'''awk -v OFS='\t' '{{if($5<=0)print $0}}' STR_summary_{bcID}.tsv > STR_summary_{bcID}_filt_failed.tsv''',shell=True)
        subprocess.run(f'grep -f high_cov_ID STR_summary_{bcID}_filt.tsv > t1_STR_summary_hcov_{bcID}.tsv',shell=True)
        subprocess.run(f'grep -f {bcID}_highgcID.txt t1_STR_summary_hcov_{bcID}.tsv > t_STR_summary_hcov_{bcID}.tsv',shell=True)
        subprocess.run(f'grep -vf {bcID}_highgcID.txt t1_STR_summary_hcov_{bcID}.tsv > STR_summary_lowgc_{bcID}.tsv',shell=True)
        subprocess.run(f'grep -vf high_cov_ID STR_summary_{bcID}_filt.tsv > STR_summary_hcov_{bcID}_failed.tsv',shell=True)
        #---Remove the read with identifier sequence between the prefix and suffix seuqences
        ttbardb = pd.read_table(f'{wkdir}/barcoding_output/{bcID}_seqgrepleng.tsv')
        ttstrdb = pd.read_table(f'{wkdir}/2.STR_result/STR_summary1_{bcID}.tsv',header=None,names=['seqID','strs','stre','strstand','strlen','strnum','qvalue','tarlen','cov','Sname'])
        nnstrdb = ttstrdb.merge(ttbardb,on='seqID')
        n1strdb = nnstrdb[(nnstrdb['start'] > nnstrdb['strs']) & (nnstrdb['start'] < nnstrdb['stre'])]
        n1strdb['seqID'].to_csv(f'STR_{bcID}_minbar_ID.txt',index=False,header=False)
        subprocess.run(f'grep -vf STR_{bcID}_minbar_ID.txt t_STR_summary_hcov_{bcID}.tsv > STR_summary_hcov_{bcID}.tsv',shell=True)
        if int(os.popen(f'cat STR_summary_hcov_{bcID}.tsv|wc -l').read()) > 0:
            print(f'Detecting the CGG repeats of {sn}_{bcID}...')
            subprocess.run(f'Rscript /home/gqw/STR_soft/STR_test/script/SummaryplotSTR.R -i {wkdir}/2.STR_result/STR_summary_hcov_{bcID}.tsv -n {sn}_{bcID} -l {lowper} --SE {SE}',shell=True)
            print(f'Finished the detecion of CGG repeats.')
            if STRmotif != 'nomotif':
                subprocess.run(f'seqkit locate -p {STRmotif} {wkdir}/1.Mapping/new.fastq -P > {STRmotif}.tsv',shell=True)
                subprocess.run(f'''sed '1d'  {STRmotif}.tsv |awk -v OFS='\t' '{{print $1,$5,$6}}' > tmp_{STRmotif}.tsv ''',shell=True)
                subprocess.run(f''' awk -v OFS='\t' '{{if($2<$3)print $0}}' ../2.STR_result/STR_summary1_{bcID}.tsv > tmp_summary1_STR_{bcID}.tsv ''',shell=True)
                subprocess.run(f'''bedtools intersect -a tmp_{STRmotif}.tsv -b tmp_summary1_STR_{bcID}.tsv -wa -wb > tmp_{STRmotif}_{bcID}_over.tsv ''',shell=True)
                tmpile = pd.read_table(f'tmp_{STRmotif}_{bcID}_over.tsv',names=['r1','s1','e1','r2','s2','e2','strand','strlen','strnum','qvalue','tarlen','cov','sam'])
                tmpile[f'{STRmotif}_pos'] = ((tmpile['s1'] - tmpile['s2'])/3)+1
                tmpile.to_csv(f'STR_sum_{STRmotif}_{bcID}.tsv',sep='\t',index=False)
                print(f'Detecting the AGG interruption of {sn}_{bcID}...')
                subprocess.run(f'Rscript /home/gqw/STR_soft/STR_test/script/SummaryplotSTR_motif.R -i {wkdir}/2.STR_result/STR_sum_{STRmotif}_{bcID}.tsv -n {sn}_{bcID} -l {lowper1} --motif {STRmotif} --Mut motif --input2 {wkdir}/2.STR_result/STR_summary_hcov_{bcID}.tsv --motif_raw {mode}',shell=True)
                print(f'Finished the detection of AGG interruption.')
            if Mutmotif != 'noMut':
                subprocess.run(f'seqkit locate -p {Mutmotif} {wkdir}/1.Mapping/new.fastq -P > {Mutmotif}.tsv',shell=True)
                subprocess.run(f'''sed '1d'  {Mutmotif}.tsv |awk -v OFS='\t' '{{print $1,$5,$6}}' > tmp_{Mutmotif}.tsv ''',shell=True)
                subprocess.run(f'''bedtools intersect -a tmp_{Mutmotif}.tsv -b tmp_summary1_STR_{bcID}.tsv -wa -wb > tmp_{Mutmotif}_{bcID}_over.tsv ''',shell=True)
                tmpile2 = pd.read_table(f'tmp_{Mutmotif}_{bcID}_over.tsv',names=['r1','s1','e1','r2','s2','e2','strand','strlen','strnum','qvalue','tarlen','cov','sam'])
                tmpile2[f'{Mutmotif}_pos'] = (tmpile2['s1'] - tmpile2['s2'])+1
                tmpile2.to_csv(f'STR_sum_{Mutmotif}_{bcID}.tsv',sep='\t',index=False)
                print(f'Detecting the {sn}_{bcID}_{Mutmotif}')
                subprocess.run(f'Rscript /home/gqw/STR_soft/STR_test/script/SummaryplotSTR_motif_Mut.R -i {wkdir}/2.STR_result/STR_sum_{Mutmotif}_{bcID}.tsv -n {sn}_{bcID} -l {lowper} --motif {Mutmotif} --Mut Mut',shell=True)
                print(f'Finished the detection of {sn}_{bcID}_{Mutmotif}.')
        else:
            print(f'The data of {bcID} is insufficient, and the detection of STR has failed.')
    else:
        print(f'The data of {bcID} is insufficient, and the detection of STR has failed.')

#-----Main process----
tarlen = tar_select(tarfile,ref,ofn,tarlen,threads)
precess_rawdata(ref,wkdir,threads,fmode)
iden_tar(threads,qvalue,tarlen,sn,Blen)
if Ifbar == 'bar':
    with open(barcodefa) as f:
        for line in f:
            if line.startswith('>'):
                nbarID = line.strip().replace('>','')
                os.makedirs(f'{wkdir}/3.Sum_result/{nbarID}')
                if os.path.isfile(f'{wkdir}/1.Mapping/{sn}.tar{tarlen}.{mode}.q{qvalue}_{nbarID}.sam'):
                    cal_STR(f'{wkdir}/1.Mapping/{sn}.tar{tarlen}.{mode}.q{qvalue}_{nbarID}.sam',nbarID,mode,f'{wkdir}/1.Mapping/All.fastq',cov,threads)
                    subprocess.run(f'mv {wkdir}/1.Mapping/*{nbarID}* {wkdir}/3.Sum_result/{nbarID} 2>/dev/null',shell=True)
                    subprocess.run(f'mv {wkdir}/2.STR_result/*{nbarID}* {wkdir}/3.Sum_result/{nbarID} 2>/dev/null',shell=True)
                    subprocess.run(f'mv {wkdir}/barcoding_output/*{nbarID}* {wkdir}/3.Sum_result/{nbarID} 2>/dev/null',shell=True)
                    #--- Summary
                    #1.STR 2.AGG
                    if os.path.isfile(f'{wkdir}/3.Sum_result/{nbarID}/{nbarID}_seqgrepleng.tsv') and os.path.isfile(f'{wkdir}/3.Sum_result/{nbarID}/STR_summary1_{nbarID}.tsv'):
                        nrdb = pd.read_table(f'{wkdir}/3.Sum_result/{nbarID}/{nbarID}_seqgrepleng.tsv')
                        tstrdb = pd.read_table(f'{wkdir}/3.Sum_result/{nbarID}/STR_summary1_{nbarID}.tsv',header=None,names=['seqID','strs','stre','strstand','strlen','strnum','qvalue','tarlen','cov','Sname'])
                        gcdb = pd.read_table(f'{wkdir}/3.Sum_result/{nbarID}/{nbarID}_gc.tsv')
                        nstrdb = tstrdb.merge(nrdb,on='seqID')
                        nstrdb =  tstrdb.merge(gcdb,on='seqID')
                        nstrdb.to_csv(f'{wkdir}/3.Sum_result/{nbarID}/{nbarID}_STR_summary_new.tsv',index=False,sep='\t')
                        if os.path.isfile(f'{wkdir}/3.Sum_result/{nbarID}/STR_sum_{STRmotif}_{nbarID}.tsv'):
                            tagdb = pd.read_table(f'{wkdir}/3.Sum_result/{nbarID}/STR_sum_{STRmotif}_{nbarID}.tsv')
                            ntagdb = tagdb.merge(nrdb,left_on='r1',right_on='seqID')
                            ntagdb.to_csv(f'{wkdir}/3.Sum_result/{nbarID}/{nbarID}_{STRmotif}_summary_new.tsv',index=False,sep='\t')

    subprocess.run(f'iconv -f UTF-8 -t GBK {wkdir}/2.STR_result/STR_statis.tsv > {wkdir}/2.STR_result/STR_statis_EN.tsv',shell=True)
    subprocess.run(f'iconv -f UTF-8 -t GBK {wkdir}/2.STR_result/STR_statis_new.tsv > {wkdir}/2.STR_result/STR_statis_new_EN.tsv',shell=True)
    subprocess.run(f'cp {wkdir}/2.STR_result/STR_statis_new.tsv  {wkdir}/3.Sum_result/STR_statis_new.tsv',shell=True)
    subprocess.run(f'cp {wkdir}/2.STR_result/STR_statis_new_EN.tsv  {wkdir}/3.Sum_result/STR_statis_new_EN.tsv',shell=True)


else:
    nbarID = 'nobar'
    cal_STR(f'{wkdir}/1.Mapping/{sn}.tar{tarlen}.{mode}.q{qvalue}_{nbarID}.sam',nbarID,mode,f'{wkdir}/1.Mapping/All.fastq',cov,threads)
