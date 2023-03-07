[TOC]


### 1.1.1 环境变量设置(每次开始分析前必须运行)

设置数据库、软件和工作目录

    # 公共数据库database(db)位置，如管理员设置/db，个人下载至~/db
    db=/db
    # Conda软件software安装目录，`conda env list`命令查看，如~/miniconda2
    soft=/conda2
    # 设置工作目录work directory(wd)，如meta
    wd=~/meta
    # 创建并进入工作目录
    mkdir -p $wd

    pwd
    
    cd $wd
    
    pwd
    
    # 方法1.加载自己安装环境
    # conda activate metagenome_env
    
    # 方法2.加载别人安装环境
    export PATH=${soft}/envs/metagenome_env/bin/:$PATH
    source ${soft}/bin/activate metagenome_env
    
    # 指定某个R语言环境
    alias Rscript="/anaconda2/bin/Rscript --vanilla"

### 1.1.2 起始文件——序列和元数据

创建3个常用子目录：序列，临时文件和结果
 
    mkdir -p seq temp result
    /bin/rm -f result/metadata.txt
    # 上传元数据metadata.txt至result目录，此处下载并重命名
    wget http://210.75.224.110/github/EasyMetagenome/result/metadata2.txt -O result/metadata.txt
    # 检查文件格式，^I为制表符，$为Linux换行，^M$为Windows回车，^M为Mac换行符
    cat -A result/metadata.txt
    # 转换Windows回车为Linux换行
    sed -i 's/\r//' result/metadata.txt
    cat -A result/metadata.txt

用户使用filezilla上传测序文件至seq目录，本次从其它位置复制，或从网络下载测试数据
    
    # 方法1. 从其它目录复制测序数据
    cp -rf /db/meta/seq/*.gz seq/
    
    # 方法2. 网络下载测试数据
    cd seq/
    awk '{system("wget -c http://210.75.224.110/github/EasyMetagenome/seq/"$1"_1.fq.gz")}' \
      <(tail -n+2 ../result/metadata.txt)
    awk '{system("wget -c http://210.75.224.110/github/EasyMetagenome/seq/"$1"_2.fq.gz")}' \
      <(tail -n+2 ../result/metadata.txt)
    cd ..

    # 查看文件大小
    ls -lsh seq
    # -l 列出详细信息 (l: list)
    # -sh 显示人类可读方式文件大小 (s: size; h: human readable)	

### 1.1.3 了解工作目录和文件

显示文件结构

    tree -L 2
    # .
    # ├── pipeline.sh
    # ├── result
    # │   └── metadata.txt
    # ├── seq
    # │   ├── C1_1.fq.gz
    # │   ├── C1_2.fq.gz
    # │   ├── N1_1.fq.gz
    # │   └── N1_2.fq.gz
    # └── temp

- pipeline.sh是分析流程代码；
- seq目录中有2个样本Illumina双端测序，4个序列文件；
- temp是临时文件夹，存储分析中间文件，结束可全部删除节约空间
- result是重要节点文件和整理化的分析结果图表，
- 实验设计metadata.txt也在此


## 1.2 (可选)FastQC质量评估

    # 第一次使用软件要记录软件版本，文章方法中必须写清楚
    fastqc --version # 0.11.8
    # time统计运行时间，fastqc质量评估
    # *.gz为原始数据，-t指定多线程
    time fastqc seq/*.gz -t 2

质控报告见`seq`目录，详细解读请阅读[《数据的质量控制软件——FastQC》](https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)。
    
multiqc将fastqc的多个报告生成单个整合报告，方法批量查看和比较

    # 记录软件版本
    multiqc --version # 1.5
    # 整理seq目录下fastqc报告，输出multiqc_report.html至result/qc目录
    multiqc -d seq/ -o result/qc

查看右侧result/qc目录中multiqc_report.html，单击，选择`View in Web Browser`查看可交互式报告。

## 1.3 KneadData质控和去宿主

kneaddata是流程，它主要依赖trimmomatic质控和去接头，bowtie2比对宿主，然后筛选非宿主序列用于下游分析 。

详细教程和常见问题，阅读：[MPB：随机宏基因组测序数据质量控制和去宿主的分析流程和常见问题](https://mp.weixin.qq.com/s/ovL4TwalqZvwx5qWb5fsYA)

    # 记录核心软件版本
    kneaddata --version # 0.6.1
    trimmomatic -version # 0.39
    bowtie2 --version # 2.3.5
    # 可只选一行中部分代码点击Run，如选中下行中#号后面命令查看程序帮助
    # kneaddata -h # 显示帮助
    
检查点：zless/zcat查看可压缩文件，检查序列质量格式(质量值大写字母为标准Phred33格式，小写字母为Phred64，需参考附录：质量值转换)；检查双端序列ID是否重复，如果重名需要在质控前改名更正。参考**附录，质控kneaddata，去宿主后双端不匹配——序列改名**。

    # 设置某个样本名为变量i，以后再无需修改
    i=C1
    # zless查看压缩文件，空格翻页，按q退出。
    zless seq/${i}_1.fq.gz | head -n4
    # zcat显示压缩文件，head指定显示行数
    zcat seq/${i}_2.fq.gz | head -n4

- "|" 为管道符，上一个命令的输出，传递给下一个命令做输入
- gzip: stdout: Broken pipe：管道断开。这里是人为断开，不是错误
- 运行过程中需要仔细阅读屏幕输出的信息

如果序列双端名称一致，改名参见附录代码：质控KneadData - 去宿主后双端不匹配——序列改名



### 1.3.2 多样品并行质控

并行队列管理软件——“parallel”。
    
    # 记录软件版本
    parallel --version # 20160222
    # 打will cite承诺引用并行软件parallel
    parallel --citation 
    
parallel软件说明和使用实例

    # 根据样本列表`:::`并行处理，并行j=2个任务，每个任务t=3个线程，2~7m
    # 运行下面这行，体会下parallel的工作原理
    # ::: 表示传递参数；第一个::: 后面为第一组参数，对应于{1};
    # 第二个::: 后面为第二组参数，对应于{2}，依次替换
    parallel -j 3 --xapply "echo {1} {2}" ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz
    # --xapply保持文件成对，否则将为两两组合，显示如下：
    parallel -j 2 "echo {1} {2}" ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz
    # 从文件列表使用
    parallel -j 3 --xapply "echo seq/{1}_1.fq.gz seq/{1}_2.fq.gz" ::: `tail -n+2 result/metadata.txt|cut -f1`
    

质控结果汇总

    # 采用kneaddata附属工具kneaddata_read_count_table
    kneaddata_read_count_table --input temp/qc \
      --output temp/kneaddata.txt
    # 筛选重点结果列
    cut -f 1,2,4,12,13 temp/kneaddata.txt | sed 's/_1_kneaddata//' > result/qc/sum.txt
    cat result/qc/sum.txt

    # 用R代码统计下质控结果
    Rscript -e "data=read.table('result/qc/sum.txt', header=T, row.names=1, sep='\t'); summary(data)"
    # R转换宽表格为长表格
    Rscript -e "library(reshape2); data=read.table('result/qc/sum.txt', header=T,row.names=1, sep='\t'); write.table(melt(data), file='result/qc/sum_long.txt',sep='\t', quote=F, col.names=T, row.names=F)"
    cat result/qc/sum_long.txt
    # 可用 http://www.ehbio.com/ImageGP/ 绘图展示



# 二、基于读长分析 Read-based (HUMAnN2)

## 2.1 准备HUMAnN2输入文件

小技巧：循环批量处理样本列表

    # 基于样本元数据提取样本列表命令解析
    # 去掉表头
    tail -n+2 result/metadata.txt
    # 提取第一列样本名
    tail -n+2 result/metadata.txt|cut -f1
    # 循环处理样本
    for i in `tail -n+2 result/metadata.txt|cut -f1`;do echo "Processing "$i; done
    # ` 反引号为键盘左上角Esc键下面的按键，一般在数字1的左边，代表运行命令返回结果

HUMAnN2要求双端序列合并的文件作为输入，for循环根据实验设计样本名批量双端序列合并。
注意星号和问号，分别代表多个和单个字符。当然大家更不能溜号~~~行分割的代码行末有一个 \ 
        
    mkdir -p temp/concat
    # 双端合并为单个文件
    for i in `tail -n+2 result/metadata.txt|cut -f1`;do 
      cat temp/qc/${i}_1_kneaddata_paired_?.fastq \
      > temp/concat/${i}.fq; 
    done
    # 查看样品数量和大小
    ls -sh temp/concat/*.fq
    # 数据太大，计算时间长，可用head对单端分析截取20M序列，即3G，则为80M行，详见附录：HUMAnN2减少输出文件加速



在通路丰度中添加分组

    ## 提取样品列表
    head -n1 result/humann2/pathabundance.tsv | sed 's/# Pathway/SampleID/' | tr '\t' '\n' > temp/header
    ## 对应分组
    awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' result/metadata.txt temp/header | tr '\n' '\t'|sed 's/\t$/\n/' > temp/group
    # 合成样本、分组+数据
    cat <(head -n1 result/humann2/pathabundance.tsv) temp/group <(tail -n+2 result/humann2/pathabundance.tsv) \
      > result/humann2/pathabundance.pcl
    head -n5 result/humann2/pathabundance.pcl


## 2.6 LEfSe差异分析物种


前面演示数据仅有2个样本，无法进行差异比较。下面使用result12目录中由12个样本生成的结果表进行演示

    # 设置结果目录，自己的数据使用result，此处演示使用result12
    result=result12

准备输入文件，修改样本品为组名(可手动修改)

    # 预览输出数据
    head -n3 $result/metaphlan2/taxonomy.tsv
    # 提取样本行，替换为每个样本一行，修改ID为SampleID
    head -n1 $result/metaphlan2/taxonomy.tsv|tr '\t' '\n'|sed '1 s/ID/SampleID/' >temp/sampleid
    head -n3 temp/sampleid
    # 提取SampleID对应的分组Group(假设为metadata.txt中第二列$2)，替换换行\n为制表符\t，再把行末制表符\t替换回换行
    awk 'BEGIN{OFS=FS="\t"}NR==FNR{a[$1]=$2}NR>FNR{print a[$1]}' $result/metadata.txt temp/sampleid|tr '\n' '\t'|sed 's/\t$/\n/' >groupid
    cat groupid
    # 合并分组和数据(替换表头)
    cat groupid <(tail -n+2 $result/metaphlan2/taxonomy.tsv) > $result/metaphlan2/lefse.txt
    head -n3 $result/metaphlan2/lefse.txt

    
    # 门水平去除脊索动物
    grep 'Chordata' result/kraken2/bracken.P.0.01
    grep -v 'Chordata' result/kraken2/bracken.P.0.01 > result/kraken2/bracken.P.0.01-H

    # 按物种名手动去除宿主污染，以人为例(需按种水平计算相关结果)
    # 种水平去除人类P:Chordata,S:Homo sapiens
    grep 'Homo sapiens' result/kraken2/bracken.S.0.01
    grep -v 'Homo sapiens' result/kraken2/bracken.S.0.01 \
      > result/kraken2/bracken.S.0.01-H
    

### metaphlan2-共有或特有物种网络图

    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=9;i<=NF;i++) a[i]=$i; print "Tax\tGroup"} \
       else {for(i=9;i<=NF;i++) if($i>0.05) print "Tax_"FNR, a[i];}}' \
       result/metaphlan2/taxonomy.spf > result/metaphlan2/taxonomy_highabundance.tsv
       
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {print "Tax\tGrpcombine";} else a[$1]=a[$1]==""?$2:a[$1]$2;}END{for(i in a) print i,a[i]}' \
       result/metaphlan2/taxonomy_highabundance.tsv > result/metaphlan2/taxonomy_group.tsv
    
    cut -f 2 result/metaphlan2/taxonomy_group.tsv | tail -n +2 | sort -u >group
    
    for i in `cat group`; do printf "#%02x%02x%02x\n" $((RANDOM%256)) $((RANDOM%256)) $((RANDOM%256)); done >colorcode
    
    paste group colorcode >group_colorcode
    
    awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2;}ARGIND==2{if(FNR==1) {print $0, "Grpcombinecolor"} else print $0,a[$2]}' \
       group_colorcode result/metaphlan2/taxonomy_group.tsv > result/metaphlan2/taxonomy_group2.tsv
    
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {print "Tax",$1,$2,$3,$4, $5, $6, $7, $8 } else print "Tax_"FNR, $1,$2,$3,$4, $5,$6, $7, $8}' \
       result/metaphlan2/taxonomy.spf > result/metaphlan2/taxonomy_anno.tsv


