#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chorus2 - 寡核苷酸FISH探针设计软件
用于设计特异性寡核苷酸探针，支持基因组大小分析和多线程处理

作者: Chorus2开发团队
版本: 2.1
"""

import argparse
import sys
import os
from multiprocessing import Pool
from math import log
from pyfasta import Fasta

# 导入自定义模块
from Choruslib import bwa
from Choruslib import jellyfish
from Choruslib import prefilter, primer3_filter, spgenome


def main():
    """主函数：执行Chorus2探针设计流程"""
    # 解析命令行参数并验证
    args = check_options(get_options())
    
    # 计算基因组大小和k-mer长度
    genomesize = int(os.path.getsize(args.genome) / 1e6)
    kmer = max(int(log(genomesize, 4) + 1), 17)
    
    # Jellyfish参数设置
    lowercount = 2
    jfsize = '100M'
    
    # 序列分割参数（大于10M的序列将被分割）
    spsize = 10000000
    
    # 计算最大k-mer得分
    maxkmerscore = int(((args.length * args.homology / 100) - kmer) * args.ploidy / 2 + 0.5)
    
    # 创建多进程池
    jfpool = Pool(args.threads)
    
    # 构建k-mer索引文件路径
    jfkmerfile = os.path.join(args.saved, f"{os.path.basename(args.genome)}_{kmer}mer.jf")
    
    # 检查是否需要重建k-mer索引
    kmerbuild = should_rebuild_index(jfkmerfile, "k-mer索引", args.docker)
    
    # 检查是否需要重建BWA索引
    bwaindexfile = os.path.basename(args.genome)
    bwatestindex = os.path.join(args.saved, bwaindexfile + '.sa')
    bwaindex = os.path.join(args.saved, bwaindexfile)
    bwabuild = should_rebuild_index(bwatestindex, "BWA索引", args.docker)
    
    print(f"基因组大小: {genomesize}MB, k-mer长度: {kmer}, k-mer文件: {jfkmerfile}")
    print(f"k-mer构建: {kmerbuild}, BWA构建: {bwabuild}, 线程数: {args.threads}")
    
    # 构建Jellyfish索引
    if kmerbuild:
        if build_jellyfish_index(args, jfkmerfile, kmer, lowercount, jfsize):
            print("JellyFish计数完成...")
        else:
            print("JellyFish计数错误!!!")
            sys.exit(1)
    else:
        print(f"使用现有文件: {jfkmerfile}")
    
    # 构建BWA索引
    if bwabuild:
        bwa.bwaindex(args.bwa, args.genome, args.saved)
        print("BWA索引构建完成...")
    else:
        print(f"使用现有索引: {bwatestindex}")
    
    # 执行探针过滤
    jffilteredprobe = execute_probe_filtering(
        args, jfpool, jfkmerfile, kmer, maxkmerscore, 
        spsize, genomesize, lowercount, jfsize
    )
    
    # 关闭进程池
    jfpool.close()
    
    # 执行BWA过滤
    bwafiltedpb = execute_bwa_filtering(args, bwaindex, jffilteredprobe)
    
    # 生成输出文件
    generate_output_files(args, bwafiltedpb, bwaindex)
    
    print("任务完成!!")


def should_rebuild_index(index_file, index_type, docker_mode):
    """
    询问用户是否重建索引文件
    
    Args:
        index_file: 索引文件路径
        index_type: 索引类型描述
        docker_mode: 是否在Docker模式下运行
    
    Returns:
        bool: 是否需要重建索引
    """
    if os.path.isfile(index_file) and not docker_mode:
        print(f"发现: {index_file}")
        message = f"发现{index_type}文件 {index_file}。是否重建？按Y或N继续:"
        print(message)
        
        while True:
            char = getch()
            if char.lower() in ("y", "n"):
                print(char)
                return char.lower() == 'y'
    
    return True


def build_jellyfish_index(args, jfkmerfile, kmer, lowercount, jfsize):
    """
    构建Jellyfish索引
    
    Args:
        args: 命令行参数
        jfkmerfile: k-mer文件路径
        kmer: k-mer长度
        lowercount: 最小计数
        jfsize: Jellyfish大小
    
    Returns:
        bool: 构建是否成功
    """
    jfcount = jellyfish.jfcount(
        jfpath=args.jellyfish, 
        mer=kmer, 
        infile=args.genome, 
        output=jfkmerfile,
        threads=args.threads, 
        lowercount=lowercount, 
        size=jfsize
    )
    return jfcount


def execute_probe_filtering(args, jfpool, jfkmerfile, kmer, maxkmerscore, spsize, genomesize, lowercount, jfsize):
    """
    执行探针过滤流程
    
    Args:
        args: 命令行参数
        jfpool: 多进程池
        jfkmerfile: k-mer文件路径
        kmer: k-mer长度
        maxkmerscore: 最大k-mer得分
        spsize: 序列分割大小
        genomesize: 基因组大小
        lowercount: 最小计数
        jfsize: Jellyfish大小
    
    Returns:
        list: 过滤后的探针列表
    """
    jffilteredprobe = []
    
    if genomesize < 1000:
        # 小基因组处理
        jffilteredprobe = process_small_genome(
            args, jfpool, jfkmerfile, kmer, maxkmerscore, spsize
        )
    else:
        # 大基因组处理（>1GB）
        print("基因组大小 > 1GB")
        jffilteredprobe = process_large_genome(
            args, jfpool, jfkmerfile, kmer, maxkmerscore, spsize
        )
    
    return jffilteredprobe


def process_small_genome(args, jfpool, jfkmerfile, kmer, maxkmerscore, spsize):
    """
    处理小基因组（<1GB）
    
    Args:
        args: 命令行参数
        jfpool: 多进程池
        jfkmerfile: k-mer文件路径
        kmer: k-mer长度
        maxkmerscore: 最大k-mer得分
        spsize: 序列分割大小
    
    Returns:
        list: 过滤后的探针列表
    """
    fastain = Fasta(args.input)
    jffpbrunerlist = create_probe_runners(
        fastain, args, jfkmerfile, kmer, maxkmerscore, spsize
    )
    
    return execute_kmer_filtering(jfpool, jffpbrunerlist, "Jellyfish过滤")


def process_large_genome(args, jfpool, jfkmerfile, kmer, maxkmerscore, spsize):
    """
    处理大基因组（>1GB）
    
    Args:
        args: 命令行参数
        jfpool: 多进程池
        jfkmerfile: k-mer文件路径
        kmer: k-mer长度
        maxkmerscore: 最大k-mer得分
        spsize: 序列分割大小
    
    Returns:
        list: 过滤后的探针列表
    """
    jffilteredprobe = []
    
    # 分割基因组文件
    subFas = spgenome.spgenome(args.input, args.saved)
    
    for subFafile in subFas:
        print(subFafile)
        fastain = Fasta(subFafile)
        
        jffpbrunerlist = create_probe_runners(
            fastain, args, jfkmerfile, kmer, maxkmerscore, spsize
        )
        
        # 执行k-mer过滤
        for curpblist in jfpool.imap_unordered(jellyfish.kmerfilterprobe, jffpbrunerlist):
            jffilteredprobe.extend(curpblist)
    
    return jffilteredprobe


def create_probe_runners(fastain, args, jfkmerfile, kmer, maxkmerscore, spsize):
    """
    创建探针运行器列表
    
    Args:
        fastain: 输入序列
        args: 命令行参数
        jfkmerfile: k-mer文件路径
        kmer: k-mer长度
        maxkmerscore: 最大k-mer得分
        spsize: 序列分割大小
    
    Returns:
        list: 探针运行器列表
    """
    jffpbrunerlist = []
    
    for seqname in fastain.keys():
        chrlen = len(fastain[seqname])
        
        if chrlen < spsize:
            # 短序列直接处理
            jffpbruner = jellyfish.JFfpbruner(
                jfpath=args.jellyfish, 
                jfkmerfile=jfkmerfile, 
                mer=kmer,
                pyfasta=fastain, 
                seqname=seqname, 
                pblength=args.length,
                maxkmerscore=maxkmerscore, 
                start=0, 
                end=chrlen - 1, 
                step=args.step
            )
            jffpbrunerlist.append(jffpbruner)
        else:
            # 长序列分割处理
            chrblock = int(chrlen / spsize) + 1
            for i in range(chrblock):
                start = i * spsize
                end = min(start + spsize - 1, chrlen - 1)
                
                jffpbruner = jellyfish.JFfpbruner(
                    jfpath=args.jellyfish, 
                    jfkmerfile=jfkmerfile, 
                    mer=kmer,
                    pyfasta=fastain, 
                    seqname=seqname, 
                    pblength=args.length,
                    maxkmerscore=maxkmerscore, 
                    start=start, 
                    end=end, 
                    step=args.step
                )
                jffpbrunerlist.append(jffpbruner)
    
    return jffpbrunerlist


def execute_kmer_filtering(jfpool, jffpbrunerlist, process_name):
    """
    执行k-mer过滤
    
    Args:
        jfpool: 多进程池
        jffpbrunerlist: 探针运行器列表
        process_name: 进程名称
    
    Returns:
        list: 过滤后的探针列表
    """
    jffilteredprobe = []
    jffinished = 0
    total_count = len(jffpbrunerlist)
    
    print(f"总任务数: {total_count}")
    
    for curpblist in jfpool.imap_unordered(jellyfish.kmerfilterprobe, jffpbrunerlist):
        jffilteredprobe.extend(curpblist)
        jffinished += 1
        print(f"{process_name}: {jffinished}/{total_count}")
    
    print(f"{process_name}完成!!")
    return jffilteredprobe


def execute_bwa_filtering(args, bwaindex, jffilteredprobe):
    """
    执行BWA过滤
    
    Args:
        args: 命令行参数
        bwaindex: BWA索引路径
        jffilteredprobe: 经过k-mer过滤的探针列表
    
    Returns:
        list: BWA过滤后的探针列表
    """
    # 创建临时探针文件
    tmppbfa = os.path.join(args.saved, f"{os.path.basename(args.input)}_tmp_probe.fa")
    
    with open(tmppbfa, 'w') as tmppbfaio:
        for seqnum, tmppb in enumerate(jffilteredprobe):
            print(f'>seq{seqnum}', file=tmppbfaio)
            print(tmppb, file=tmppbfaio)
    
    # 清理内存
    del jffilteredprobe
    
    print("运行BWA过滤")
    maxxs = int(args.length * args.homology / 100)
    print(f"最大错配数: {maxxs}")
    
    # 执行BWA过滤
    bwafiltedpb = bwa.bwafilter(
        bwabin=args.bwa, 
        reffile=bwaindex, 
        inputfile=tmppbfa, 
        minas=args.length,
        maxxs=maxxs, 
        threadnumber=args.threads
    )
    
    print(f"BWA过滤后探针数量: {len(bwafiltedpb)}")
    print(f"前10个探针: {bwafiltedpb[0:10]}")
    
    return bwafiltedpb


def generate_output_files(args, bwafiltedpb, bwaindex):
    """
    生成输出文件
    
    Args:
        args: 命令行参数
        bwafiltedpb: BWA过滤后的探针列表
        bwaindex: BWA索引路径
    """
    # 准备探针过滤数据
    oligobefortmf = prepare_probe_filtering_data(bwafiltedpb, args)
    
    # 执行探针过滤
    keepedprobe = execute_probe_filtering_final(oligobefortmf)
    
    # 按染色体组织探针
    pbdictbychr = organize_probes_by_chromosome(keepedprobe)
    
    # 生成BED格式输出文件
    generate_bed_files(args, pbdictbychr, bwaindex)


def prepare_probe_filtering_data(bwafiltedpb, args):
    """
    准备探针过滤数据
    
    Args:
        bwafiltedpb: BWA过滤后的探针列表
        args: 命令行参数
    
    Returns:
        list: 准备过滤的探针数据
    """
    oligobefortmf = []
    
    for pbtmp in bwafiltedpb:
        nowpbcounter = {
            'seq': pbtmp,
            'dTm': args.dtm,
            'rprimer': args.primer
        }
        oligobefortmf.append(nowpbcounter)
    
    return oligobefortmf


def execute_probe_filtering_final(oligobefortmf):
    """
    执行最终的探针过滤
    
    Args:
        oligobefortmf: 准备过滤的探针数据
    
    Returns:
        list: 保留的探针列表
    """
    keepedprobe = []
    ctedpb = 0
    oligobefortmflen = len(oligobefortmf)
    
    print(f"探针过滤总数: {oligobefortmflen}")
    
    # 使用多进程进行探针过滤
    pbftpool = Pool()
    
    for (pb, keep) in pbftpool.imap_unordered(probefilter, oligobefortmf):
        if keep:
            keepedprobe.append(pb)
        ctedpb += 1
        
        if ctedpb % 10000 == 0:
            print(f"{ctedpb}/{oligobefortmflen}")
    
    pbftpool.close()
    return keepedprobe


def organize_probes_by_chromosome(keepedprobe):
    """
    按染色体组织探针
    
    Args:
        keepedprobe: 保留的探针列表
    
    Returns:
        dict: 按染色体组织的探针字典
    """
    pbdictbychr = {}
    
    for pb in keepedprobe:
        seq, chro, start = pb.split('\t')
        start = int(start)
        
        if chro not in pbdictbychr:
            pbdictbychr[chro] = {}
        
        pbdictbychr[chro][start] = seq
    
    return pbdictbychr


def generate_bed_files(args, pbdictbychr, bwaindex):
    """
    生成BED格式输出文件
    
    Args:
        args: 命令行参数
        pbdictbychr: 按染色体组织的探针字典
        bwaindex: BWA索引路径
    """
    # 输出文件路径
    tmpbwaftlist = os.path.join(args.saved, f"{os.path.basename(args.input)}.bed")
    alltmpbwaftlist = os.path.join(args.saved, f"{os.path.basename(args.input)}_all.bed")
    seqlenfile = os.path.join(args.saved, f"{os.path.basename(args.input)}.len")
    
    # 生成序列长度文件
    seqlength = bwa.bwareflength(bwabin=args.bwa, reffile=bwaindex)
    with open(seqlenfile, 'w') as seqlenio:
        for seqname, length in seqlength.items():
            print(f"{seqname}\t{length}", file=seqlenio)
    
    # 计算滑动窗口参数
    lenrprimer = len(args.primer) if args.primer else 5
    slidwindow = lenrprimer + args.length
    
    # 生成BED文件
    with open(tmpbwaftlist, 'w') as tmpbwaftlistio, \
         open(alltmpbwaftlist, 'w') as allbwaftlistio:
        
        for chro in pbdictbychr:
            startn = 0
            
            for startnow in sorted(pbdictbychr[chro]):
                endnow = startnow + args.length - 1
                probe_seq = pbdictbychr[chro][startnow]
                
                # 写入所有探针
                print(f"{chro}\t{startnow}\t{endnow}\t{probe_seq}", 
                      file=allbwaftlistio, sep='\t')
                
                # 根据滑动窗口筛选探针
                if startnow > startn + slidwindow:
                    startn = startnow
                    print(f"{chro}\t{startnow}\t{endnow}\t{probe_seq}", 
                          file=tmpbwaftlistio, sep='\t')


def probefilter(nowpbcounter):
    """
    探针过滤器函数
    
    Args:
        nowpbcounter: 包含探针信息的字典
    
    Returns:
        tuple: (探针序列, 是否保留)
    """
    sequence, seqname, start = nowpbcounter['seq'].split('\t')
    dTm = nowpbcounter['dTm']
    rPrimer = nowpbcounter['rprimer']
    
    keep = False
    
    # ATCG含量过滤
    if prefilter.atcg_filter(sequence, 5):
        pass  # 未通过ATCG过滤
    elif dTm == -1:
        # 跳过dTm计算
        keep = True
    else:
        # 执行引物3过滤
        if rPrimer:
            # 使用反向引物进行过滤
            if not primer3_filter.primer3_filter_withRprimer(sequence, rPrimer, dtm=dTm):
                keep = True
        else:
            # 标准引物3过滤
            if not primer3_filter.primer3_filter(sequence, dtm=dTm):
                keep = True
    
    return (nowpbcounter['seq'], keep)


def check_options(parser):
    """
    检查命令行参数的有效性
    
    Args:
        parser: 参数解析器
    
    Returns:
        Namespace: 验证后的参数
    """
    args = parser.parse_args()
    
    # 检查BWA
    args.bwa = validate_bwa_path(args.bwa)
    
    # 检查Jellyfish
    args.jellyfish = validate_jellyfish_path(args.jellyfish)
    
    # 检查输入文件
    validate_input_files(args, parser)
    
    # 检查参数范围
    validate_parameter_ranges(args, parser)
    
    # 检查输出目录
    validate_output_directory(args)
    
    # 打印验证信息
    print_validation_info(args)
    
    return args


def validate_bwa_path(bwa_path):
    """
    验证BWA路径
    
    Args:
        bwa_path: BWA路径
    
    Returns:
        str: 有效的BWA路径
    """
    if bwa_path:
        if not os.path.exists(bwa_path):
            raise ValueError("无法定位BWA，请输入BWA的完整路径")
        
        bwaversion = bwa.bwaversion(bwa_path)
        if bwaversion == 'None':
            raise ValueError("无法定位BWA，请输入BWA的完整路径")
    else:
        # 尝试在PATH中查找BWA
        bwa_paths = which('bwa')
        if bwa_paths:
            bwaversion = bwa.bwaversion(bwa_paths[0])
            if bwaversion != 'None':
                bwa_path = bwa_paths[0]
            else:
                raise ValueError("无法定位BWA，请输入BWA的完整路径")
        else:
            raise ValueError("无法定位BWA，请输入BWA的完整路径")
    
    return bwa_path


def validate_jellyfish_path(jellyfish_path):
    """
    验证Jellyfish路径
    
    Args:
        jellyfish_path: Jellyfish路径
    
    Returns:
        str: 有效的Jellyfish路径
    """
    if jellyfish_path:
        if not os.path.exists(jellyfish_path):
            raise ValueError("无法定位Jellyfish，请输入Jellyfish的完整路径")
        
        jellyfishversion = jellyfish.jfversion(jellyfish_path)
        if jellyfishversion == 'None':
            raise ValueError("无法定位Jellyfish，请输入Jellyfish的完整路径")
    else:
        # 尝试在PATH中查找Jellyfish
        jellyfish_paths = which('jellyfish')
        if jellyfish_paths:
            jellyfishversion = jellyfish.jfversion(jellyfish_paths[0])
            if jellyfishversion != 'None':
                jellyfish_path = jellyfish_paths[0]
            else:
                raise ValueError("无法定位Jellyfish，请输入Jellyfish的完整路径")
        else:
            raise ValueError("无法定位Jellyfish，请输入Jellyfish的完整路径")
    
    return jellyfish_path


def validate_input_files(args, parser):
    """
    验证输入文件
    
    Args:
        args: 命令行参数
        parser: 参数解析器
    """
    if not os.path.exists(args.genome):
        print("无法定位基因组文件，请输入基因组文件。")
        parser.print_help()
        sys.exit(1)
    
    if not os.path.exists(args.input):
        print("无法定位输入文件，请输入输入文件。")
        parser.print_help()
        sys.exit(1)


def validate_parameter_ranges(args, parser):
    """
    验证参数范围
    
    Args:
        args: 命令行参数
        parser: 参数解析器
    """
    if args.homology > 100 or args.homology < 50:
        print("请将同源性设置在50到100之间。")
        parser.print_help()
        sys.exit(1)
    
    if args.dtm > 37 or args.dtm < -1:
        print("请将dTm设置在0到37之间。")
        parser.print_help()
        sys.exit(1)
    
    if args.skipdtm:
        args.dtm = -1
        print("跳过dTm计算")
    
    if args.length > 50:
        args.dtm = -1
        print("寡核苷酸长度超过50nt，跳过dTm计算")
    
    if args.step < 1:
        print("步长必须大于1。")
        parser.print_help()
        sys.exit(1)


def validate_output_directory(args):
    """
    验证输出目录
    
    Args:
        args: 命令行参数
    """
    if os.path.exists(args.saved):
        if not args.docker:
            print(f"输出文件夹 {args.saved} 已存在。")
            print('''
                按Y使用此文件夹输出文件。此文件夹中的所有内容将被删除。
                按N并设置-s/--saved为不同的文件夹:
                ''')
            
            while True:
                char = getch()
                if char.lower() in ("y", "n"):
                    print(char)
                    if char.lower() == 'n':
                        sys.exit(1)
                    break
    else:
        os.mkdir(args.saved)


def print_validation_info(args):
    """
    打印验证信息
    
    Args:
        args: 命令行参数
    """
    print("#" * 40)
    print(f"BWA版本: {args.bwa}")
    print(f"Jellyfish版本: {args.jellyfish}")
    print(f"基因组文件: {args.genome}")
    print(f"输入文件: {args.input}")
    print(f"5'标记反向引物: {args.primer}")
    print(f"结果输出文件夹: {os.path.realpath(args.saved)}")
    print(f"线程数: {args.threads}")
    print(f"同源性: {args.homology}")
    print(f"dTm: {args.dtm}")
    
    if args.dtm == -1:
        print("警告: 跳过dTm计算。")
    
    print("#" * 40)


def getch():
    """
    获取用户输入（用于是/否选择）
    
    Returns:
        str: 用户输入的字符
    """
    import sys, tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch


def which(filename):
    """
    在PATH中查找可执行文件
    
    Args:
        filename: 文件名
    
    Returns:
        list: 找到的文件路径列表
    """
    locations = os.environ.get("PATH").split(os.pathsep)
    candidates = []
    
    for location in locations:
        candidate = os.path.join(location, filename)
        if os.path.isfile(candidate):
            candidates.append(candidate)
    
    return candidates


def get_options():
    """
    创建命令行参数解析器
    
    Returns:
        ArgumentParser: 参数解析器
    """
    parser = argparse.ArgumentParser(
        description="Chorus2 寡核苷酸FISH探针设计软件", 
        prog='Chorus2',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例:\n"
               "  Chorus2 -i TAIR10_chr_all.fas -g TAIR10_chr_all.fas -t 4 \\ \n"
               "          -j /opt/software/jellyfish/bin/jellyfish -b /opt/software/bwa/bwa -s sample"
    )
    
    # 版本信息
    parser.add_argument('--version', action='version', version='%(prog)s 2.1')
    
    # 软件路径参数
    parser.add_argument('-j', '--jellyfish', dest='jellyfish', 
                       help='Jellyfish软件安装路径')
    parser.add_argument('-b', '--bwa', dest='bwa', 
                       help='BWA软件安装路径')
    
    # 必需参数
    parser.add_argument('-g', '--genome', dest='genome', required=True,
                       help='FASTA格式基因组文件，应包含基因组的所有序列')
    parser.add_argument('-i', '--input', dest='input', required=True,
                       help='FASTA格式输入文件，可以是全基因组、染色体或基因组的一个区域')
    
    # 可选参数
    parser.add_argument('-s', '--save', dest='saved', default='probes',
                       help='保存结果的输出文件夹')
    parser.add_argument('-p', '--primer', dest='primer', default='',
                       help='PCR反应的特异性5\'标记反向引物。例如: CGTGGTCGCGTCTCA (默认无)')
    parser.add_argument('-t', '--threads', dest='threads', default=1, type=int,
                       help='使用的线程或CPU数量 (默认: 1)')
    parser.add_argument('-l', '--length', dest='length', default=45, type=int,
                       help='探针长度 (默认: 45)')
    parser.add_argument('--homology', dest='homology', default=75, type=float,
                       help='目标序列与探针之间的最大同源性(%%)，范围50-100 (默认: 75)')
    parser.add_argument('-d', '--dtm', dest='dtm', default=10, type=float,
                       help='dTm的最小值（杂交Tm - 发夹Tm），范围0-37 (默认: 10)')
    parser.add_argument('--skipdtm', dest='skipdtm', default=False, type=bool,
                       help='跳过dTm计算，适用于长度超过50的寡核苷酸')
    parser.add_argument('--step', dest='step', default=5, type=int,
                       help='滑动窗口中k-mer搜索的步长，步长>=1 (默认: 5)')
    parser.add_argument('--docker', help='仅用于Chorus的Docker版本', default=False)
    parser.add_argument('--ploidy', dest='ploidy', default=2, type=int,
                       help='给定基因组的倍性（测试版本）(默认: 2)')
    
    return parser


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("用户中断\n")
        sys.exit(0)
    except Exception as e:
        print(f"程序执行错误: {e}")
        sys.exit(1)
