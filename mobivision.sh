#!/bin/bash  

#SBATCH --job-name=mobivision.fmt_gtf.mkindex.quantify   
#SBATCH --partition=debug  
#SBATCH --nodes=1  
#SBATCH -n 4  
#SBATCH --cpus-per-task=3  
#SBATCH --time=7-00:00:00  
#SBATCH --mem-per-cpu=2G  
#SBATCH --output=/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/mobivision.fmt_gtf.quantify.out  
#SBATCH --error=/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/mobivision.fmt_gtf.quantify.err  

# 定义输入和输出的基础路径  
input_base="/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/seqRawdata"  
output_base="/mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/QuantifyRawdata"  

# 定义要遍历的文件夹  
folders=("PJ2410100101-6例-150G")  

# 循环遍历每个文件夹  
for folder in "${folders[@]}"; do  
    # 获取当前文件夹的完整路径  
    rawdata_folder="${input_base}/${folder}/Rawdata"  

    # 确保 Rawdata 文件夹存在  
    if [ -d "$rawdata_folder" ]; then  
        # 循环遍历 Rawdata 文件夹中的样本文件夹  
        for sample in "$rawdata_folder"/*; do  
            if [ -d "$sample" ]; then  
                sample_name=$(basename "$sample")  
                output_file="${output_base}/${sample_name}"  

                # 运行 mobivision quantify 命令  
                mobivision quantify -i /mnt2/wanggd_group/zjj/BGCscRNA/mobivision/RefData/CanFam3.1 \
                -t 12 \
                -f "$sample" \
                -o "$output_file"  

                # 检查 mobivision 的执行状态  
                if [ $? -ne 0 ]; then  
                    echo "mobivision quantify 失败，样本: $sample" >> /mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/mobivision.fmt_gtf.quantify.err  
                fi  
            fi  
        done  
    else  
        echo "目录不存在: $rawdata_folder" >> /mnt2/wanggd_group/zjj/BGCscRNA/ZengMin/mobivision.fmt_gtf.quantify.err  
    fi  
done



#mobivision fmt_gtf -i Canis_lupus_familiaris.CanFam3.1.104.gtf.gz \
#  -o Canis_lupus_familiaris.CanFam3.1.104.filtered.gtf \
#  -t protein_coding \
# -t lncRNA \
# -t antisense \
# -t IG_LV_gene \
# -t IG_V_gene \
# -t IG_V_pseudogene \
# -t IG_D_gene \
# -t IG_J_gene \
# -t IG_J_pseudogene \
# -t IG_C_gene \
# -t IG_C_pseudogene \
# -t TR_V_gene \
# -t TR_V_pseudogene \
# -t TR_D_gene \
# -t TR_J_gene \
# -t TR_J_pseudogene \
# -t TR_C_gene

#mobivision mkindex -n CanFam3.1 \
#-f Canis_lupus_familiaris.CanFam3.1.dna_sm.toplevel.fa.gz \
#-g Canis_lupus_familiaris.CanFam3.1.104.filtered.gtf


