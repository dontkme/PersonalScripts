#!/usr/bin/env python3
import sys

def rename_fasta(input_file, output_file, prefix="PEP"):
    count = 0
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            # 去除行尾换行符
            line = line.rstrip('\n')
            
            if line.startswith('>'):
                count += 1
                # 生成7位数字编号，例如 0000001
                new_id = f"{prefix}{count:07d}"
                
                # 提取原名称（去掉 > 号，保留整行作为描述）
                # 原名称可能本身包含空格和描述，我们整体保留
                original_name = line[1:] 
                
                # 组合新行：>PEP0000001 原名称
                fout.write(f">{new_id} {original_name}\n")
            else:
                # 序列行原样输出
                fout.write(f"{line}\n")
                
    print(f"✅ 处理完成！共重命名 {count} 条序列。")
    print(f"输出文件: {output_file}")

if __name__ == "__main__":
    # 使用方法：python rename_fasta.py 输入.fasta 输出.fasta
    if len(sys.argv) != 3:
        print("Usage: python rename_fasta.py <input.fasta> <output.fasta>")
        sys.exit(1)
        
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    rename_fasta(input_fasta, output_fasta)
