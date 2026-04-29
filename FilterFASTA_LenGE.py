#!/usr/bin/env python3
import sys

def filter_fasta(input_file, output_file, min_length):
    count_total = 0
    count_passed = 0
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        current_id = None
        current_seq = []
        
        for line in fin:
            line = line.rstrip('\n')
            
            if line.startswith('>'):
                # 遇到新的序列头，先处理上一条
                if current_id is not None:
                    seq_len = len(''.join(current_seq))
                    count_total += 1
                    if seq_len >= min_length:
                        count_passed += 1
                        fout.write(f">{current_id}\n")
                        fout.write(''.join(current_seq) + '\n')
                
                # 开始记录新序列
                # 去掉 > 符号，保留完整描述行
                current_id = line[1:]
                current_seq = []
            else:
                # 累积序列行（忽略空行）
                if line.strip():
                    current_seq.append(line.strip())
        
        # 处理最后一条序列（文件末尾没有 > 触发）
        if current_id is not None:
            seq_len = len(''.join(current_seq))
            count_total += 1
            if seq_len >= min_length:
                count_passed += 1
                fout.write(f">{current_id}\n")
                fout.write(''.join(current_seq) + '\n')
    
    print(f"✅ Done!")
    print(f"   Total count: {count_total}")
    print(f"   Min length:     ≥ {min_length} aa")
    print(f"   Passed count:   {count_passed}")
    print(f"   Filtered out:   {count_total - count_passed}")
    print(f"   Output file: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_fasta_by_length.py <input.fasta> <output.fasta> <min_length>")
        print("Example: python filter_fasta_by_length.py input.fasta output.fasta 9")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    min_len = int(sys.argv[3])
    
    filter_fasta(input_fasta, output_fasta, min_len)
