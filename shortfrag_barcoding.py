#!/usr/bin/env python
import time
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import subprocess
import os
from shutil import copyfile
from shutil import move
from glob import glob
from multiprocessing import Pool
import multiprocessing
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


parser = OptionParser()

parser.add_option("-5", "--i5_map", dest = "i5_map", type = str, help = "i5 index map in tab-delimited format: <plate_name><\\t><i5_index><\\t><locus_name>")
parser.add_option("-7", "--i7_map", dest = "i7_map", type = str, help = "i7 index map in tab-delimited format: <well_name><\\t><i7_index>")
parser.add_option("-x", "--i5_index_file", dest = "i5_index_file", type = str, help = "gzipped i5 index read file (read 3)")
parser.add_option("-s", "--i7_index_file", dest = "i7_index_file", type = str, help = "gzipped i7 index read file (read 2)")
parser.add_option("-1", "--R1_file", dest = "r1_file", type = str, help = "gzipped read 1 file")
parser.add_option("-2", "--R2_file", dest = "r2_file", type = str, help = "gzipped read 2 file", default = "None")
parser.add_option("-d", "--adapter_file", dest = "adapter_file", type = str, help = "fasta file of adapter sequences")
parser.add_option("-p", "--num_threads", dest = "num_threads", type = int, default = 1, help = "Number of threads to use for sequencing trimming and assembly [default = 1]")
parser.add_option("-c", "--config_file", dest = "config_file", default = "config.txt", type = str, help = "Config file with program paths [default = config.txt]")
parser.add_option("-o", "--output_db", dest = "output_db", type = str, default = "output.txt", help = "Output database file name [default = output.txt]")
parser.add_option("-e", "--email", dest = "email", type = str, help = "Email address. Required if you want to submit sequences to NCBI BLAST. If you don't set this then BLAST will not be run. If you do set it, then it will.", default = None)
parser.add_option("-g", "--collections_file", dest = "collections_file", type = str, help = "Three column file with <site><\\t><plate><\\t><well>", default = None)

(options, args) = parser.parse_args()

USE_NCBI_BLAST_SERVER = False
ADD_SITE_CODE = False

if options.email != None:
    Entrez.email = options.email
    USE_NCBI_BLAST_SERVER = True
if options.collections_file != None:
    ADD_SITE_CODE = True

def read_config():
    reader = open(options.config_file, 'rU')
    path_dic = {}
    for line in reader:
        cur_line = line.split()
        path_dic[cur_line[0]] = cur_line[1]
    return path_dic

def reformat_i5_map(i5_map):
    reader = open(i5_map, 'rU')
    i5_map_name = "temp_i5.txt"
    outfile = open(i5_map_name, 'w')
    for line in reader:
        cur_line = line.split()
        outfile.write("%s_%s\t%s\n" % (cur_line[0], cur_line[2], cur_line[1]))
    outfile.close()
    return i5_map_name

def demult_i5(i5_index_file, i7_index_file, r1_file, r2_file, i5_map, path_dic):
    new_i5_map = reformat_i5_map(i5_map)
    cmd = [path_dic["multx_path"], "-B", new_i5_map, "%s" % i5_index_file, "%s" % i7_index_file, "%s" % r1_file, "%s" % r2_file, "-o", "demult_i5/%_i5.fastq.gz", "-o", "demult_i5/%_i7.fastq.gz", "-o", "demult_i5/%_R1.fastq.gz", "-o", "demult_i5/%_R2.fastq.gz"]
    if r2_file == "None":
        cmd = [path_dic["multx_path"], "-B", new_i5_map, "%s" % i5_index_file, "%s" % i7_index_file, "%s" % r1_file, "-o", "demult_i5/%_i5.fastq.gz", "-o", "demult_i5/%_i7.fastq.gz", "-o", "demult_i5/%_R1.fastq.gz"]
    if not os.path.exists("demult_i5"):
        os.mkdir("demult_i5")
    FNULL = open(os.devnull, 'w')
    with open("demult_i5/fastq_multx.log", 'w') as outfile:
        subprocess.call(cmd, stderr = FNULL, stdout = outfile)
    outfile.close()
    reader = open(i5_map, 'rU')
    for line in reader:
        cur_line = line.split()
        if not os.path.exists("./plates"):
            os.mkdir("plates")
        if not os.path.exists("./plates/%s" % cur_line[0]):
            os.mkdir("./plates/%s" % cur_line[0])
        if not os.path.exists("./plates/%s/%s" % (cur_line[0], cur_line[2])):
            os.mkdir("./plates/%s/%s" % (cur_line[0], cur_line[2]))
        if r2_file == "None":
            for suf in ["i7", "R1"]:
                move("demult_i5/%s_%s_%s.fastq.gz" % (cur_line[0], cur_line[2], suf), "./plates/%s/%s" % (cur_line[0], cur_line[2]))
        else:
            for suf in ["i7", "R1", "R2"]:
                move("demult_i5/%s_%s_%s.fastq.gz" % (cur_line[0], cur_line[2], suf), "./plates/%s/%s" % (cur_line[0], cur_line[2]))
            
def demult_i7(i5_file, well_map, path_dic):
    reader = open(i5_file, 'rU')
    for line in reader:
        cur_line = line.split()
        if not os.path.exists("./plates/%s/%s/demult" % (cur_line[0], cur_line[2])):
            os.mkdir("./plates/%s/%s/demult" % (cur_line[0], cur_line[2]))
        if options.r2_file == "None":
            cmd = [path_dic["multx_path"], "-B", well_map, "./plates/%s/%s/%s_%s_i7.fastq.gz" % (cur_line[0], cur_line[2], cur_line[0], cur_line[2]), "./plates/%s/%s/%s_%s_R1.fastq.gz" % (cur_line[0], cur_line[2], cur_line[0], cur_line[2]), "-o", "./plates/%s/%s/demult/%s_i7.fastq" % (cur_line[0], cur_line[2], "%"), "-o", "./plates/%s/%s/demult/%s_R1.fastq" % (cur_line[0], cur_line[2], "%")]
        else:
            cmd = [path_dic["multx_path"], "-B", well_map, "./plates/%s/%s/%s_%s_i7.fastq.gz" % (cur_line[0], cur_line[2], cur_line[0], cur_line[2]), "./plates/%s/%s/%s_%s_R1.fastq.gz" % (cur_line[0], cur_line[2], cur_line[0], cur_line[2]), "./plates/%s/%s/%s_%s_R2.fastq.gz" % (cur_line[0], cur_line[2], cur_line[0], cur_line[2]), "-o", "./plates/%s/%s/demult/%s_i7.fastq" % (cur_line[0], cur_line[2], "%"), "-o", "./plates/%s/%s/demult/%s_R1.fastq" % (cur_line[0], cur_line[2], "%"), "-o", "./plates/%s/%s/demult/%s_R2.fastq" % (cur_line[0], cur_line[2], "%")]
        FNULL = open(os.devnull, 'w')        
        with open("./plates/%s/%s/demult/fastq_multx.log" % (cur_line[0], cur_line[2]), 'w') as outfile:
            subprocess.call(cmd, stdout = outfile, stderr = FNULL)
        outfile.close()

def merge_pairs_worker(param_list):
    plate = param_list[0]
    locus = param_list[1]
    sample = param_list[2]
    path_dic = param_list[3]
    cmd = [path_dic["usearch_path"], "-fastq_mergepairs", "./plates/%s/%s/demult/%s_R1.fastq" % (plate, locus, sample), "-reverse", "./plates/%s/%s/demult/%s_R2.fastq" % (plate, locus, sample), "-fastqout", "./plates/%s/%s/merged/%s_merged.fq" % (plate, locus, sample), "-relabel", "%s_" % sample, "-fastqout_notmerged_fwd", "./plates/%s/%s/merged/%s_notmerged_R1.fq" % (plate, locus, sample),"-fastqout_notmerged_rev", "./plates/%s/%s/merged/%s_notmerged_R2.fq" % (plate, locus, sample), "-report", "./plates/%s/%s/merged/%s_report.txt" % (plate, locus, sample), "-fastq_minmergelen", "360", "-fastq_maxdiffs", "20"]
    subprocess.call(cmd)
    seq_count = 0
    reader = SeqIO.parse("./plates/%s/%s/merged/%s_merged.fq" % (plate, locus, sample), format = 'fastq')
    for rec in reader:
        seq_count += 1
    if seq_count == 0:
        return
    cmd = [path_dic["usearch_path"], "-fastx_truncate", "./plates/%s/%s/merged/%s_merged.fq" % (plate, locus, sample), "-stripleft", "26", "-stripright", "26", "-fastqout", "./plates/%s/%s/merged/%s_stripped.fq" % (plate, locus, sample)]
    subprocess.call(cmd)
    cmd = [path_dic["usearch_path"], "-fastq_filter", "./plates/%s/%s/merged/%s_stripped.fq" % (plate, locus, sample), "-fastq_maxee", "0.5", "-fastaout", "./plates/%s/%s/merged/%s_filtered.fa" % (plate, locus, sample)]
    subprocess.call(cmd)
    FNULL = open(os.devnull, 'w')
    cmd = [path_dic["usearch_path"], "-fastx_uniques", "./plates/%s/%s/merged/%s_filtered.fa" % (plate, locus, sample), "-fastaout", "./plates/%s/%s/derep/%s_derep.fa" % (plate, locus, sample), "-minuniquesize", "2", "-sizeout", "-relabel", "%s-%s_" % (plate, sample)]
    subprocess.call(cmd, stderr = FNULL)
    if not os.path.exists("./plates/%s/%s/derep/%s_derep.fa" % (plate, locus, sample)):
        return
    seq_count = 0
    reader = SeqIO.parse("./plates/%s/%s/derep/%s_derep.fa" % (plate, locus, sample), format = 'fasta')
    for rec in reader:
        seq_count += 1
    if seq_count == 0:
        return
    cmd = [path_dic["usearch_path"], "-cluster_smallmem", "./plates/%s/%s/derep/%s_derep.fa" % (plate, locus, sample), "-id", "1.0", "-sizein", "-sizeout", "-sortedby", "other", "-relabel", "%s-%s_" % (plate, sample), "-centroids", "./plates/%s/%s/derep/%s_cluster.fa" % (plate, locus, sample)]
    subprocess.call(cmd, stderr = FNULL)
    most_abundant_read("./plates/%s/%s/derep/%s_cluster.fa" % (plate, locus, sample), "./plates/%s/%s/derep/%s_common.fa" % (plate, locus, sample), path_dic)

def most_abundant_read(infile, outfile, path_dic):
    reader = SeqIO.parse(infile, format = 'fasta')
    seq_dic = {}
    size_dic = {}
    for rec in reader:
        seq_dic[rec.id] = str(rec.seq)
        cur_size = int(rec.id.split("size=")[1][0:-1])
        size_dic[rec.id] = cur_size
    biggest = ""
    commonest = 0
    for k, v in size_dic.items():
        if v > commonest:
            commonest = v
            biggest = k
    common_seqs = {}
    for k, v in seq_dic.items():
        if size_dic[k] >= 0.5 * commonest:
            common_seqs[k] = v
    output = open(outfile, 'w')
    for k, v in common_seqs.items():
        output.write(">%s\n%s\n" % (k, v))
    output.close()

def merge_pairs(i5_file, i7_file, path_dic):
    reader = open(i5_file, 'rU')
    plate_list = []
    locus_list = []
    id_list = []
    work_list = []
    pool = multiprocessing.Pool(processes = options.num_threads)

    for line in reader:
        cur_line = line.split()
        reader = open(i7_file, 'rU')
        for i7_line in reader:
            cur_i7_line = i7_line.split()
            plate_list.append(cur_line[0])
            locus_list.append(cur_line[2])
            id_list.append(cur_i7_line[0])
            if not os.path.exists("./plates/%s/%s/merged" % (cur_line[0], cur_line[2])):
                os.mkdir("./plates/%s/%s/merged" % (cur_line[0], cur_line[2]))
            if not os.path.exists("./plates/%s/%s/derep" % (cur_line[0], cur_line[2])):
                os.mkdir("./plates/%s/%s/derep" % (cur_line[0], cur_line[2]))
            if not os.path.exists("./plates/temp"):
                os.mkdir("./plates/temp")
            work_list.append([cur_line[0], cur_line[2], cur_i7_line[0], path_dic])
    pool.map_async(merge_pairs_worker, work_list).get(9999999)

def compile_seqs(i5_file, i7_file, output_db):
    outfile = open(output_db, 'w')
    reader = open(i5_file, 'rU')
    plate_list = []
    locus_list = []
    id_list = []
    for line in reader:
        cur_line = line.split()
        reader = open(i7_file, 'rU')
        for i7_line in reader:
            cur_i7_line = i7_line.split()
            plate = cur_line[0]
            locus = cur_line[2]
            sample = cur_i7_line[0]
            if not os.path.exists("./plates/%s/%s/derep/%s_common.fa" % (plate, locus, sample)):
                continue
            reader = SeqIO.parse("./plates/%s/%s/derep/%s_common.fa" % (plate, locus, sample), format = 'fasta')
            for rec in reader:
                outfile.write(">%s\n%s\n" % (rec.id, str(rec.seq)))
    outfile.close()

def get_best_blast_hit(infile, path_dic):
    reader = SeqIO.parse(infile, format = 'fasta')
    seq_dic = {}
    for rec in reader:
        seq_dic[rec.id[0:-1]] = str(rec.seq)
    inbase = infile.rsplit(".", 1)[0]
    blastn_cline = NcbiblastnCommandline(query=infile, db = path_dic["ncbi_nt_path"], outfmt = 5, out = "%s_blastn.xml" % inbase, num_threads = options.num_threads, max_target_seqs = 1)
    stdout, stderr = blastn_cline()
    result_handle = open("%s_blastn.xml" % inbase)
    outfile = open("%s_bestblast.fa" % inbase, 'w')
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        cur_query = blast_record.query
        if len(blast_record.alignments) < 1:
            outfile.write(">%s;no_blast_hit\n%s\n" % (cur_query, seq_dic[cur_query]))
            continue
        besthit = blast_record.alignments[0]
        blast_seq = besthit.hsps[0].sbjct
        hitlen = besthit.hsps[0].align_length
        idents = float(besthit.hsps[0].identities)
        percid = round(idents / hitlen, 4)
        percid = 100*percid
        gi = besthit.title.split("|")[1]
        ncbiResult = Entrez.efetch(db="nucleotide", id=gi, rettype="gb", retmode="text")
        time.sleep(5) #don't overwhelm NCBI
        ncbiResultRec = SeqIO.read(ncbiResult, "gb")
        cur_species = ncbiResultRec.annotations["organism"].replace(" ", "_")
        
        if "Hymenoptera" in ncbiResultRec.annotations["taxonomy"]:
            outfile.write(">%s;%s;%s;%s;%s;%s\n%s\n" % (cur_query, gi, cur_species, percid, hitlen, len(seq_dic[cur_query]), seq_dic[cur_query]))
        else:
            outfile.write(">%s;%s;%s;%s;%s;%s;not_hymenoptera\n%s\n" % (cur_query, gi, cur_species, percid, hitlen, len(seq_dic[cur_query]), seq_dic[cur_query]))
        outfile.flush()
    outfile.close()

def add_site_code(infile, collections_file):
    reader = open(collections_file, 'rU')
    site_dic = {}
    for line in reader:
        if line.startswith("#"):
            continue
        cur_line = line.strip().split("\t")
        cur_id = cur_line[1] + "-" + cur_line[2]
        cur_site = cur_line[0]
        site_dic[cur_id] = cur_site
    outfile = open(infile.rsplit(".", 1)[0] + "_codes.fa", 'w')
    reader = SeqIO.parse(infile, format = 'fasta')
    for rec in reader:
        cur_info = rec.id.split(";")
        cur_id = cur_info[0].split("_")[0]

        outfile.write(">%s;%s;%s\n%s\n" % (cur_info[0], site_dic[cur_id], ";".join(cur_info[1:]), str(rec.seq)))
    outfile.close()
    
def main():
    path_dic = read_config()
    demult_i5(options.i5_index_file, options.i7_index_file, options.r1_file, options.r2_file, options.i5_map, path_dic)
    demult_i7(options.i5_map, options.i7_map, path_dic)
    merge_pairs(options.i5_map, options.i7_map, path_dic)
    compile_seqs(options.i5_map, options.i7_map, options.output_db)
    if USE_NCBI_BLAST_SERVER:
        get_best_blast_hit(options.output_db, path_dic)
        if ADD_SITE_CODE:
            inbase = options.output_db.rsplit(".", 1)[0]
            add_site_code(inbase + "_bestblast.fa", options.collections_file)
    else:
        if ADD_SITE_CODE:
              add_site_code(option.output_db, options.collections_file)
        
if __name__ == '__main__':
    main()
