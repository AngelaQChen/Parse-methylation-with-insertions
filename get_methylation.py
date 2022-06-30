import pysam
import argparse

parser = argparse.ArgumentParser(description = 
                                 "Parse sorted and indexed BAM file for methylation modification in specified region.")
parser.add_argument("-b", "--bam", metavar='BAM', help = 
                    "indexed and sorted BAM file with methylation data from primrose")
parser.add_argument("-r", "--region", metavar='STR', help = 
                    "CHROM:START-END format for region")
parser.add_argument("-o", "--out", metavar='STR', default="output", help = 
                    "prefix of the output files [output]")
args = parser.parse_args()

'''
slc9a3 VNTR
Here is the region: chr5:393,462-696,129
and here are the two VNTRs we are interested:
chr5:474,089-474,863
chr5:662,841-663,226
'''

chrom = args.region[:args.region.index(":")]
start = int(args.region[args.region.index(":") +1 : args.region.index("-")])
end  = int(args.region[args.region.index("-") +1 :])

out_prefix = args.out

#dictionary of each position and associated methylation quality scores
hap1_modpos = {} #{ref_pos : [scores] }
hap2_modpos = {}
combined_modpos = {}

#read spans, for calculating coverage 
#list of tuples
hap1_span = [] #[(start, end)...]
hap2_span = []
combined_span = []

hap1_ins_span = []
hap2_ins_span = []
combined_ins_span = []

#ins start sight based on original sequence (not strand flipped)
#aligned_pairs is a list of tuples [(qpos, rpos)...]
def ins_start_site(aligned_pairs, ref_start, hap):
  ins_start_site = [[ref_start, 0, 0]] #[(rpos, qpos)]
  prev_ref_pos = ref_start
  prev_query_pos = 0
  ins_span_start = None
  inside_ins = False
  ins_length = 0
  for p in aligned_pairs:
    #fresh insertion
    #p = (qpos, rpos)
    if p[1] == None and prev_ref_pos != None:
      ins_start_site.append([prev_ref_pos,prev_query_pos])
      ins_span_start = prev_ref_pos
      inside_ins = True
      ins_length = 0
    #end of ins
    elif inside_ins and p[1] != None :
      ins_start_site[len(ins_start_site)-1].append(ins_length + 1) #append the length of the insertion
      if hap == 1 :
        hap1_ins_span.append((ins_span_start, p[1]))
      elif hap == 2:
        hap2_ins_span.append((ins_span_start, p[1]))
      combined_ins_span.append((ins_span_start, p[1]))
      inside_ins = False
    #continuing inside insertion
    elif inside_ins and p[1] == None:
      ins_length += 1
    #continuing within aligned regions
    elif inside_ins == False and p[1] != None:
      pass
    

    prev_ref_pos = p[1]
    if p[0] != None:
      prev_query_pos = p[0]
     
  ins_start_site[len(ins_start_site)-1].append(ins_length + 1) #sometimes the last one is soft clipped
  return(ins_start_site) #[[rpos, qpos, ins_length]...]

#for forward_strand
#modified_bases_forward is [(qpos, qual)...]
def count_cpg_forward(modified_bases_forward, aligned_pairs, hap, ins_start, ref_start):
  ins_start = sorted(ins_start)
  prev_ref_mod_site = 0
  pairs = dict(aligned_pairs) #{qpos: rpos...}
  for base in modified_bases_forward: #base is (qpos, qual)

    if hap == 1:
      if pairs[base[0]] != None:
        try: hap1_modpos[pairs[base[0]]].append(base[1])
        except: hap1_modpos[pairs[base[0]]] = [base[1]]
      else:
        i = 0
        prev_ref_base = None
        while prev_ref_base ==  None:
          try: prev_ref_base = pairs[base[0]-i] #start at pairs[base[0]] then move back 
          except: prev_ref_base = ref_start #means the it was at the start of the read
          i += 1
        for i in ins_start: #ins_start is [rpos, qpos, length]
          #when the start is bigger stop
          if i[0] >= prev_ref_base:
            #start pos of the ins: ins_start+extra_bp
            try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1])
            #error most likely i[1] = None
            except: 
              print(i, base)
              #ins_pos_relative = str(i[0]) + "+"
            if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
            try: hap1_modpos[ins_pos_relative].append(base[1])
            except: hap1_modpos[ins_pos_relative] = [base[1]]
            break

    elif hap == 2:
      if pairs[base[0]] != None:
        try: hap2_modpos[pairs[base[0]]].append(base[1])
        except: hap2_modpos[pairs[base[0]]] = [base[1]]
      else:
        i = 0
        prev_ref_base = None
        while prev_ref_base ==  None:
          try: prev_ref_base = pairs[base[0]-i]
          except: prev_ref_base = ref_start
          i += 1
        for i in ins_start:
            #when the start is bigger stop
            if i[0] >= prev_ref_base: 
              #start pos of the ins: ins_start+extra_bp
              try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1])
              except: ins_pos_relative = str(i[0]) + "+" 
              if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
              try: hap2_modpos[ins_pos_relative].append(base[1])
              except: hap2_modpos[ins_pos_relative] = [base[1]]
              break

    #do this for all reads phased or unphased
    if pairs[base[0]] != None:
      try: combined_modpos[pairs[base[0]]].append(base[1])
      except: combined_modpos[pairs[base[0]]] = [base[1]]
      #remember prev modification site in case next is within an insertion
      #assign new prev_ref_mod_site after done with all the insertion locations
      prev_ref_mod_site = pairs[base[0]] 
    else:
      i = 0
      prev_ref_base = None
      while prev_ref_base ==  None:
        try: prev_ref_base = pairs[base[0]-i]
        except: prev_ref_base = ref_start
        i += 1
      for i in ins_start:
            #when the start is bigger stop
            if i[0] >= prev_ref_base: 
              #start pos of the ins: ins_start+extra_bp
              try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1])
              except: ins_pos_relative = str(i[0]) + "+"
              if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
              try: combined_modpos[ins_pos_relative].append(base[1])
              except: combined_modpos[ins_pos_relative] = [base[1]]
              break
              
#for reverse_strand, minus 1bp to the postion for strand flip
def count_cpg_reverse(modified_bases_reverse, aligned_pairs, hap, ins_start, ref_start):
  ins_start = sorted(ins_start)
  prev_ref_mod_site = 0
  pairs = dict(aligned_pairs)
  for base in modified_bases_reverse: #base is (qpos, qual)

    if hap == 1:
      if pairs[base[0]] != None:
        ins_start_rpos = [x[0] for x in ins_start]
        if pairs[base[0]] == ref_start or pairs[base[0]] - 1 not in ins_start_rpos:
          try: hap1_modpos[pairs[base[0]] -1].append(base[1])
          except: hap1_modpos[pairs[base[0]] -1] = [base[1]] 
        #when the CpG overlaps the last pos of ins and next pos with ref base
        elif pairs[base[0]] -1 in ins_start_rpos:
          ins_start_lengths = [x[2] for x in ins_start]
          ins_pos_relative = str(pairs[base[0]] - 1) + "+" + str(ins_start_lengths[ins_start_rpos.index(pairs[base[0]] - 1)]) # (base-1) + (len(ins))
          try: hap1_modpos[ins_pos_relative].append(base[1])
          except: hap1_modpos[ins_pos_relative] = [base[1]]
      else:
        i = 0
        prev_ref_base = None
        while prev_ref_base ==  None:
          try: prev_ref_base = pairs[base[0]-i]
          except: prev_ref_base = ref_start #when you can't go back in position
          i += 1
        for i in ins_start:
          #when the start is bigger stop
          if i[0] >= prev_ref_base: # +1 bc we -1 from prev_red_mod_site
            #start pos of the ins: ins_start+extra_bp
            
            try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1] -1) 
            except: ins_pos_relative = str(i[0]) + "+"
            if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
            try: hap1_modpos[ins_pos_relative].append(base[1])
            except: hap1_modpos[ins_pos_relative] = [base[1]]
            break

    elif hap == 2:
      if pairs[base[0]] != None:
        ins_start_rpos = [x[0] for x in ins_start]
        if pairs[base[0]] == ref_start or pairs[base[0]] - 1 not in ins_start_rpos:
          try: hap2_modpos[pairs[base[0]] -1].append(base[1])
          except: hap2_modpos[pairs[base[0]] -1] = [base[1]]
        #when the CpG overlaps the last pos of ins and next pos with ref base
        elif pairs[base[0]] -1 in ins_start_rpos: 
          ins_start_lengths = [x[2] for x in ins_start]
          ins_pos_relative = str(pairs[base[0]] - 1) + "+" + str(ins_start_lengths[ins_start_rpos.index(pairs[base[0]] - 1)]) # (base-1) + (len(ins))
          try: hap2_modpos[ins_pos_relative].append(base[1])
          except: hap2_modpos[ins_pos_relative] = [base[1]]
      else:
        i = 0
        prev_ref_base = None
        while prev_ref_base ==  None:
          try: prev_ref_base = pairs[base[0]-i]
          except: prev_ref_base = ref_start
          i += 1
        for i in ins_start:
            #when the start is bigger stop
            if i[0] >= prev_ref_base: 
              #start pos of the ins: ins_start+extra_bp
              try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1] -1 ) 
              except: ins_pos_relative = str(i[0]) + "+"
              
              if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
              try: hap2_modpos[ins_pos_relative].append(base[1])
              except: hap2_modpos[ins_pos_relative] = [base[1]]
              break
           
    #do this for all reads phased or unphased (hap == None)
    if pairs[base[0]] != None:
      ins_start_rpos = [x[0] for x in ins_start]
      if pairs[base[0]] == ref_start or pairs[base[0]] - 1 not in ins_start_rpos:
        try: combined_modpos[pairs[base[0]] -1].append(base[1])
        except: combined_modpos[pairs[base[0]] -1] = [base[1]]
      #when the CpG overlaps the last pos of ins and next pos with ref base
      elif pairs[base[0]] -1 in ins_start_rpos: 
        ins_start_lengths = [x[2] for x in ins_start]
        ins_pos_relative = str(pairs[base[0]] - 1) + "+" + str(ins_start_lengths[ins_start_rpos.index(pairs[base[0]] - 1)]) # (base-1) + (len(ins))
        try: combined_modpos[ins_pos_relative].append(base[1])
        except: combined_modpos[ins_pos_relative] = [base[1]]
    else:
      i = 0
      prev_ref_base = None
      while prev_ref_base ==  None:
        try: prev_ref_base = pairs[base[0]-i]
        except: prev_ref_base = ref_start
        i += 1
      for i in ins_start:
            #when the start of the next insertion is bigger stop
            if i[0] >= prev_ref_base:  
              #start pos of the ins: ins_start+extra_bp
              try: ins_pos_relative = str(i[0]) + "+" + str(base[0] - i[1] -1 )
              except: ins_pos_relative = str(i[0]) + "+"
              if ins_pos_relative[ins_pos_relative.index("+") +1:] == "0":
                 ins_pos_relative = i[0]
              try: combined_modpos[ins_pos_relative].append(base[1])
              except: combined_modpos[ins_pos_relative] = [base[1]]
              break
    
  
#count coverage, write output
def write_output(read_spans, ins_read_spans, mod_pos_dict, hap):
  
  out = open("{}_{}_{}_{}-{}.txt".format(out_prefix, hap, chrom, start, end), "w+")
  out.write("#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            "chrom", "start", "ref_pos", "read_coverage", "num_mod_calls", 
            "num_nonzero_mod_calls", "num_high_qual_mod", "avg_qual", 
            "avg_qual_percent", "percent_high_qual_calls", "mod_qual"))
  
  for mod_pos in mod_pos_dict.keys():
    qual = mod_pos_dict[mod_pos]
    if str(mod_pos).isdigit():
      pos = mod_pos
    else:
      pos = int(mod_pos[:mod_pos.index("+")])
      
    if pos < start or pos > end:
      continue
    
    #get coverage
    coverage = 0
    if str(mod_pos).isdigit():
      for span in read_spans:
        if pos in range(span[0], span[1]):
          coverage += 1
      if coverage == 0: #coverage = 0 are soft clipped 
        continue
    elif "+" in mod_pos:
      for span in ins_read_spans:
        if pos in range(span[0], span[1]):
          coverage += 1
      if coverage == 0:
        continue
    #remove zeros from quality 
    
    qual_nonzero = [i for i in qual if i != 0]
    #if qual_nonzero == []: continue #clear blank rows
    
    #mod prob >= 50%
    high_qual_mod = [i for i in qual if i >= (255+1)/2] 
    avg_qual = sum(qual)/len(qual)
    avg_qual_percent = (avg_qual + 1)/256
    percent_high_qual_calls = len(high_qual_mod)/len(qual)
    
    out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        chrom, mod_pos, pos, coverage, len(qual), len(qual_nonzero), len(high_qual_mod),  
        avg_qual, avg_qual_percent, percent_high_qual_calls, str(qual)[1:-1]))
  out.close()
      

def main(bam_file):
  bam = pysam.AlignmentFile(bam_file, "rb")
  
  for read in bam.fetch(chrom, start, end):
    combined_span.append((read.reference_start, read.reference_end))
    for base, pos_qual in read.modified_bases.items():
      if read.has_tag("HP"):
        hap = read.get_tag("HP")
        if base[1] == 0:
          count_cpg_forward(pos_qual, read.get_aligned_pairs(), hap, 
                            ins_start_site(read.get_aligned_pairs(), 
                                           read.reference_start, hap),
                            read.reference_start)
        elif base[1] == 1:
          count_cpg_reverse(pos_qual, read.get_aligned_pairs(), hap, 
                            ins_start_site(read.get_aligned_pairs(), 
                                           read.reference_start, hap),
                            read.reference_start)
        if hap == 1:
          hap1_span.append((read.reference_start, read.reference_end))
        elif hap == 2:
          hap2_span.append((read.reference_start, read.reference_end))
      else: 
        if base[1] == 0:
          count_cpg_forward(pos_qual, read.get_aligned_pairs(), None, 
                            ins_start_site(read.get_aligned_pairs(), 
                                           read.reference_start, None),
                            read.reference_start)
        elif base[1] == 1:
          count_cpg_reverse(pos_qual, read.get_aligned_pairs(), None, 
                            ins_start_site(read.get_aligned_pairs(), 
                                           read.reference_start, None),
                            read.reference_start)

  write_output(hap1_span, hap1_ins_span, hap1_modpos, "hap1")
  write_output(hap2_span, hap2_ins_span, hap2_modpos, "hap2")
  write_output(combined_span, combined_ins_span, combined_modpos, "combined")
  return(0)
      
              

main(args.bam)