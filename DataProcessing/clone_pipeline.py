import sys
import os

#input files: 1 tumor tcr filtered contig annotations; 2 blood tcr filtered contig annotations; 3 tumor projection; 4 blood projection; 5 tumor clonotypes; 6 blood clonotypes
#output files: 7 new tumor projection; 8 new blood projection

btcr_fh = open(sys.argv[2],'r')
ttcr_fh = open(sys.argv[1],'r')
bclonotypes_fh = open(sys.argv[6],'r')
tclonotypes_fh = open(sys.argv[5],'r')

blood_clono_to_barcode = {}
tumor_clono_to_barcode = {}
blood_clonotypes = set()
tumor_clonotypes = set()


#logic to handle multiple barcodes assigned to same clonotype
def make_clono_to_barcode_dicts(fh,dictionary):
    for lines in fh:
        line = lines.rstrip().split(',')
        if line[16] != "None":
            if line[10] == "True":
                if dictionary.has_key(line[16]):
                    existing_barcodes = dictionary[line[16]]
                    existing_barcodes.append(line[0])
                    dictionary[line[16]] = existing_barcodes
                else:
                    dictionary[line[16]] = [line[0]]
    return dictionary

blood_clono_to_barcode = make_clono_to_barcode_dicts(btcr_fh,blood_clono_to_barcode)
tumor_clono_to_barcode = make_clono_to_barcode_dicts(ttcr_fh,tumor_clono_to_barcode)


btcr_fh.close()
ttcr_fh.close()

class Clonotype:
    def __init__(self,chains_aa,clonotype,clono_to_barcode,frequency):
        self.barcodes = clono_to_barcode[clonotype]
        #self.frequency = frequency
        self.frequency = str(len(set(self.barcodes)))
        self.tra = []
        self.trb = []
        chains_arr = chains_aa.split(';')
        for elem in chains_arr:
            chain = elem.split(':')
            if chain[0] == "TRA":
                self.tra.append(chain[1])
            if chain[0] == "TRB":
                self.trb.append(chain[1])
    def geta(self):
        return self.tra
    def getb(self):
        return self.trb
    def getname(self):
        return '-'.join(self.tra)+'|'+'-'.join(self.trb)

def make_clonotype_sets(fh,clono_set,clono_to_barcode):
    for lines in fh:
        line = lines.rstrip().split(',')
        clonotype_id = line[0]
        chains = line[3]
        frequency = line[1]
        if clonotype_id != "clonotype_id":
            if clono_to_barcode.has_key(clonotype_id):
                clonotype = Clonotype(chains,clonotype_id,clono_to_barcode,frequency)
            clono_set.add(clonotype)
    return clono_set

blood_clonotypes = make_clonotype_sets(bclonotypes_fh,blood_clonotypes,blood_clono_to_barcode)
tumor_clonotypes = make_clonotype_sets(tclonotypes_fh,tumor_clonotypes,tumor_clono_to_barcode)

bclonotypes_fh.close()
tclonotypes_fh.close()
tp_fh = open(sys.argv[3],'r')
bp_fh = open(sys.argv[4],'r')
output_tp_fh = open(sys.argv[7],'w')
output_bp_fh = open(sys.argv[8],'w')

def find_a_chain_matches(clonotype,other_set):
    if len(clonotype.tra) > 0:
        other_matches = set()
        for tra_chain in clonotype.tra:
            matches = filter(lambda ctype: tra_chain in ctype.tra, other_set)
            for match in matches:
                other_matches.add(match)
        if len(other_matches) == 0:
            return (matches,"")
        perfect_hits = filter(lambda ctype: clonotype.tra == ctype.tra, other_set)
        if len(perfect_hits) > 0:
            return (perfect_hits,'f')
        return (matches,'')
    return (set(),'')

def find_b_chain_matches(clonotype,other_set):
    if len(clonotype.trb) > 0:
        other_matches = set()
        for trb_chain in clonotype.trb:
            matches = filter(lambda ctype: trb_chain in ctype.trb, other_set)
            for match in matches:
                other_matches.add(match)
        if len(other_matches) == 0:
            return (matches,"")
        perfect_hits = filter(lambda ctype: clonotype.trb == ctype.trb, other_set)
        if len(perfect_hits) > 0:
            return (perfect_hits,'f')
        return (matches,'')
    return (set(),'')



def make_tsne_output(input_fh,output_fh,matching_set,other_set):
    clonotype_dict = {}
    count = 0
    for lines in input_fh:
        line = lines.rstrip().split(',')
        barcode = line[0]
        clonotype_hits = filter(lambda ctype: barcode in ctype.barcodes, matching_set)
        assert len(clonotype_hits) < 2
        if len(clonotype_hits) > 0:
            clonotype = clonotype_hits[0]
            #find all a chain matches
            a_matches = find_a_chain_matches(clonotype,other_set)
            #find all b chain matches
            b_matches = find_b_chain_matches(clonotype,other_set)
            if b_matches[1] == a_matches[1] == 'f':
                state = False
                #perfect matches
                for x in b_matches[0]:
                    for y in a_matches[0]:
                        if x.geta() == y.geta() == clonotype.geta():
                            if x.getb() == y.getb() == clonotype.getb():
                                #check that a clonotype exists with both a perfect a chain and b chain match
                                assert len(set(b_matches[0]).intersection(set(a_matches[0]))) > 0
                                state = True
                                strkey = ",".join(y.geta())+'|'+",".join(y.getb())
                                if clonotype_dict.has_key(strkey):
                                    t = clonotype_dict[strkey]
                                    t.append([lines.rstrip(),"matching",clonotype.frequency,clonotype.getname()])
                                    clonotype_dict[strkey] = t
                                else:
                                    clonotype_dict[strkey] = [[lines.rstrip(),"matching",clonotype.frequency,clonotype.getname()]]
                if not state:
                    #a and b chains both present and from different clonotypes --- not counted as hit
                    aset = set(sum(map(lambda x:x.geta(), a_matches[0]),[])) #geta gives list of clonotypes, so map returns list of lists that needs to be flattened
                    bset = set(sum(map(lambda x:x.getb(), b_matches[0]),[]))
                    #check that both a and b chains have a clonotype match
                    assert len(set(clonotype.geta()).intersection(aset)) > 0 and len(set(clonotype.getb()).intersection(bset)) > 0
                    #confirm that separate a and b chain matches are assigned to different clonotypes
                    assert len(set(b_matches[0]).intersection(set(a_matches[0]))) == 0
                    output_fh.write(lines.rstrip()+",not_matching,"+clonotype.frequency+','+clonotype.getname()+'\n')
            else:
                # one matching chain, missing second chain --- counted as hit
                state1 = False
                if bool(len(b_matches[0]) > 0) ^ bool(len(a_matches[0]) > 0):
                    if len(b_matches[0]) > 0:
                        if clonotype.getb() in map(lambda x: x.trb, b_matches[0]):
                            #confirm that there are no a matches
                            assert len(map(lambda x: x.tra, a_matches[0])) == 0
                            for x in b_matches[0]:
                                if len(x.tra) == 0:
                                    #confirm that full set of chains match
                                    assert clonotype.getb() == x.getb()
                                    if len(clonotype.tra) == 0:
                                        strkey = ",".join(x.geta())+'|'+",".join(x.getb())
                                        state1 = True
                                        if clonotype_dict.has_key(strkey):
                                            t = clonotype_dict[strkey]
                                            t.append([lines.rstrip(),"beta_matching",clonotype.frequency,clonotype.getname()])
                                            clonotype_dict[strkey] = t
                                        else:
                                            clonotype_dict[strkey] = [[lines.rstrip(),"beta_matching",clonotype.frequency,clonotype.getname()]]
                    if len(a_matches[0]) > 0:
                        if clonotype.geta() in map(lambda x: x.tra, a_matches[0]):
                            #confirm that there are no b matches
                            assert len(map(lambda x: x.trb, b_matches[0])) == 0
                            for x in a_matches[0]:
                                if len(x.trb) == 0:
                                    #confirm that full set of chains match
                                    assert clonotype.geta() == x.geta()
                                    if len(clonotype.trb) == 0:
                                        strkey = ",".join(x.geta())+'|'+",".join(x.getb())
                                        state1 = True
                                        if clonotype_dict.has_key(strkey):
                                            t = clonotype_dict[strkey]
                                            t.append([lines.rstrip(),"matching",clonotype.frequency,clonotype.getname()])
                                            clonotype_dict[strkey] = t
                                        else:
                                            clonotype_dict[strkey] = [[lines.rstrip(),"matching",clonotype.frequency,clonotype.getname()]]
                    if not state1:
                        #NEED
                        assert(len(clonotype.getname()) > 1)
                        output_fh.write(lines.rstrip()+",not_matching,"+clonotype.frequency+','+clonotype.getname()+'\n')
                else:
                    #NEED
                    assert(len(clonotype.getname()) > 1)
                    output_fh.write(lines.rstrip()+",not_matching,"+clonotype.frequency+','+clonotype.getname()+'\n')
        else:
            if count == 0:
                output_fh.write(lines.rstrip()+",Matching_pre_filter,Frequency_pre_filter,TCR\n")
                count = 1
            else:
                #assert(len(clonotype.getname()) > 1)
                output_fh.write(lines.rstrip()+",notcr,1,notcr\n")
    return clonotype_dict

tumors = make_tsne_output(tp_fh,output_tp_fh,tumor_clonotypes,blood_clonotypes)
bloods = make_tsne_output(bp_fh,output_bp_fh,blood_clonotypes,tumor_clonotypes)

intersection = set(bloods.keys()).intersection(set(tumors.keys()))

def write_out_matches(d,intersection,fh):
    for key in d.keys():
        if key in intersection:
            for x in d[key]:
                fh.write(",".join(x)+'\n')
        else:
            print "hit"
            for x in d[key]:
                x[1] = "not_matching"
                #fh.write(",".join(x)+'\n')




write_out_matches(bloods,intersection,output_bp_fh)
write_out_matches(tumors,intersection,output_tp_fh)

output_tp_fh.close()
output_bp_fh.close()



tp_fh.close()
bp_fh.close()
output_tp_fh.close()
output_bp_fh.close()

